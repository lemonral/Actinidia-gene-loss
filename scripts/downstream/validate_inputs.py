#!/usr/bin/env python3
"""Validate CD-HIT and loss-master inputs before copy-number analysis."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.name.endswith((".tsv", ".tab", ".tsv.gz", ".tab.gz")) else ","
    return pd.read_csv(path, sep=sep, compression="infer", dtype=str, keep_default_na=False)


def parse_clstr(path: Path) -> pd.DataFrame:
    """Parse CD-HIT membership while preserving periods in gene identifiers."""
    current_cluster: str | None = None
    rows: list[dict[str, str]] = []
    member = re.compile(r">(.+?)\.\.\.")
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                current_cluster = line[1:].replace(" ", "_")
                continue
            match = member.search(line)
            if match and current_cluster is not None:
                rows.append({"cluster_id": current_cluster, "gene_id": match.group(1)})
    if not rows:
        raise ValueError(f"No CD-HIT members parsed from {path}")
    return pd.DataFrame(rows)


def cluster_members(args: argparse.Namespace) -> pd.DataFrame:
    path = Path(args.clusters)
    if path.name.endswith(".clstr"):
        members = parse_clstr(path)
    else:
        members = read_table(path)
        required = [args.cluster_column, args.cluster_gene_column]
        missing = [column for column in required if column not in members.columns]
        if missing:
            raise ValueError(f"Cluster membership table missing: {', '.join(missing)}")
        members = members[required].rename(columns={
            args.cluster_column: "cluster_id", args.cluster_gene_column: "gene_id"
        })
    members["cluster_id"] = members["cluster_id"].astype(str).str.strip()
    members["gene_id"] = members["gene_id"].astype(str).str.strip()
    if (members["cluster_id"] == "").any() or (members["gene_id"] == "").any():
        raise ValueError("Cluster membership contains an empty cluster or gene ID")
    duplicated = members["gene_id"].duplicated(keep=False)
    if duplicated.any():
        examples = ", ".join(members.loc[duplicated, "gene_id"].head(5))
        raise ValueError(f"A reference gene belongs to multiple clusters: {examples}")
    sizes = members.groupby("cluster_id", as_index=False).size().rename(columns={"size": "copy_number"})
    return members.merge(sizes, on="cluster_id", validate="many_to_one")


def add_row(rows: list[dict[str, object]], level: str, check: str, value: object, detail: str) -> None:
    rows.append({"level": level, "check": check, "value": value, "detail": detail})


def read_gene_id_list(path: Path | None) -> set[str]:
    if path is None:
        return set()
    values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not values:
        raise ValueError(f"{path}: exclusion list is empty")
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: exclusion list contains duplicate gene IDs")
    if any(len(value.split()) != 1 for value in values):
        raise ValueError(f"{path}: exclusion IDs must be single whitespace-free fields")
    return set(values)


def parse_copy_values(value: str) -> list[int]:
    if not value.strip():
        return []
    try:
        return sorted({int(item.strip()) for item in value.split(",") if item.strip()})
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Expected comma-separated integer copy numbers") from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clusters", required=True, help="CD-HIT .clstr or a cluster_id/gene_id membership TSV/CSV")
    parser.add_argument("--cluster-column", default="cluster_id")
    parser.add_argument("--cluster-gene-column", default="gene_id")
    parser.add_argument("--loss-table", required=True)
    parser.add_argument("--gene-column", default="reference_gene_id")
    parser.add_argument("--sample-column", default="target_haplotype")
    parser.add_argument("--ploidy-column", default="ploidy")
    parser.add_argument("--class-column", default="classification")
    parser.add_argument("--unresolved-values", default="uncertain,unassessed,unassessed_no_candidate")
    parser.add_argument("--allow-unresolved", action="store_true")
    parser.add_argument("--min-genes", type=int, default=101, help="Use 101 for the manuscript's '>100 genes' rule")
    parser.add_argument("--cdhit-identity", type=float, default=0.90)
    parser.add_argument("--expected-samples", type=int, default=None)
    parser.add_argument("--expected-copy-numbers", default="1,2,3,4,5,6,7",
                        help="Comma-separated eligible classes expected for the Figure 3/Table S16 run; use '' to disable")
    parser.add_argument("--exclude-gene-list", type=Path,
                        help="Exact gene IDs removed after original CD-HIT copy numbers are assigned")
    parser.add_argument("--manuscript-strict", action="store_true")
    parser.add_argument("--output", type=Path, help="Optional TSV QC report")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows: list[dict[str, object]] = []
    failures: list[str] = []
    try:
        excluded = read_gene_id_list(args.exclude_gene_list)
        members_all = cluster_members(args)
        members_all["excluded_from_analysis"] = members_all["gene_id"].isin(excluded)
        members = members_all.loc[~members_all["excluded_from_analysis"]].copy()
        if members.empty:
            raise ValueError("No clustered genes remain after applying --exclude-gene-list")
        losses = read_table(Path(args.loss_table))
        needed = [args.gene_column, args.sample_column, args.class_column]
        if args.ploidy_column:
            needed.append(args.ploidy_column)
        missing = [column for column in needed if column not in losses.columns]
        if missing:
            raise ValueError(f"Loss table missing required columns: {', '.join(missing)}")
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise SystemExit(f"ERROR: {exc}")

    losses = losses.rename(columns={
        args.gene_column: "gene_id",
        args.sample_column: "sample_id",
        args.class_column: "classification",
    }).copy()
    if args.ploidy_column:
        losses = losses.rename(columns={args.ploidy_column: "ploidy"})
    else:
        losses["ploidy"] = "unspecified"
    for column in ["gene_id", "sample_id", "ploidy", "classification"]:
        losses[column] = losses[column].astype(str).str.strip()

    samples_before_exclusion = set(losses["sample_id"])
    missing_loss_exclusions = excluded.difference(set(losses["gene_id"]))
    if missing_loss_exclusions:
        examples = ", ".join(sorted(missing_loss_exclusions)[:5])
        raise SystemExit(
            f"ERROR: Loss table lacks {len(missing_loss_exclusions)} requested exclusion IDs "
            f"(for example: {examples})"
        )
    if excluded:
        observed_excluded_pairs = len(
            losses.loc[losses["gene_id"].isin(excluded), ["gene_id", "sample_id"]].drop_duplicates()
        )
        expected_excluded_pairs = len(excluded) * len(samples_before_exclusion)
        add_row(rows, "OK" if observed_excluded_pairs == expected_excluded_pairs else "ERROR",
                "excluded_gene_sample_pairs", observed_excluded_pairs,
                f"Expected {expected_excluded_pairs} before exclusion.")
        if observed_excluded_pairs != expected_excluded_pairs:
            failures.append("Loss table lacks a complete exclusion-gene × sample grid")
        losses = losses.loc[~losses["gene_id"].isin(excluded)].copy()

    unresolved = {value.strip().lower() for value in args.unresolved_values.split(",") if value.strip()}
    unresolved_rows = int(losses["classification"].str.lower().isin(unresolved).sum())
    add_row(rows, "WARN" if args.allow_unresolved else "ERROR" if unresolved_rows else "OK",
            "unresolved_loss_calls", unresolved_rows,
            "Resolve uncertain/unassessed calls before a primary rate analysis; sensitivity runs must be labelled.")
    if unresolved_rows and not args.allow_unresolved:
        failures.append("Loss table contains unresolved classifications")

    add_row(rows, "OK", "original_cluster_count", int(members_all["cluster_id"].nunique()), "Clusters before gene exclusion.")
    add_row(rows, "OK", "original_clustered_reference_genes", len(members_all), "Genes before exclusion.")
    add_row(rows, "OK", "requested_exclusion_gene_ids", len(excluded), "Exact IDs requested for exclusion.")
    add_row(rows, "OK", "excluded_clustered_gene_ids", int(members_all["excluded_from_analysis"].sum()),
            "Copy numbers were assigned before these genes were removed.")
    add_row(rows, "WARN" if excluded.difference(set(members_all["gene_id"])) else "OK",
            "exclusion_ids_absent_from_cluster_map", len(excluded.difference(set(members_all["gene_id"]))),
            "These genes were already outside the CD-HIT denominator.")
    add_row(rows, "OK", "analysis_cluster_count", int(members["cluster_id"].nunique()), "Clusters retaining at least one analysis gene.")
    add_row(rows, "OK", "analysis_clustered_reference_genes", len(members), "Each must have a loss call in every sample.")
    add_row(rows, "OK", "maximum_original_copy_number", int(members["copy_number"].max()), "Original CD-HIT cluster size.")
    add_row(rows, "OK", "cdhit_identity_declared", args.cdhit_identity, "Full CD-HIT command/version must also be recorded in run metadata.")

    duplicate_rows = int(losses.duplicated(["gene_id", "sample_id"], keep=False).sum())
    add_row(rows, "ERROR" if duplicate_rows else "OK", "loss_duplicate_gene_sample_rows", duplicate_rows,
            "One final classification per reference gene × target haplotype is required.")
    if duplicate_rows:
        failures.append("Loss table has duplicate gene × sample rows")
    empty_samples = int((losses["sample_id"] == "").sum())
    if empty_samples:
        failures.append("Loss table has empty sample IDs")
    add_row(rows, "ERROR" if empty_samples else "OK", "loss_empty_sample_ids", empty_samples, "")

    samples = sorted(value for value in losses["sample_id"].unique() if value)
    add_row(rows, "OK", "sample_count", len(samples), ", ".join(samples))
    if args.expected_samples is not None and len(samples) != args.expected_samples:
        failures.append("Unexpected number of samples")
        add_row(rows, "ERROR", "expected_sample_count", len(samples), f"Expected {args.expected_samples}.")

    cluster_genes = set(members["gene_id"])
    relevant = losses.loc[losses["gene_id"].isin(cluster_genes), ["gene_id", "sample_id"]].drop_duplicates()
    expected_pairs = len(cluster_genes) * len(samples)
    observed_pairs = len(relevant)
    missing_pairs = expected_pairs - observed_pairs
    unclustered_loss_genes = len(set(losses["gene_id"]) - cluster_genes)
    add_row(rows, "ERROR" if missing_pairs else "OK", "missing_clustered_gene_sample_calls", missing_pairs,
            f"Expected {expected_pairs}; observed {observed_pairs} classified pairs.")
    add_row(rows, "WARN" if unclustered_loss_genes else "OK", "loss_genes_without_copy_number", unclustered_loss_genes,
            "A nonzero count needs an explicit exclusion/mapping decision.")
    if missing_pairs:
        failures.append("Not all clustered genes have one loss classification per sample")

    original_per_copy = (
        members_all.groupby("copy_number", as_index=False)
        .agg(original_clusters=("cluster_id", "nunique"), original_reference_genes=("gene_id", "nunique"))
    )
    analysis_per_copy = (
        members.groupby("copy_number", as_index=False)
        .agg(clusters=("cluster_id", "nunique"), reference_genes=("gene_id", "nunique"))
    )
    per_copy = original_per_copy.merge(analysis_per_copy, on="copy_number", how="left")
    for column in ["clusters", "reference_genes"]:
        per_copy[column] = per_copy[column].fillna(0).astype(int)
    per_copy = per_copy.sort_values("copy_number")
    for _, record in per_copy.iterrows():
        number = int(record["copy_number"])
        n_genes = int(record["reference_genes"])
        eligible = n_genes >= args.min_genes
        add_row(rows, "OK", f"copy_number_{number}", n_genes,
                f"original_genes={int(record['original_reference_genes'])}; analysis_clusters={int(record['clusters'])}; "
                f"{'eligible' if eligible else 'excluded'} at min_genes={args.min_genes}.")

    expected_copy_values = parse_copy_values(args.expected_copy_numbers)
    eligible_values = per_copy.loc[per_copy["reference_genes"] >= args.min_genes, "copy_number"].astype(int).tolist()
    if expected_copy_values and eligible_values != expected_copy_values:
        failures.append("Eligible copy-number classes differ from the declared Figure/Table contract")
        add_row(rows, "ERROR", "eligible_copy_numbers", ",".join(map(str, eligible_values)),
                f"Expected {','.join(map(str, expected_copy_values))}.")
    else:
        add_row(rows, "OK", "eligible_copy_numbers", ",".join(map(str, eligible_values)),
                "Classes retained by the denominator-based filter.")

    if args.manuscript_strict:
        if abs(args.cdhit_identity - 0.90) > 1e-12:
            failures.append("The manuscript contract requires CD-HIT identity 0.90")
            add_row(rows, "ERROR", "method_cdhit_identity", args.cdhit_identity, "Expected 0.90.")
        if args.min_genes != 101:
            failures.append("The manuscript contract requires more than 100 genes per class")
            add_row(rows, "ERROR", "method_min_genes", args.min_genes, "Use 101 for a strict >100 rule.")
        if args.allow_unresolved:
            failures.append("The Methods-compatible copy-number analysis cannot allow unresolved classifications")
            add_row(rows, "ERROR", "method_unresolved_calls", "allowed", "Resolve calls before strict analysis.")

    report = pd.DataFrame(rows, columns=["level", "check", "value", "detail"])
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        report.to_csv(args.output, sep="\t", index=False)
    print(report.to_string(index=False))
    if failures:
        print("\nFAILED: " + "; ".join(sorted(set(failures))), file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
