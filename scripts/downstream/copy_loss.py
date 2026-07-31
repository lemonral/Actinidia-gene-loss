#!/usr/bin/env python3
"""Build a reproducible CD-HIT copy-number versus gene-loss-rate table.

This is a replacement candidate for `scripts/copy_loss/copy_loss.py`.  It
parses either CD-HIT `.clstr` output or a tidy membership table, requires an
explicit reference-gene denominator, implements the manuscript's >100-gene
filter as `--min-genes 101`, and writes Figure 3/Table S16-ready tidy data.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.name.endswith((".tsv", ".tab", ".tsv.gz", ".tab.gz")) else ","
    return pd.read_csv(path, sep=sep, compression="infer", dtype=str, keep_default_na=False)


def parse_clstr(path: Path) -> pd.DataFrame:
    """Read CD-HIT members without truncating periods in gene identifiers."""
    current_cluster: str | None = None
    records: list[dict[str, str]] = []
    pattern = re.compile(r">(.+?)\.\.\.")
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                current_cluster = line[1:].replace(" ", "_")
                continue
            match = pattern.search(line)
            if current_cluster is not None and match:
                records.append({"cluster_id": current_cluster, "gene_id": match.group(1)})
    if not records:
        raise ValueError(f"No cluster members parsed from {path}")
    return pd.DataFrame(records)


def cluster_map(path: Path, cluster_column: str, gene_column: str) -> pd.DataFrame:
    if path.name.endswith(".clstr"):
        members = parse_clstr(path)
    else:
        members = read_table(path)
        missing = [name for name in [cluster_column, gene_column] if name not in members.columns]
        if missing:
            raise ValueError(f"Cluster table missing columns: {', '.join(missing)}")
        members = members[[cluster_column, gene_column]].rename(columns={
            cluster_column: "cluster_id", gene_column: "gene_id"
        })
    members["cluster_id"] = members["cluster_id"].astype(str).str.strip()
    members["gene_id"] = members["gene_id"].astype(str).str.strip()
    if (members[["cluster_id", "gene_id"]] == "").any().any():
        raise ValueError("Cluster table contains empty cluster or gene IDs")
    duplicated = members["gene_id"].duplicated(keep=False)
    if duplicated.any():
        examples = ", ".join(members.loc[duplicated, "gene_id"].head(5))
        raise ValueError(f"A reference gene belongs to multiple clusters: {examples}")
    sizes = members.groupby("cluster_id", as_index=False).size().rename(columns={"size": "copy_number"})
    return members.merge(sizes, on="cluster_id", validate="many_to_one")


def require_columns(frame: pd.DataFrame, names: list[str], label: str) -> None:
    missing = [name for name in names if name and name not in frame.columns]
    if missing:
        raise ValueError(f"{label} missing columns: {', '.join(missing)}")


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


def parse_ints(text: str) -> list[int]:
    if not text.strip():
        return []
    try:
        return sorted({int(item.strip()) for item in text.split(",") if item.strip()})
    except ValueError as exc:
        raise ValueError("Expected comma-separated integer values") from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clusters", required=True, help="CD-HIT .clstr or two-column membership TSV/CSV")
    parser.add_argument("--loss-table", required=True, help="Complete long gene-loss master table")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--gene-column", default="reference_gene_id")
    parser.add_argument("--sample-column", default="target_haplotype")
    parser.add_argument("--ploidy-column", default="ploidy", help="Set to '' only when ploidy is unavailable")
    parser.add_argument("--class-column", default="classification")
    parser.add_argument("--cluster-column", default="cluster_id")
    parser.add_argument("--cluster-gene-column", default="gene_id")
    parser.add_argument("--lost-values", default="pseudogenized,deleted,decayed,degraded,loss")
    parser.add_argument("--unresolved-values", default="uncertain,unassessed,unassessed_no_candidate",
                        help="Comma-separated classifications that must not silently enter a rate denominator")
    parser.add_argument("--allow-unresolved", action="store_true",
                        help="Permit unresolved calls only for a clearly labelled sensitivity analysis")
    parser.add_argument("--min-genes", type=int, default=1, help="Use 101 for the manuscript's strict >100 rule")
    parser.add_argument("--cdhit-identity", type=float, default=0.90, help="Metadata/QC declaration; clustering is performed upstream")
    parser.add_argument("--expected-copy-numbers", default="", help="Optional eligible classes, e.g. 1,2,3,4,5,6,7")
    parser.add_argument("--exclude-gene-list", type=Path,
                        help="Exact gene IDs removed after original CD-HIT copy numbers are assigned")
    parser.add_argument("--collapse-duplicates", action="store_true")
    parser.add_argument("--allow-missing-copy", action="store_true")
    parser.add_argument("--allow-incomplete-loss-coverage", action="store_true")
    parser.add_argument("--manuscript-strict", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args()


def build_loss_table(args: argparse.Namespace, excluded: set[str]) -> tuple[pd.DataFrame, dict[str, object]]:
    losses = read_table(Path(args.loss_table))
    required = [args.gene_column, args.sample_column, args.class_column]
    if args.ploidy_column:
        required.append(args.ploidy_column)
    require_columns(losses, required, "Loss table")
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
    if (losses[["gene_id", "sample_id", "classification"]] == "").any().any():
        raise ValueError("Loss table contains empty gene/sample/classification fields")
    samples_before_exclusion = set(losses["sample_id"])
    missing_exclusions = excluded.difference(set(losses["gene_id"]))
    if missing_exclusions:
        examples = ", ".join(sorted(missing_exclusions)[:5])
        raise ValueError(
            f"Loss table lacks {len(missing_exclusions)} requested exclusion IDs (for example: {examples})"
        )
    if excluded:
        observed_excluded_pairs = len(
            losses.loc[losses["gene_id"].isin(excluded), ["gene_id", "sample_id"]].drop_duplicates()
        )
        expected_excluded_pairs = len(excluded) * len(samples_before_exclusion)
        if observed_excluded_pairs != expected_excluded_pairs:
            raise ValueError(
                "Loss table does not contain a complete exclusion-gene × sample grid: "
                f"observed {observed_excluded_pairs}; expected {expected_excluded_pairs}"
            )
        losses = losses.loc[~losses["gene_id"].isin(excluded)].copy()
    if losses.empty:
        raise ValueError("No loss-call rows remain after applying --exclude-gene-list")
    lost_values = {value.strip().lower() for value in args.lost_values.split(",") if value.strip()}
    if not lost_values:
        raise ValueError("--lost-values cannot be empty")
    unresolved_values = {value.strip().lower() for value in args.unresolved_values.split(",") if value.strip()}
    unresolved_rows = int(losses["classification"].str.lower().isin(unresolved_values).sum())
    if unresolved_rows and not args.allow_unresolved:
        example = ", ".join(losses.loc[losses["classification"].str.lower().isin(unresolved_values), "classification"].head(3))
        raise ValueError(
            f"Loss table contains {unresolved_rows} unresolved calls ({example}); resolve them or use --allow-unresolved for a documented sensitivity analysis"
        )
    losses["is_lost"] = losses["classification"].str.lower().isin(lost_values)
    key = ["gene_id", "sample_id", "ploidy"]
    duplicate = losses.duplicated(key, keep=False)
    n_duplicate = int(duplicate.sum())
    if n_duplicate:
        if not args.collapse_duplicates:
            examples = losses.loc[duplicate, key].head(5).to_dict("records")
            raise ValueError(f"Duplicate gene × sample × ploidy loss calls: {examples}")
        losses = (
            losses.groupby(key, as_index=False)["is_lost"].max()
            .assign(classification=lambda frame: np.where(frame["is_lost"], "collapsed_lost", "collapsed_not_lost"))
        )
    ploidy_per_sample = losses.groupby("sample_id")["ploidy"].nunique()
    if (ploidy_per_sample != 1).any():
        bad = ", ".join(ploidy_per_sample.loc[ploidy_per_sample != 1].index[:5])
        raise ValueError(f"A sample has more than one ploidy label: {bad}")
    return losses, {
        "lost_values": ",".join(sorted(lost_values)),
        "unresolved_values": ",".join(sorted(unresolved_values)),
        "unresolved_rows": unresolved_rows,
        "loss_rows": len(losses),
        "excluded_loss_genes": len(excluded),
        "loss_samples": int(losses["sample_id"].nunique()),
        "loss_duplicate_rows_collapsed": n_duplicate if args.collapse_duplicates else 0,
    }


def fit_models(rates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    usable = rates.loc[rates["passes_min_genes"]].copy()
    if len(usable) < 3 or usable["copy_number"].nunique() < 2:
        return pd.DataFrame(columns=["n_points", "slope", "intercept", "r_value", "r_squared", "p_value", "stderr"]), pd.DataFrame(columns=["copy_number", "lowess_loss_rate"])
    try:
        from scipy.stats import linregress
        fit = linregress(usable["copy_number"], usable["loss_rate"])
        ols = pd.DataFrame([{
            "n_points": len(usable), "slope": fit.slope, "intercept": fit.intercept,
            "r_value": fit.rvalue, "r_squared": fit.rvalue ** 2, "p_value": fit.pvalue, "stderr": fit.stderr,
        }])
    except ImportError:
        ols = pd.DataFrame([{"n_points": len(usable), "note": "scipy_not_available"}])
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        result = lowess(usable["loss_rate"], usable["copy_number"], frac=0.6, return_sorted=True)
        smooth = pd.DataFrame(result, columns=["copy_number", "lowess_loss_rate"])
    except ImportError:
        smooth = pd.DataFrame(columns=["copy_number", "lowess_loss_rate", "note"])
    return ols, smooth


def plot_rates(rates: pd.DataFrame, ols: pd.DataFrame, lowess: pd.DataFrame, output: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    usable = rates.loc[rates["passes_min_genes"]].copy()
    if usable.empty:
        return
    palette = {"diploid": "#1f77b4", "tetraploid": "#ff7f0e", "hexaploid": "#2ca02c"}
    fig, ax = plt.subplots(figsize=(8, 5), dpi=220)
    for ploidy, subset in usable.groupby("ploidy", sort=True):
        ax.scatter(subset["copy_number"], subset["loss_rate"], s=24, alpha=0.78,
                   color=palette.get(str(ploidy).lower()), label=ploidy)
    if not ols.empty and {"slope", "intercept"}.issubset(ols.columns):
        values = np.linspace(usable["copy_number"].min(), usable["copy_number"].max(), 200)
        ax.plot(values, ols.iloc[0]["intercept"] + ols.iloc[0]["slope"] * values, color="black", linewidth=1.2, label="OLS")
    if not lowess.empty and {"copy_number", "lowess_loss_rate"}.issubset(lowess.columns):
        ax.plot(lowess["copy_number"], lowess["lowess_loss_rate"], color="grey", linestyle="--", linewidth=1.0, label="LOWESS")
    ax.set_xlabel("Reference gene-family copy number")
    ax.set_ylabel("Gene-loss rate")
    ax.set_xticks(sorted(usable["copy_number"].unique()))
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, title="Ploidy")
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"ERROR: Output directory already exists: {args.output_dir}")
    if args.min_genes < 1:
        raise SystemExit("ERROR: --min-genes must be at least one")
    if args.manuscript_strict:
        if abs(args.cdhit_identity - 0.90) > 1e-12:
            raise SystemExit("ERROR: --manuscript-strict requires --cdhit-identity 0.90")
        if args.min_genes != 101:
            raise SystemExit("ERROR: --manuscript-strict requires --min-genes 101")
        if args.allow_unresolved:
            raise SystemExit("ERROR: --manuscript-strict cannot allow unresolved loss calls")
    try:
        excluded = read_gene_id_list(args.exclude_gene_list)
        clusters_all = cluster_map(Path(args.clusters), args.cluster_column, args.cluster_gene_column)
        # Copy number is an immutable property of the original CD-HIT cluster.
        # Exclusion happens only after cluster sizes have been assigned.
        clusters_all["excluded_from_analysis"] = clusters_all["gene_id"].isin(excluded)
        clusters = clusters_all.loc[~clusters_all["excluded_from_analysis"]].copy()
        if clusters.empty:
            raise ValueError("No clustered genes remain after applying --exclude-gene-list")
        losses, loss_meta = build_loss_table(args, excluded)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise SystemExit(f"ERROR: {exc}")

    samples = sorted(losses["sample_id"].unique())
    cluster_genes = set(clusters["gene_id"])
    observed_pairs = len(losses.loc[losses["gene_id"].isin(cluster_genes), ["gene_id", "sample_id"]].drop_duplicates())
    expected_pairs = len(cluster_genes) * len(samples)
    missing_pairs = expected_pairs - observed_pairs
    if missing_pairs and not args.allow_incomplete_loss_coverage:
        raise SystemExit(
            f"ERROR: Loss master lacks {missing_pairs} clustered-gene × sample calls; a complete denominator is required."
        )

    merged = losses.merge(clusters[["gene_id", "cluster_id", "copy_number"]], on="gene_id", how="left", validate="many_to_one")
    missing_copy = merged["copy_number"].isna()
    missing_copy_rows = int(missing_copy.sum())
    missing_copy_unique_genes = int(merged.loc[missing_copy, "gene_id"].nunique())
    if missing_copy.any() and not args.allow_missing_copy:
        examples = ", ".join(merged.loc[missing_copy, "gene_id"].head(5))
        raise SystemExit(f"ERROR: {int(missing_copy.sum())} loss rows lack a CD-HIT copy number (for example: {examples})")
    merged = merged.loc[~missing_copy].copy()
    merged["copy_number"] = merged["copy_number"].astype(int)
    rates = (
        merged.groupby(["sample_id", "ploidy", "copy_number"], as_index=False)
        .agg(total_genes=("gene_id", "nunique"), lost_genes=("is_lost", "sum"))
        .assign(loss_rate=lambda frame: frame["lost_genes"] / frame["total_genes"])
        .sort_values(["sample_id", "copy_number"])
    )
    original_class_qc = (
        clusters_all.groupby("copy_number", as_index=False)
        .agg(original_n_clusters=("cluster_id", "nunique"), original_reference_genes=("gene_id", "nunique"))
    )
    analysis_class_qc = (
        clusters.groupby("copy_number", as_index=False)
        .agg(n_clusters=("cluster_id", "nunique"), reference_genes=("gene_id", "nunique"))
    )
    class_qc = original_class_qc.merge(analysis_class_qc, on="copy_number", how="left")
    for column in ["n_clusters", "reference_genes"]:
        class_qc[column] = class_qc[column].fillna(0).astype(int)
    class_qc["excluded_reference_genes"] = (
        class_qc["original_reference_genes"] - class_qc["reference_genes"]
    )
    class_qc["passes_min_genes"] = class_qc["reference_genes"] >= args.min_genes
    class_qc = class_qc.sort_values("copy_number")
    eligible = class_qc.loc[class_qc["passes_min_genes"], "copy_number"].astype(int).tolist()
    # The Methods threshold is defined on the immutable reference copy class,
    # not on the number of resolved calls left for an individual assembly.
    # With denominator-aware filtering, per-sample totals legitimately differ
    # because uncertain calls are excluded rather than counted as retained.
    rates["passes_min_genes"] = rates["copy_number"].isin(set(eligible))
    expected_copy = parse_ints(args.expected_copy_numbers)
    if expected_copy and eligible != expected_copy:
        raise SystemExit(f"ERROR: Eligible classes are {eligible}; expected {expected_copy}")
    if not args.allow_incomplete_loss_coverage:
        expected_rate_rows = len(samples) * int((class_qc["reference_genes"] > 0).sum())
        if len(rates) != expected_rate_rows:
            raise SystemExit(f"ERROR: Produced {len(rates)} rows; expected {expected_rate_rows} for complete coverage")
        totals_per_copy = rates.groupby("copy_number")["total_genes"].nunique()
        if (totals_per_copy != 1).any():
            bad = ", ".join(map(str, totals_per_copy.loc[totals_per_copy != 1].index.tolist()))
            raise SystemExit(f"ERROR: Copy-class denominators vary across samples for: {bad}")

    ols, lowess = fit_models(rates)
    eligible_rates = rates.loc[rates["passes_min_genes"]].copy()
    qc = pd.DataFrame([
        ("expected_clustered_gene_sample_pairs", expected_pairs),
        ("observed_clustered_gene_sample_pairs", observed_pairs),
        ("missing_clustered_gene_sample_pairs", missing_pairs),
        ("loss_rows_without_copy_number", missing_copy_rows),
        ("unique_loss_genes_without_copy_number", missing_copy_unique_genes),
        ("requested_exclusion_gene_ids", len(excluded)),
        ("excluded_clustered_gene_ids", int(clusters_all["excluded_from_analysis"].sum())),
        ("exclusion_ids_absent_from_cluster_map", len(excluded.difference(set(clusters_all["gene_id"])))),
        ("eligible_copy_numbers", ",".join(map(str, eligible))),
        ("rate_rows", len(rates)),
        ("eligible_rate_rows", len(eligible_rates)),
    ], columns=["check", "value"])
    args.output_dir.mkdir(parents=True, exist_ok=False)
    clusters_all.to_csv(args.output_dir / "gene_copy_map.tsv", sep="\t", index=False)
    rates.to_csv(args.output_dir / "copy_loss_rates.tsv", sep="\t", index=False)
    eligible_rates.to_csv(args.output_dir / "copy_loss_rates_eligible.tsv", sep="\t", index=False)
    class_qc.to_csv(args.output_dir / "copy_class_qc.tsv", sep="\t", index=False)
    qc.to_csv(args.output_dir / "input_and_coverage_qc.tsv", sep="\t", index=False)
    ols.to_csv(args.output_dir / "copy_loss_ols.tsv", sep="\t", index=False)
    lowess.to_csv(args.output_dir / "copy_loss_lowess.tsv", sep="\t", index=False)
    if not args.no_plot:
        plot_rates(rates, ols, lowess, args.output_dir / "copy_loss_scatter.png")
    metadata = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(shlex.quote(value) for value in sys.argv),
        "clusters_sha256": sha256(Path(args.clusters)),
        "loss_table_sha256": sha256(Path(args.loss_table)),
        "exclude_gene_list": str(args.exclude_gene_list.resolve()) if args.exclude_gene_list else "",
        "exclude_gene_list_sha256": sha256(args.exclude_gene_list) if args.exclude_gene_list else "",
        "cdhit_identity_declared": args.cdhit_identity,
        "min_genes": args.min_genes,
        "expected_copy_numbers": args.expected_copy_numbers,
        "allow_incomplete_loss_coverage": args.allow_incomplete_loss_coverage,
        "allow_missing_copy": args.allow_missing_copy,
        **loss_meta,
    }
    pd.DataFrame(sorted(metadata.items()), columns=["key", "value"]).to_csv(
        args.output_dir / "run_metadata.tsv", sep="\t", index=False
    )
    print(f"Wrote reusable copy-number/loss outputs to {args.output_dir}")


if __name__ == "__main__":
    main()
