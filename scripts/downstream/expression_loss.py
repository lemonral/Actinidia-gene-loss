#!/usr/bin/env python3
"""Build a reproducible expression-bin versus gene-loss-rate table.

This is a drop-in replacement candidate for
`scripts/expression_loss/expression_loss.py`.  It keeps the original tidy
input model but makes the expression measurement, transformation, tie rule,
and coverage denominator explicit.  Render Figure 3 from its rate table with
`make_figure3_panels.py`, rather than fitting in this data-preparation step.
"""

from __future__ import annotations

import argparse
import hashlib
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if path.name.endswith((".tsv", ".tab", ".tsv.gz", ".tab.gz")) else ","
    return pd.read_csv(path, sep=separator, compression="infer", dtype=str, keep_default_na=False)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True, help="One expression value per C. scandens reference gene")
    parser.add_argument("--gene-column", default="reference_gene_id")
    parser.add_argument("--expression-column", required=True)
    parser.add_argument("--loss-table", required=True, help="Long, complete gene-loss master table")
    parser.add_argument("--loss-gene-column", default="reference_gene_id")
    parser.add_argument("--sample-column", default="target_haplotype")
    parser.add_argument("--ploidy-column", default="ploidy", help="Set to '' only when ploidy is unavailable")
    parser.add_argument("--class-column", default="classification")
    parser.add_argument("--lost-values", default="pseudogenized,deleted,decayed,degraded,loss")
    parser.add_argument("--unresolved-values", default="uncertain,unassessed,unassessed_no_candidate",
                        help="Comma-separated classifications that must not silently enter a rate denominator")
    parser.add_argument("--allow-unresolved", action="store_true",
                        help="Permit unresolved calls only for a clearly labelled sensitivity analysis")
    parser.add_argument("--bins", type=int, default=14)
    parser.add_argument("--binning-method", choices=["rank_first", "value_quantile"], default="rank_first")
    parser.add_argument("--pseudocount", type=float, default=0.1)
    parser.add_argument("--no-log-transform", action="store_true", help="Use only for a clearly labelled legacy/non-Methods run")
    parser.add_argument("--measurement", choices=["fpkm", "tpm", "raw_count", "other", "unspecified"], default="unspecified")
    parser.add_argument("--tissue", default="unspecified")
    parser.add_argument("--dataset-id", default="unspecified")
    parser.add_argument("--exclude-gene-list", type=Path,
                        help="Exact gene IDs removed from both the expression universe and loss-call denominator")
    parser.add_argument("--collapse-duplicates", action="store_true", help="Collapse duplicate loss calls by 'any loss'; document why")
    parser.add_argument("--allow-incomplete-loss-coverage", action="store_true")
    parser.add_argument("--manuscript-strict", action="store_true", help="Require FPKM + log(FPKM+0.1) + 14 rank-first bins")
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def build_expression_table(args: argparse.Namespace, excluded: set[str]) -> tuple[pd.DataFrame, dict[str, object]]:
    expression = read_table(Path(args.expression))
    require_columns(expression, [args.gene_column, args.expression_column], "Expression table")
    expression = expression[[args.gene_column, args.expression_column]].copy()
    expression.columns = ["gene_id", "expression_value"]
    expression["gene_id"] = expression["gene_id"].astype(str).str.strip()
    if (expression["gene_id"] == "").any():
        raise ValueError("Expression table contains empty gene IDs")
    if expression["gene_id"].duplicated(keep=False).any():
        examples = ", ".join(expression.loc[expression["gene_id"].duplicated(keep=False), "gene_id"].head(5))
        raise ValueError(f"Expression table has duplicate gene IDs (for example: {examples})")
    expression_ids = set(expression["gene_id"])
    missing_exclusions = excluded.difference(expression_ids)
    if missing_exclusions:
        examples = ", ".join(sorted(missing_exclusions)[:5])
        raise ValueError(
            f"Expression table lacks {len(missing_exclusions)} requested exclusion IDs (for example: {examples})"
        )
    n_expression_before_exclusion = len(expression)
    if excluded:
        expression = expression.loc[~expression["gene_id"].isin(excluded)].copy()
    if expression.empty:
        raise ValueError("No expression genes remain after applying --exclude-gene-list")
    expression["expression_value"] = pd.to_numeric(expression["expression_value"], errors="coerce")
    if expression["expression_value"].isna().any() or (expression["expression_value"] < 0).any():
        raise ValueError("Expression values must be complete, numeric, and non-negative")
    if args.pseudocount <= 0:
        raise ValueError("--pseudocount must be positive")
    if args.no_log_transform:
        expression["transformed_expression"] = expression["expression_value"]
        transform = "identity"
    else:
        expression["transformed_expression"] = np.log(expression["expression_value"] + args.pseudocount)
        transform = f"log(value+{args.pseudocount:g})"
    if len(expression) < args.bins:
        raise ValueError(f"Only {len(expression)} expression genes are available for {args.bins} bins")
    if args.binning_method == "rank_first":
        expression["expression_rank"] = expression["transformed_expression"].rank(method="first")
        expression["expression_rank_bin"] = pd.qcut(expression["expression_rank"], q=args.bins, labels=False)
    else:
        expression["expression_rank"] = expression["transformed_expression"].rank(method="average")
        expression["expression_rank_bin"] = pd.qcut(
            expression["transformed_expression"], q=args.bins, labels=False, duplicates="drop"
        )
    actual_bins = int(expression["expression_rank_bin"].nunique())
    if actual_bins != args.bins:
        raise ValueError(
            f"Binning produced {actual_bins}, not {args.bins}, groups. Use --binning-method rank_first for tied values."
        )
    expression["expression_rank_bin"] = expression["expression_rank_bin"].astype(int)
    metadata = {
        "n_expression_genes_before_exclusion": n_expression_before_exclusion,
        "n_expression_genes": int(len(expression)),
        "excluded_expression_genes": len(excluded),
        "expression_measurement": args.measurement,
        "expression_tissue": args.tissue,
        "expression_dataset_id": args.dataset_id,
        "expression_transform": transform,
        "expression_binning_method": args.binning_method,
        "expression_bins": args.bins,
        "expression_bin_min_size": int(expression["expression_rank_bin"].value_counts().min()),
        "expression_bin_max_size": int(expression["expression_rank_bin"].value_counts().max()),
    }
    return expression, metadata


def build_loss_table(args: argparse.Namespace, excluded: set[str]) -> tuple[pd.DataFrame, dict[str, object]]:
    losses = read_table(Path(args.loss_table))
    required = [args.loss_gene_column, args.sample_column, args.class_column]
    if args.ploidy_column:
        required.append(args.ploidy_column)
    require_columns(losses, required, "Loss table")
    losses = losses.rename(columns={
        args.loss_gene_column: "gene_id",
        args.sample_column: "sample_id",
        args.class_column: "classification",
    }).copy()
    if args.ploidy_column:
        losses = losses.rename(columns={args.ploidy_column: "ploidy"})
    else:
        losses["ploidy"] = "unspecified"
    for name in ["gene_id", "sample_id", "ploidy", "classification"]:
        losses[name] = losses[name].astype(str).str.strip()
    if (losses[["gene_id", "sample_id", "classification"]] == "").any().any():
        raise ValueError("Loss table contains empty gene, sample, or classification fields")
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
            raise ValueError(f"Duplicate gene × sample × ploidy loss rows: {examples}")
        losses = (
            losses.groupby(key, as_index=False)["is_lost"].max()
            .assign(classification=lambda frame: np.where(frame["is_lost"], "collapsed_lost", "collapsed_not_lost"))
        )
    return losses, {
        "lost_values": ",".join(sorted(lost_values)),
        "unresolved_values": ",".join(sorted(unresolved_values)),
        "unresolved_rows": unresolved_rows,
        "loss_rows_input": int(len(losses)),
        "excluded_loss_genes": len(excluded),
        "loss_duplicate_rows_collapsed": n_duplicate if args.collapse_duplicates else 0,
        "loss_samples": int(losses["sample_id"].nunique()),
    }


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"ERROR: Output directory already exists: {args.output_dir}")
    if args.manuscript_strict:
        if args.measurement != "fpkm":
            raise SystemExit("ERROR: --manuscript-strict requires --measurement fpkm")
        if args.no_log_transform or args.pseudocount != 0.1:
            raise SystemExit("ERROR: --manuscript-strict requires log(FPKM + 0.1)")
        if args.bins != 14 or args.binning_method != "rank_first":
            raise SystemExit("ERROR: --manuscript-strict requires 14 rank-first bins")
        if args.allow_unresolved:
            raise SystemExit("ERROR: --manuscript-strict cannot allow unresolved loss calls")
    try:
        excluded = read_gene_id_list(args.exclude_gene_list)
        expression, expression_meta = build_expression_table(args, excluded)
        losses, loss_meta = build_loss_table(args, excluded)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise SystemExit(f"ERROR: {exc}")

    ploidy_per_sample = losses.groupby("sample_id")["ploidy"].nunique()
    if (ploidy_per_sample != 1).any():
        bad = ", ".join(ploidy_per_sample.loc[ploidy_per_sample != 1].index[:5])
        raise SystemExit(f"ERROR: Each sample must have one ploidy label (for example: {bad})")
    samples = sorted(losses["sample_id"].unique())
    expected_pairs = len(expression) * len(samples)
    observed_pairs = len(losses.loc[losses["gene_id"].isin(set(expression["gene_id"])), ["gene_id", "sample_id"]].drop_duplicates())
    missing_pairs = expected_pairs - observed_pairs
    if missing_pairs and not args.allow_incomplete_loss_coverage:
        raise SystemExit(
            "ERROR: Loss master lacks "
            f"{missing_pairs} expression-gene × sample calls; use a complete master table or explicitly opt in to incomplete coverage."
        )
    merged = losses.merge(expression, on="gene_id", how="inner", validate="many_to_one")
    rates = (
        merged.groupby(["sample_id", "ploidy", "expression_rank_bin"], as_index=False)
        .agg(
            total_genes=("gene_id", "nunique"),
            lost_genes=("is_lost", "sum"),
            expression_median=("expression_value", "median"),
            transformed_expression_median=("transformed_expression", "median"),
        )
        .assign(loss_rate=lambda frame: frame["lost_genes"] / frame["total_genes"])
        .sort_values(["sample_id", "expression_rank_bin"])
    )
    expected_rate_rows = len(samples) * args.bins
    if len(rates) != expected_rate_rows:
        raise SystemExit(f"ERROR: Produced {len(rates)} rate rows; expected {expected_rate_rows} for a complete grid")

    qc = pd.DataFrame([
        ("expected_expression_gene_sample_pairs", expected_pairs),
        ("observed_expression_gene_sample_pairs", observed_pairs),
        ("missing_expression_gene_sample_pairs", missing_pairs),
        ("loss_genes_excluded_for_no_expression", len(set(losses["gene_id"]) - set(expression["gene_id"]))),
        ("rate_rows", len(rates)),
        ("expected_rate_rows", expected_rate_rows),
    ], columns=["check", "value"])
    args.output_dir.mkdir(parents=True, exist_ok=False)
    expression.to_csv(args.output_dir / "expression_gene_bins.tsv", sep="\t", index=False)
    rates.to_csv(args.output_dir / "expression_loss_rates.tsv", sep="\t", index=False)
    qc.to_csv(args.output_dir / "input_and_coverage_qc.tsv", sep="\t", index=False)
    metadata = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(shlex.quote(value) for value in sys.argv),
        "expression_sha256": sha256(Path(args.expression)),
        "loss_table_sha256": sha256(Path(args.loss_table)),
        "exclude_gene_list": str(args.exclude_gene_list.resolve()) if args.exclude_gene_list else "",
        "exclude_gene_list_sha256": sha256(args.exclude_gene_list) if args.exclude_gene_list else "",
        "allow_incomplete_loss_coverage": args.allow_incomplete_loss_coverage,
        **expression_meta,
        **loss_meta,
    }
    pd.DataFrame(sorted(metadata.items()), columns=["key", "value"]).to_csv(
        args.output_dir / "run_metadata.tsv", sep="\t", index=False
    )
    print(f"Wrote reusable expression-loss outputs to {args.output_dir}")


if __name__ == "__main__":
    main()
