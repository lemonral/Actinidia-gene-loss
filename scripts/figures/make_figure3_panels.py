#!/usr/bin/env python3
"""Create Figure 3-ready tidy tables, fits, and panels from expression/copy rate tables.

Both panels use a single global OLS line over the retained classes. The script
deliberately writes every plotted point and statistic before rendering matching
PNG and PDF files.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PALETTE = {"diploid": "#1f77b4", "tetraploid": "#ff7f0e", "hexaploid": "#2ca02c"}
PLOIDY_ALIASES = {
    "2x": "diploid",
    "4x": "tetraploid",
    "6x": "hexaploid",
    **{name: name for name in PALETTE},
}


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.name.endswith((".tsv", ".tab", ".tsv.gz", ".tab.gz")) else ","
    return pd.read_csv(path, sep=sep, compression="infer", dtype=str, keep_default_na=False)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_integer_values(text: str) -> list[int]:
    if not text.strip():
        return []
    try:
        return sorted({int(item.strip()) for item in text.split(",") if item.strip()})
    except ValueError as exc:
        raise ValueError("Expected a comma-separated list of integers") from exc


def require_and_convert(frame: pd.DataFrame, value_column: str, label: str) -> pd.DataFrame:
    required = ["sample_id", "ploidy", value_column, "total_genes", "lost_genes", "loss_rate"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing: {', '.join(missing)}")
    output = frame[required].copy()
    for column in [value_column, "total_genes", "lost_genes", "loss_rate"]:
        output[column] = pd.to_numeric(output[column], errors="coerce")
    invalid = int(output[[value_column, "total_genes", "lost_genes", "loss_rate"]].isna().any(axis=1).sum())
    if invalid:
        raise ValueError(f"{label} has {invalid} rows with non-numeric values")
    output["sample_id"] = output["sample_id"].astype(str).str.strip()
    output["ploidy"] = (
        output["ploidy"]
        .astype(str)
        .str.strip()
        .str.lower()
        .map(PLOIDY_ALIASES)
        .fillna("")
    )
    if (output["ploidy"] == "").any():
        raise ValueError(f"{label} contains unsupported ploidy values")
    if (output["total_genes"] <= 0).any() or (output["lost_genes"] < 0).any():
        raise ValueError(f"{label} has invalid lost/total gene counts")
    if ((output["lost_genes"] - output["total_genes"]) > 1e-9).any():
        raise ValueError(f"{label} has lost_genes greater than total_genes")
    recomputed = output["lost_genes"] / output["total_genes"]
    if not np.allclose(recomputed, output["loss_rate"], rtol=0, atol=1e-10):
        raise ValueError(f"{label} loss_rate is inconsistent with lost_genes / total_genes")
    if ((output["loss_rate"] < 0) | (output["loss_rate"] > 1)).any():
        raise ValueError(f"{label} has loss rates outside [0, 1]")
    if output.duplicated(["sample_id", value_column], keep=False).any():
        raise ValueError(f"{label} has duplicate sample × {value_column} rows")
    ploidy_per_sample = output.groupby("sample_id")["ploidy"].nunique()
    if (ploidy_per_sample != 1).any():
        raise ValueError(f"{label} assigns more than one ploidy to a sample")
    return output


def require_complete_grid(frame: pd.DataFrame, value_column: str, expected: list[int], label: str) -> None:
    observed = sorted(frame[value_column].astype(int).unique().tolist())
    if expected and observed != expected:
        raise ValueError(f"{label} has values {observed}, expected {expected}")
    if expected:
        expected_set = set(expected)
        for sample, subset in frame.groupby("sample_id"):
            values = set(subset[value_column].astype(int))
            if values != expected_set:
                raise ValueError(f"{label} sample {sample} has {sorted(values)}, expected {expected}")


def linear_statistics(x: pd.Series, y: pd.Series) -> dict[str, float | int | str]:
    from scipy.stats import linregress, spearmanr

    fit = linregress(x, y)
    spearman = spearmanr(x, y)
    return {
        "n_points": int(len(x)),
        "slope": float(fit.slope),
        "intercept": float(fit.intercept),
        "r_value": float(fit.rvalue),
        "r_squared": float(fit.rvalue ** 2),
        "p_value": float(fit.pvalue),
        "stderr": float(fit.stderr),
        "spearman_rho": float(spearman.statistic),
        "spearman_p_value": float(spearman.pvalue),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression-rates", required=True)
    parser.add_argument("--copy-rates", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--expected-samples",
        type=int,
        help="Optional frozen assembly-unit count; omit when QC changes the accepted cohort.",
    )
    parser.add_argument("--expected-expression-bins", type=int, default=14)
    parser.add_argument("--expected-copy-numbers", default="1,2,3,4,5,6,7")
    parser.add_argument("--expression-label", default="Expression level bin (ranked)")
    parser.add_argument("--figure-basename", default="Figure3")
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args()


def p_text(value: float) -> str:
    return "p < 0.001" if value < 0.001 else f"p = {value:.3f}"


def plot_panels(expression: pd.DataFrame, copies: pd.DataFrame,
                expression_stats: dict[str, float | int | str],
                copy_stats: dict[str, float | int | str], output_stem: Path,
                expression_label: str) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 8.0,
            "axes.linewidth": 0.8,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(
        2, 1, figsize=(7.2, 8.4), dpi=220, constrained_layout=True
    )
    ax_a, ax_b = axes
    for ploidy in ("diploid", "tetraploid", "hexaploid"):
        subset = expression.loc[expression["ploidy"].str.lower() == ploidy]
        if subset.empty:
            continue
        ax_a.scatter(subset["expression_rank_bin"], subset["loss_rate"], s=20, alpha=0.78,
                     color=PALETTE[ploidy], label=ploidy)
    x = np.linspace(expression["expression_rank_bin"].min(), expression["expression_rank_bin"].max(), 400)
    y = expression_stats["intercept"] + expression_stats["slope"] * x
    ax_a.plot(x, y, color="black", linewidth=1.1)
    ax_a.set_xlabel(expression_label)
    ax_a.set_ylabel("Decayed-gene loss rate")
    ax_a.set_xticks(sorted(expression["expression_rank_bin"].astype(int).unique()))
    ax_a.set_ylim(bottom=0)
    ax_a.legend(frameon=False, title=None, loc="upper right")
    ax_a.text(-0.09, 1.02, "(a)", transform=ax_a.transAxes, va="bottom", ha="left",
              fontsize=10.5, fontweight="bold")
    ax_a.text(
        0.57,
        0.67,
        rf"Global OLS: $R^2$={expression_stats['r_squared']:.3f}, "
        f"{p_text(float(expression_stats['p_value']))}",
        transform=ax_a.transAxes,
        fontsize=8.5,
    )

    for ploidy in ("diploid", "tetraploid", "hexaploid"):
        subset = copies.loc[copies["ploidy"].str.lower() == ploidy]
        if subset.empty:
            continue
        ax_b.scatter(subset["copy_number"], subset["loss_rate"], s=22, alpha=0.78,
                     color=PALETTE[ploidy], label=ploidy)
    x_copy = np.linspace(copies["copy_number"].min(), copies["copy_number"].max(), 200)
    ax_b.plot(x_copy, copy_stats["intercept"] + copy_stats["slope"] * x_copy, color="black", linewidth=1.1)
    ax_b.set_xlabel("Reference gene-family size")
    ax_b.set_ylabel("Decayed-gene loss rate")
    ax_b.set_xticks(sorted(copies["copy_number"].astype(int).unique()))
    ax_b.set_ylim(bottom=0)
    ax_b.text(-0.09, 1.02, "(b)", transform=ax_b.transAxes, va="bottom", ha="left",
              fontsize=10.5, fontweight="bold")
    ax_b.text(0.60, 0.60, rf"$R^2$={copy_stats['r_squared']:.3f}, {p_text(float(copy_stats['p_value']))}",
              transform=ax_b.transAxes, fontsize=8.5)
    for axis in axes:
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.55)
        axis.set_axisbelow(True)
        axis.spines[["top", "right"]].set_visible(False)
    png = output_stem.with_suffix(".png")
    pdf = output_stem.with_suffix(".pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def main() -> None:
    args = parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.figure_basename):
        raise SystemExit("ERROR: --figure-basename must be a simple file basename")
    if args.output_dir.exists():
        raise SystemExit(f"ERROR: Output directory already exists: {args.output_dir}")
    try:
        expression = require_and_convert(read_table(Path(args.expression_rates)), "expression_rank_bin", "Expression rate table")
        copies = require_and_convert(read_table(Path(args.copy_rates)), "copy_number", "Copy-loss rate table")
        for label, frame in [("Expression rate table", expression), ("Copy-loss rate table", copies)]:
            observed_ploidy = set(frame["ploidy"].str.lower())
            unknown_ploidy = sorted(observed_ploidy - set(PALETTE))
            if unknown_ploidy:
                raise ValueError(f"{label} contains unsupported ploidy values: {', '.join(unknown_ploidy)}")
        expected_expression = list(range(args.expected_expression_bins))
        expected_copies = parse_integer_values(args.expected_copy_numbers)
        if expected_copies:
            copies = copies.loc[
                copies["copy_number"].astype(int).isin(expected_copies)
            ].copy()
        require_complete_grid(expression, "expression_rank_bin", expected_expression, "Expression rate table")
        require_complete_grid(copies, "copy_number", expected_copies, "Copy-loss rate table")
        expression_samples = set(expression["sample_id"])
        copy_samples = set(copies["sample_id"])
        if expression_samples != copy_samples:
            raise ValueError("Expression-loss and copy-number/loss tables have different sample IDs")
        if args.expected_samples is not None and len(expression_samples) != args.expected_samples:
            raise ValueError(f"Observed {len(expression_samples)} samples; expected {args.expected_samples}")
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        raise SystemExit(f"ERROR: {exc}")

    expression["expression_rank_bin"] = expression["expression_rank_bin"].astype(int)
    copies["copy_number"] = copies["copy_number"].astype(int)
    try:
        expression_global = linear_statistics(expression["expression_rank_bin"], expression["loss_rate"])
        copy_global = linear_statistics(copies["copy_number"], copies["loss_rate"])
    except ImportError as exc:
        raise SystemExit(f"ERROR: Figure 3 fitting needs scipy: {exc}")

    global_rows = [
        {"panel": "a", "model": "global_ols_and_spearman", "scope": "all_assembly_units", **expression_global},
        {"panel": "b", "model": "global_ols_and_spearman", "scope": "all_assembly_units", **copy_global},
    ]
    per_sample_rows = []
    for sample, subset in expression.groupby("sample_id", sort=True):
        per_sample_rows.append({"sample_id": sample, "ploidy": subset["ploidy"].iloc[0],
                                **linear_statistics(subset["expression_rank_bin"], subset["loss_rate"])})
    fit_table = pd.DataFrame(global_rows)

    args.output_dir.mkdir(parents=True, exist_ok=False)
    expression_points = expression.assign(panel="a", x_value=expression["expression_rank_bin"])
    copy_points = copies.assign(panel="b", x_value=copies["copy_number"])
    expression_points.to_csv(args.output_dir / "figure3a_points.tsv", sep="\t", index=False)
    copy_points.to_csv(args.output_dir / "figure3b_points.tsv", sep="\t", index=False)
    fit_table.to_csv(args.output_dir / "figure3_fit_statistics.tsv", sep="\t", index=False)
    pd.DataFrame(per_sample_rows).to_csv(args.output_dir / "figure3a_per_sample_statistics.tsv", sep="\t", index=False)
    pd.DataFrame([
        ("timestamp_utc", datetime.now(timezone.utc).isoformat()),
        ("expression_rates_basename", Path(args.expression_rates).name),
        ("expression_rates_sha256", sha256(Path(args.expression_rates))),
        ("copy_rates_basename", Path(args.copy_rates).name),
        ("copy_rates_sha256", sha256(Path(args.copy_rates))),
        ("expected_samples", args.expected_samples),
        ("expected_expression_bins", args.expected_expression_bins),
        ("expected_copy_numbers", args.expected_copy_numbers),
        ("expression_fit", "single_global_ols"),
        ("figure_basename", args.figure_basename),
    ], columns=["key", "value"]).to_csv(args.output_dir / "run_metadata.tsv", sep="\t", index=False)
    if not args.no_plot:
        png, pdf = plot_panels(
            expression,
            copies,
            expression_global,
            copy_global,
            args.output_dir / args.figure_basename,
            args.expression_label,
        )
        pd.DataFrame(
            [
                {"format": "png", "path": png.name, "bytes": png.stat().st_size, "sha256": sha256(png)},
                {"format": "pdf", "path": pdf.name, "bytes": pdf.stat().st_size, "sha256": sha256(pdf)},
            ]
        ).to_csv(args.output_dir / "figure_manifest.tsv", sep="\t", index=False)
    print(f"Wrote Figure 3-ready tables and statistics to {args.output_dir}")


if __name__ == "__main__":
    main()
