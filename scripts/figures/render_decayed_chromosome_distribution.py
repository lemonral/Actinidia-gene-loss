#!/usr/bin/env python3
"""Render publication figures for decayed-only chromosome-position analyses."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.figure_bundle import _register_arial_fonts, write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label


GROUP_LABELS = {
    "all_decayed": "All decayed",
    "strict_pseudogenized": "Strict pseudogenized",
    "non_strict_decayed": "Other decayed",
    "frameshift_supported": "Frameshift",
    "inframe_stop_supported": "In-frame stop",
    "frameshift_and_stop_supported": "Frameshift + stop",
}
ZONE_ORDER = (
    "Z1_terminal",
    "Z2_subterminal",
    "Z3_intermediate_outer",
    "Z4_intermediate_inner",
    "Z5_central",
)
ZONE_LABELS = {
    "Z1_terminal": "Terminal",
    "Z2_subterminal": "Subterminal",
    "Z3_intermediate_outer": "Outer-intermediate",
    "Z4_intermediate_inner": "Inner-intermediate",
    "Z5_central": "Central",
}
PLOT_COLUMNS = (
    "figure",
    "panel",
    "analysis_group",
    "assembly_unit_id",
    "chromosome_hy4a",
    "five_zone",
    "estimate",
    "ci_lower",
    "ci_upper",
    "bh_q_value",
)


class DecayedPositionFigureError(RuntimeError):
    """Raised when the decayed-position figure inputs do not close."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    if not reader.fieldnames or not rows:
        raise DecayedPositionFigureError(f"{path.name}: empty or invalid TSV")
    return rows


def f(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise DecayedPositionFigureError(f"{key}: non-finite value")
    return value


def chromosome_key(value: str) -> int:
    if not value.startswith("Chr") or not value[3:].isdigit():
        raise DecayedPositionFigureError(f"invalid chromosome label: {value}")
    return int(value[3:])


def plot_row(
    *,
    figure: str,
    panel: str,
    analysis_group: str,
    assembly_unit_id: str = "",
    chromosome_hy4a: str = "",
    five_zone: str = "",
    estimate: float,
    ci_lower: float | None = None,
    ci_upper: float | None = None,
    bh_q_value: float | None = None,
) -> dict[str, object]:
    return {
        "figure": figure,
        "panel": panel,
        "analysis_group": analysis_group,
        "assembly_unit_id": assembly_unit_id,
        "chromosome_hy4a": chromosome_hy4a,
        "five_zone": five_zone,
        "estimate": estimate,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "bh_q_value": bh_q_value,
    }


def load_inputs(analysis_dir: Path) -> dict[str, object]:
    paths = {
        "manifest": analysis_dir / "run_manifest.json",
        "unit_summary": analysis_dir / "unit_decayed_localization_summary.tsv",
        "unit_chromosome": analysis_dir / "unit_chromosome_decayed_burden.tsv",
        "chromosome_rr": analysis_dir / "chromosome_adjusted_rate_ratios.tsv",
        "chromosome_model": analysis_dir / "chromosome_model_summary.tsv",
        "pooled_zone": analysis_dir / "pooled_five_zone_decayed_burden.tsv",
        "zone_rr": analysis_dir / "five_zone_adjusted_rate_ratios.tsv",
        "zone_model": analysis_dir / "five_zone_model_summary.tsv",
        "chromosome_zone": analysis_dir
        / "pooled_chromosome_five_zone_decayed_burden.tsv",
        "zone_tests": analysis_dir / "five_zone_opportunity_tests.tsv",
    }
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_DECAYED_CHROMOSOME_POSITION_ANALYSIS":
        raise DecayedPositionFigureError("analysis manifest is not PASS")
    tables = {
        key: read_tsv(path)
        for key, path in paths.items()
        if key != "manifest"
    }
    return {"paths": paths, "manifest": manifest, **tables}


def unit_order_and_labels(
    unit_summary: list[dict[str, str]],
) -> tuple[list[str], dict[str, str]]:
    units = [row["assembly_unit_id"] for row in unit_summary]
    if len(units) != 23 or len(set(units)) != 23:
        raise DecayedPositionFigureError("expected 23 unique assembly units")
    labels = {
        row["assembly_unit_id"]: format_downstream_taxon_label(
            row["biological_species"],
            (row["haplotype_or_subgenome"],),
            abbreviate_genus=True,
            separator=" ",
        )
        for row in unit_summary
    }
    return units, labels


def style_matplotlib() -> None:
    import matplotlib as mpl

    _register_arial_fonts()
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 10.0,
            "axes.linewidth": 0.7,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 8.0,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def panel_letter(axis: object, letter: str) -> None:
    axis.text(
        -0.075,
        1.035,
        f"({letter})",
        transform=axis.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
    )


def render_between(
    data: dict[str, object],
    output_dir: Path,
    dpi: int,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import Normalize
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    style_matplotlib()
    unit_rows = data["unit_summary"]
    units, labels = unit_order_and_labels(unit_rows)
    chromosomes = [f"Chr{index:02d}" for index in range(1, 30)]
    raw_rows = [
        row
        for row in data["unit_chromosome"]
        if row["analysis_group"] == "all_decayed"
    ]
    if len(raw_rows) != 23 * 29:
        raise DecayedPositionFigureError("all-decayed unit × chromosome grid is incomplete")
    raw_lookup = {
        (row["assembly_unit_id"], row["chromosome_hy4a"]): f(
            row, "decayed_loci_per_1000_genes"
        )
        for row in raw_rows
    }
    matrix = np.asarray(
        [[raw_lookup[(unit, chromosome)] for chromosome in chromosomes] for unit in units]
    )
    rr_rows = sorted(
        [
            row
            for row in data["chromosome_rr"]
            if row["analysis_group"] == "all_decayed"
        ],
        key=lambda row: chromosome_key(row["chromosome_hy4a"]),
    )
    if [row["chromosome_hy4a"] for row in rr_rows] != chromosomes:
        raise DecayedPositionFigureError("adjusted chromosome estimates are incomplete")
    model = next(
        row
        for row in data["chromosome_model"]
        if row["analysis_group"] == "all_decayed"
    )

    figure = plt.figure(figsize=(7.2, 8.2), constrained_layout=True)
    grid = figure.add_gridspec(2, 1, height_ratios=(1.0, 2.35), hspace=0.08)
    axis_a = figure.add_subplot(grid[0, 0])
    x = np.arange(29)
    adjusted = np.asarray([f(row, "adjusted_decayed_loci_per_1000_genes") for row in rr_rows])
    lower = np.asarray([f(row, "adjusted_rate_95ci_lower_per_1000") for row in rr_rows])
    upper = np.asarray([f(row, "adjusted_rate_95ci_upper_per_1000") for row in rr_rows])
    q_values = np.asarray([f(row, "bh_q_value") for row in rr_rows])
    grand_mean = float(
        np.median(
            [
                f(row, "adjusted_decayed_loci_per_1000_genes")
                / f(row, "rate_ratio")
                for row in rr_rows
            ]
        )
    )
    colors = [
        "#B24A4A" if q < 0.05 and value > grand_mean else
        "#3D6E9E" if q < 0.05 and value < grand_mean else
        "#A0A0A0"
        for value, q in zip(adjusted, q_values)
    ]
    axis_a.errorbar(
        x,
        adjusted,
        yerr=np.vstack([adjusted - lower, upper - adjusted]),
        fmt="none",
        ecolor="#555555",
        elinewidth=0.65,
        capsize=1.5,
        zorder=1,
    )
    axis_a.scatter(x, adjusted, c=colors, s=20, edgecolor="white", linewidth=0.35, zorder=2)
    axis_a.axhline(grand_mean, color="#555555", linestyle="--", linewidth=0.8)
    axis_a.set_xticks(x, chromosomes, rotation=90)
    axis_a.set_ylabel("Adjusted decayed loci\nper 1,000 target genes")
    axis_a.grid(axis="y", color="#D8D8D8", linewidth=0.45)
    axis_a.set_axisbelow(True)
    axis_a.text(
        0.995,
        1.02,
        rf"NB LRT: $\chi^2_{{28}}$ = {f(model, 'likelihood_ratio_chi_square'):.1f}, "
        r"$P < 0.001$",
        ha="right",
        va="bottom",
        transform=axis_a.transAxes,
        fontsize=7.5,
    )
    from matplotlib.lines import Line2D

    axis_a.legend(
        handles=(
            Line2D([], [], marker="o", linestyle="none", color="#B24A4A", label="Higher"),
            Line2D([], [], marker="o", linestyle="none", color="#3D6E9E", label="Lower"),
            Line2D(
                [],
                [],
                marker="o",
                linestyle="none",
                color="#A0A0A0",
                label="Not significant",
            ),
        ),
        frameon=False,
        ncol=3,
        loc="upper left",
        handletextpad=0.35,
        columnspacing=0.9,
    )
    panel_letter(axis_a, "a")

    axis_b = figure.add_subplot(grid[1, 0])
    low, high = np.quantile(matrix, (0.02, 0.98))
    image = axis_b.imshow(
        matrix,
        aspect="auto",
        cmap="viridis",
        norm=Normalize(vmin=float(low), vmax=float(high)),
        interpolation="nearest",
    )
    axis_b.set_xticks(np.arange(29), chromosomes, rotation=90)
    axis_b.set_yticks(np.arange(23), [labels[unit] for unit in units])
    axis_b.set_xlabel("HY4A-standardized chromosome")
    colorbar = figure.colorbar(image, ax=axis_b, fraction=0.025, pad=0.012)
    colorbar.set_label("Decayed loci per 1,000 target genes", fontsize=8)
    colorbar.ax.tick_params(labelsize=7)
    axis_b.set_xticks(np.arange(-0.5, 29, 1), minor=True)
    axis_b.set_yticks(np.arange(-0.5, 23, 1), minor=True)
    axis_b.grid(which="minor", color="white", linewidth=0.25, alpha=0.45)
    axis_b.tick_params(which="minor", bottom=False, left=False)
    panel_letter(axis_b, "b")

    plot_rows = [
        plot_row(
            figure="between_chromosomes",
            panel="a",
            analysis_group="all_decayed",
            chromosome_hy4a=row["chromosome_hy4a"],
            estimate=f(row, "adjusted_decayed_loci_per_1000_genes"),
            ci_lower=f(row, "adjusted_rate_95ci_lower_per_1000"),
            ci_upper=f(row, "adjusted_rate_95ci_upper_per_1000"),
            bh_q_value=f(row, "bh_q_value"),
        )
        for row in rr_rows
    ]
    plot_rows.extend(
        plot_row(
            figure="between_chromosomes",
            panel="b",
            analysis_group="all_decayed",
            assembly_unit_id=row["assembly_unit_id"],
            chromosome_hy4a=row["chromosome_hy4a"],
            estimate=f(row, "decayed_loci_per_1000_genes"),
        )
        for row in raw_rows
    )
    caption = (
        "Decayed-only gene-loss burden among HY4A-standardized chromosomes. "
        "(a) Negative-binomial adjusted rates and 95% confidence intervals; red and "
        "blue points are respectively above and below the adjusted grand mean at "
        "Benjamini-Hochberg q < 0.05, and grey points are not significant. "
        "(b) Unadjusted burden for each of 23 independently retained assembly units. "
        "Numerators contain only article-method decayed calls with a target-genome "
        "coordinate; denominators are annotated genes on the same target chromosome. "
        "Deleted calls are excluded because they have no observed target-genome locus."
    )
    validation = {
        "status": "PASS_DECAYED_BETWEEN_CHROMOSOME_FIGURE",
        "assembly_units": 23,
        "chromosomes": 29,
        "numerator": "spatially placed article-method decayed calls only",
        "denominator": "target annotated genes on the same chromosome",
        "deleted_included": False,
        "chromosome_model_lrt_p_value": f(model, "p_value"),
    }
    paths = data["paths"]
    try:
        bundle = write_figure_bundle(
            figure=figure,
            output_dir=output_dir,
            basename="decayed_between_chromosomes",
            plot_rows=plot_rows,
            plot_columns=PLOT_COLUMNS,
            caption=caption,
            validation=validation,
            input_paths=(
                paths["manifest"],
                paths["unit_summary"],
                paths["unit_chromosome"],
                paths["chromosome_rr"],
                paths["chromosome_model"],
            ),
            dpi=dpi,
        )
    finally:
        plt.close(figure)
    return bundle.directory


def render_within(
    data: dict[str, object],
    output_dir: Path,
    dpi: int,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import Normalize

    style_matplotlib()
    group_styles = {
        "all_decayed": ("#1F4E79", "o", "-"),
        "non_strict_decayed": ("#D47A2C", "s", "--"),
        "strict_pseudogenized": ("#8F4A8B", "^", "-."),
    }
    strict_styles = {
        "frameshift_supported": ("#3D7A57", "o"),
        "inframe_stop_supported": ("#B56A2D", "s"),
        "frameshift_and_stop_supported": ("#7A4E9D", "^"),
    }
    pooled = data["pooled_zone"]
    pooled_lookup = {
        (row["analysis_group"], row["five_zone"]): row for row in pooled
    }
    zone_rr_lookup = {
        (row["analysis_group"], row["five_zone"]): row
        for row in data["zone_rr"]
    }
    for group in (*group_styles, *strict_styles):
        if any((group, zone) not in pooled_lookup for zone in ZONE_ORDER):
            raise DecayedPositionFigureError(f"incomplete pooled five-zone grid: {group}")
    all_chromosome_zone = [
        row
        for row in data["chromosome_zone"]
        if row["analysis_group"] == "all_decayed"
    ]
    if len(all_chromosome_zone) != 29 * 5:
        raise DecayedPositionFigureError("pooled chromosome × five-zone grid is incomplete")
    chromosomes = [f"Chr{index:02d}" for index in range(1, 30)]
    heat_lookup = {
        (row["chromosome_hy4a"], row["five_zone"]): f(
            row, "decayed_loci_per_1000_genes"
        )
        for row in all_chromosome_zone
    }
    heat = np.asarray(
        [[heat_lookup[(chromosome, zone)] for zone in ZONE_ORDER] for chromosome in chromosomes]
    )
    chromosome_q = {
        row["chromosome_hy4a"]: f(row, "bh_q_value_across_29_chromosomes")
        for row in data["zone_tests"]
        if row["analysis_group"] == "all_decayed"
        and row["test_scope"] == "single_chromosome"
    }
    if set(chromosome_q) != set(chromosomes):
        raise DecayedPositionFigureError("single-chromosome five-zone tests are incomplete")
    zone_model = next(
        row for row in data["zone_model"] if row["analysis_group"] == "all_decayed"
    )

    figure = plt.figure(figsize=(7.2, 8.25), constrained_layout=True)
    grid = figure.add_gridspec(
        2,
        1,
        height_ratios=(1.0, 1.38),
        hspace=0.08,
    )
    bottom_grid = grid[1, 0].subgridspec(
        1,
        3,
        width_ratios=(1.16, 0.24, 1.16),
        wspace=0.02,
    )
    x = np.arange(5)
    axis_a = figure.add_subplot(grid[0, 0])
    for group, (color, marker, linestyle) in group_styles.items():
        rows = [zone_rr_lookup[(group, zone)] for zone in ZONE_ORDER]
        values = np.asarray([f(row, "rate_ratio") for row in rows])
        lower = np.asarray([f(row, "rate_ratio_95ci_lower") for row in rows])
        upper = np.asarray([f(row, "rate_ratio_95ci_upper") for row in rows])
        axis_a.errorbar(
            x,
            values,
            yerr=np.vstack([values - lower, upper - values]),
            color=color,
            marker=marker,
            linestyle=linestyle,
            linewidth=1.25,
            markersize=4.5,
            capsize=2,
            label=GROUP_LABELS[group],
        )
    axis_a.set_xticks(x, [ZONE_LABELS[zone] for zone in ZONE_ORDER])
    axis_a.axhline(1.0, color="#777777", linestyle=":", linewidth=0.8)
    axis_a.set_ylabel("Adjusted rate ratio vs central zone")
    axis_a.grid(axis="y", color="#D8D8D8", linewidth=0.45)
    axis_a.set_axisbelow(True)
    axis_a.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.17))
    axis_a.text(
        0.995,
        0.96,
        rf"All decayed, NB LRT: $\chi^2_4$ = "
        rf"{f(zone_model, 'likelihood_ratio_chi_square'):.1f}, "
        r"$P < 0.001$",
        ha="right",
        va="top",
        transform=axis_a.transAxes,
        fontsize=7.5,
    )
    panel_letter(axis_a, "a")

    axis_b = figure.add_subplot(bottom_grid[0, 0])
    low, high = np.quantile(heat, (0.02, 0.98))
    image = axis_b.imshow(
        heat,
        aspect="auto",
        cmap="YlGnBu",
        norm=Normalize(vmin=float(low), vmax=float(high)),
        interpolation="nearest",
    )
    y_labels = [
        chromosome + (" *" if chromosome_q[chromosome] < 0.05 else "")
        for chromosome in chromosomes
    ]
    axis_b.set_xticks(x, [ZONE_LABELS[zone] for zone in ZONE_ORDER], rotation=35, ha="right")
    axis_b.set_yticks(np.arange(29), y_labels)
    axis_b.set_xlabel("Chromosome zone (end → center)")
    axis_b.set_xticks(np.arange(-0.5, 5, 1), minor=True)
    axis_b.set_yticks(np.arange(-0.5, 29, 1), minor=True)
    axis_b.grid(which="minor", color="white", linewidth=0.3, alpha=0.55)
    axis_b.tick_params(which="minor", bottom=False, left=False)
    colorbar_host = figure.add_subplot(bottom_grid[0, 1])
    colorbar_host.set_axis_off()
    colorbar_axis = inset_axes(
        colorbar_host,
        width="18%",
        height="60%",
        loc="center left",
        bbox_to_anchor=(0.0, 0.0, 1.0, 1.0),
        bbox_transform=colorbar_host.transAxes,
        borderpad=0,
    )
    colorbar = figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_label("All decayed per 1,000 target genes", fontsize=8, labelpad=2)
    colorbar.ax.tick_params(labelsize=7, pad=1)
    axis_b.text(
        0.995,
        1.01,
        r"* five-zone heterogeneity, BH $q < 0.05$",
        transform=axis_b.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
    )
    panel_letter(axis_b, "b")

    axis_c = figure.add_subplot(bottom_grid[0, 2])
    for group, (color, marker) in strict_styles.items():
        rows = [pooled_lookup[(group, zone)] for zone in ZONE_ORDER]
        values = np.asarray([f(row, "decayed_loci_per_1000_genes") for row in rows])
        lower = np.asarray([f(row, "poisson_95ci_lower_per_1000") for row in rows])
        upper = np.asarray([f(row, "poisson_95ci_upper_per_1000") for row in rows])
        axis_c.errorbar(
            x,
            values,
            yerr=np.vstack([values - lower, upper - values]),
            color=color,
            marker=marker,
            linewidth=1.15,
            markersize=4.3,
            capsize=2,
            label=GROUP_LABELS[group],
        )
    axis_c.set_xticks(x, ("Terminal", "Subterminal", "Outer-int.", "Inner-int.", "Central"))
    axis_c.tick_params(axis="x", rotation=55)
    axis_c.set_ylabel("Strict subtype loci per 1,000 genes")
    axis_c.grid(axis="y", color="#D8D8D8", linewidth=0.45)
    axis_c.set_axisbelow(True)
    axis_c.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    panel_letter(axis_c, "c")

    plot_rows: list[dict[str, object]] = []
    for group in group_styles:
        for zone in ZONE_ORDER:
            row = pooled_lookup[(group, zone)]
            rr = zone_rr_lookup[(group, zone)]
            plot_rows.append(
                plot_row(
                    figure="within_chromosomes",
                    panel="a",
                    analysis_group=group,
                    five_zone=zone,
                    estimate=f(rr, "rate_ratio"),
                    ci_lower=f(rr, "rate_ratio_95ci_lower"),
                    ci_upper=f(rr, "rate_ratio_95ci_upper"),
                    bh_q_value=f(rr, "bh_q_value"),
                )
            )
    plot_rows.extend(
        plot_row(
            figure="within_chromosomes",
            panel="b",
            analysis_group="all_decayed",
            chromosome_hy4a=row["chromosome_hy4a"],
            five_zone=row["five_zone"],
            estimate=f(row, "decayed_loci_per_1000_genes"),
            bh_q_value=chromosome_q[row["chromosome_hy4a"]],
        )
        for row in all_chromosome_zone
    )
    for group in strict_styles:
        for zone in ZONE_ORDER:
            row = pooled_lookup[(group, zone)]
            rr = zone_rr_lookup[(group, zone)]
            plot_rows.append(
                plot_row(
                    figure="within_chromosomes",
                    panel="c",
                    analysis_group=group,
                    five_zone=zone,
                    estimate=f(row, "decayed_loci_per_1000_genes"),
                    ci_lower=f(row, "poisson_95ci_lower_per_1000"),
                    ci_upper=f(row, "poisson_95ci_upper_per_1000"),
                    bh_q_value=f(rr, "bh_q_value"),
                )
            )
    caption = (
        "Decayed-only burden within chromosomes divided into five equal, "
        "orientation-independent zones from chromosome ends toward the center. "
        "(a) Negative-binomial adjusted rate ratios and 95% confidence intervals "
        "relative to the central zone for all decayed calls, its strict "
        "pseudogenized subset, and the remaining decayed calls. "
        "(b) All-decayed burden by chromosome and zone; an asterisk marks a "
        "chromosome with five-zone heterogeneity at BH q < 0.05. "
        "(c) Strictly supported frameshift, in-frame-stop, and combined mechanisms. "
        "Each denominator is the number of target annotated genes in the same zone. "
        "Strict categories are subsets of article-method decayed calls and do not "
        "replace the all-decayed primary positional analysis."
    )
    validation = {
        "status": "PASS_DECAYED_WITHIN_CHROMOSOME_FIGURE",
        "chromosomes": 29,
        "five_zones": list(ZONE_ORDER),
        "orientation_independent": True,
        "zone_direction": "chromosome_end_to_center",
        "numerator": "spatially placed article-method decayed calls only",
        "denominator": "target annotated genes in the same chromosome zone",
        "deleted_included": False,
        "all_decayed_zone_model_lrt_p_value": f(zone_model, "p_value"),
    }
    paths = data["paths"]
    try:
        bundle = write_figure_bundle(
            figure=figure,
            output_dir=output_dir,
            basename="decayed_within_chromosomes",
            plot_rows=plot_rows,
            plot_columns=PLOT_COLUMNS,
            caption=caption,
            validation=validation,
            input_paths=(
                paths["manifest"],
                paths["pooled_zone"],
                paths["zone_rr"],
                paths["zone_model"],
                paths["chromosome_zone"],
                paths["zone_tests"],
            ),
            dpi=dpi,
        )
    finally:
        plt.close(figure)
    return bundle.directory


def main() -> int:
    args = parse_args()
    try:
        data = load_inputs(args.analysis_dir)
        if args.output_root.exists() and any(args.output_root.iterdir()):
            raise DecayedPositionFigureError(
                f"refusing nonempty output root: {args.output_root}"
            )
        args.output_root.mkdir(parents=True, exist_ok=True)
        between = render_between(data, args.output_root / "between", args.dpi)
        within = render_within(data, args.output_root / "within", args.dpi)
        collection = {
            "status": "PASS_DECAYED_CHROMOSOME_DISTRIBUTION_FIGURE_COLLECTION",
            "bundles": [between.name, within.name],
            "primary_numerator": "article-method decayed",
            "deleted_included": False,
        }
        (args.output_root / "collection_manifest.json").write_text(
            json.dumps(collection, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    except (OSError, ValueError, TypeError, DecayedPositionFigureError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"figure_collection\t{args.output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
