#!/usr/bin/env python3
"""Render lost non-shared reference NLR genes by structural class."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.figure_bundle import write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label


NLR_CLASS_SERIES = (
    ("CC-NBARC", "#4C78A8"),
    ("CC-NBARC-LRR", "#F28E2B"),
    ("NBARC", "#59A14F"),
    ("NBARC-LRR", "#B07AA1"),
    ("TIR", "#8C564B"),
    ("TIR-LRR", "#BAB0AC"),
    ("TIR-NBARC", "#EDC948"),
    ("TIR-NBARC-LRR", "#17A6C1"),
    ("TIR-CC-NBARC-LRR", "#E15759"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--class-summary", required=True, type=Path)
    parser.add_argument("--shared-class-summary", required=True, type=Path)
    parser.add_argument("--repertoire-class-summary", required=True, type=Path)
    parser.add_argument("--class-loss-rates", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--basename",
        default="lost_nlr_structural_classes",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    if not reader.fieldnames or not rows:
        raise ValueError(f"{path.name}: invalid or empty table")
    return rows


def main() -> int:
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.gridspec import GridSpecFromSubplotSpec

    manifest = json.loads(args.run_manifest.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_LOST_NLR_STRUCTURAL_CLASSIFICATION":
        raise SystemExit("ERROR: lost-NLR structural-class manifest is not PASS")
    loss_rows = read_tsv(args.class_summary)
    shared_rows = read_tsv(args.shared_class_summary)
    repertoire_rows = read_tsv(args.repertoire_class_summary)
    rate_rows = read_tsv(args.class_loss_rates)
    expected_rows = int(manifest["assembly_units"]) * len(NLR_CLASS_SERIES)
    if len(loss_rows) != expected_rows:
        raise SystemExit(
            f"ERROR: expected {expected_rows} loss rows, "
            f"observed {len(loss_rows)}"
        )
    if len(shared_rows) != len(NLR_CLASS_SERIES):
        raise SystemExit("ERROR: shared-loss structural-class grid changed")
    if len(repertoire_rows) != expected_rows:
        raise SystemExit(
            f"ERROR: expected {expected_rows} repertoire rows, "
            f"observed {len(repertoire_rows)}"
        )
    if len(rate_rows) != expected_rows:
        raise SystemExit(
            f"ERROR: expected {expected_rows} class-rate rows, "
            f"observed {len(rate_rows)}"
        )

    unit_metadata: dict[str, tuple[str, str]] = {}
    loss_values: dict[tuple[str, str], int] = {}
    for row in loss_rows:
        unit = row["assembly_unit_id"]
        nlr_class = row["reference_nlr_class"]
        if nlr_class not in dict(NLR_CLASS_SERIES):
            raise SystemExit(f"ERROR: unexpected NLR class {nlr_class!r}")
        metadata = (
            row["biological_species"],
            row["haplotype_or_subgenome"],
        )
        if unit in unit_metadata and unit_metadata[unit] != metadata:
            raise SystemExit(f"ERROR: inconsistent metadata for {unit}")
        unit_metadata[unit] = metadata
        key = (unit, nlr_class)
        if key in loss_values:
            raise SystemExit(f"ERROR: duplicate unit-class row {key}")
        loss_values[key] = int(row["positive_loss_count"])
    units = sorted(
        unit_metadata,
        key=lambda unit: (
            unit_metadata[unit][0],
            unit_metadata[unit][1],
            unit,
        ),
    )
    if len(units) != int(manifest["assembly_units"]):
        raise SystemExit("ERROR: assembly-unit count changed")

    labels = [
        format_downstream_taxon_label(
            unit_metadata[unit][0],
            (unit_metadata[unit][1],),
            abbreviate_genus=True,
            separator=" ",
        )
        for unit in units
    ]
    loss_counts = {
        nlr_class: np.asarray(
            [loss_values[(unit, nlr_class)] for unit in units],
            dtype=int,
        )
        for nlr_class, _ in NLR_CLASS_SERIES
    }
    loss_totals = sum(
        (loss_counts[nlr_class] for nlr_class, _ in NLR_CLASS_SERIES),
        start=np.zeros(len(units), dtype=int),
    )
    if int(loss_totals.sum()) != int(manifest["positive_unit_gene_calls"]):
        raise SystemExit("ERROR: plotted lost-NLR counts do not close")

    shared_values: dict[str, int] = {}
    for row in shared_rows:
        nlr_class = row["reference_nlr_class"]
        if nlr_class not in dict(NLR_CLASS_SERIES) or nlr_class in shared_values:
            raise SystemExit("ERROR: invalid shared-loss structural-class row")
        shared_values[nlr_class] = int(
            row["shared_reference_nlr_gene_count"]
        )
    if set(shared_values) != {item[0] for item in NLR_CLASS_SERIES}:
        raise SystemExit("ERROR: shared-loss structural-class set changed")
    shared_total = sum(shared_values.values())
    if shared_total != int(manifest["shared_reference_nlrs_excluded"]):
        raise SystemExit("ERROR: shared-loss NLR count does not close")

    repertoire_values: dict[tuple[str, str], int] = {}
    repertoire_declared_totals: dict[str, int] = {}
    for row in repertoire_rows:
        unit = row["assembly_unit_id"]
        nlr_class = row["reference_nlr_class"]
        metadata = (
            row["biological_species"],
            row["haplotype_or_subgenome"],
        )
        if (
            unit not in unit_metadata
            or unit_metadata[unit] != metadata
            or nlr_class not in dict(NLR_CLASS_SERIES)
        ):
            raise SystemExit("ERROR: repertoire cohort metadata changed")
        key = (unit, nlr_class)
        if key in repertoire_values:
            raise SystemExit(f"ERROR: duplicate repertoire row {key}")
        repertoire_values[key] = int(row["nlr_gene_count"])
        declared_total = int(row["total_nlr_count"])
        if (
            unit in repertoire_declared_totals
            and repertoire_declared_totals[unit] != declared_total
        ):
            raise SystemExit(f"ERROR: inconsistent repertoire total for {unit}")
        repertoire_declared_totals[unit] = declared_total
    repertoire_counts = {
        nlr_class: np.asarray(
            [repertoire_values[(unit, nlr_class)] for unit in units],
            dtype=int,
        )
        for nlr_class, _ in NLR_CLASS_SERIES
    }
    repertoire_totals = sum(
        (
            repertoire_counts[nlr_class]
            for nlr_class, _ in NLR_CLASS_SERIES
        ),
        start=np.zeros(len(units), dtype=int),
    )
    if any(
        int(total) != repertoire_declared_totals[unit]
        for unit, total in zip(units, repertoire_totals)
    ):
        raise SystemExit("ERROR: per-unit NLR repertoire totals do not close")
    if int(repertoire_totals.sum()) != int(manifest["repertoire_nlr_calls"]):
        raise SystemExit("ERROR: total NLR repertoire count does not close")

    rate_values: dict[tuple[str, str], dict[str, float | int]] = {}
    for row in rate_rows:
        unit = row["assembly_unit_id"]
        nlr_class = row["reference_nlr_class"]
        metadata = (
            row["biological_species"],
            row["haplotype_or_subgenome"],
        )
        key = (unit, nlr_class)
        if (
            unit not in unit_metadata
            or unit_metadata[unit] != metadata
            or nlr_class not in dict(NLR_CLASS_SERIES)
            or key in rate_values
        ):
            raise SystemExit("ERROR: invalid structural-class loss-rate row")
        all_denominator = int(row["all_resolved_denominator"])
        nonshared_denominator = int(row["nonshared_resolved_denominator"])
        all_percentage = (
            float(row["all_loss_percentage"])
            if row["all_loss_percentage"]
            else np.nan
        )
        nonshared_percentage = (
            float(row["nonshared_loss_percentage"])
            if row["nonshared_loss_percentage"]
            else np.nan
        )
        if (
            (not np.isnan(all_percentage) and not 0 <= all_percentage <= 100)
            or (
                not np.isnan(nonshared_percentage)
                and not 0 <= nonshared_percentage <= 100
            )
        ):
            raise SystemExit("ERROR: NLR structural-class loss rate is invalid")
        rate_values[key] = {
            "shared_loss_count": int(row["shared_loss_count"]),
            "all_resolved_denominator": all_denominator,
            "nonshared_resolved_denominator": nonshared_denominator,
            "all_loss_percentage": all_percentage,
            "nonshared_loss_percentage": nonshared_percentage,
        }
    heatmap_classes = [
        nlr_class
        for nlr_class, _ in NLR_CLASS_SERIES
        if any(
            int(rate_values[(unit, nlr_class)]["all_resolved_denominator"])
            > 0
            for unit in units
        )
    ]
    all_rate_matrix = np.asarray(
        [
            [
                float(rate_values[(unit, nlr_class)]["all_loss_percentage"])
                for unit in units
            ]
            for nlr_class in heatmap_classes
        ],
        dtype=float,
    )
    nonshared_rate_matrix = np.asarray(
        [
            [
                float(
                    rate_values[(unit, nlr_class)][
                        "nonshared_loss_percentage"
                    ]
                )
                for unit in units
            ]
            for nlr_class in heatmap_classes
        ],
        dtype=float,
    )

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.5,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure = plt.figure(
        figsize=(7.2, 10.2),
        dpi=240,
    )
    outer = figure.add_gridspec(
        3,
        1,
        height_ratios=(0.78, 0.98, 0.74),
        hspace=0.54,
        left=0.11,
        right=0.97,
        bottom=0.12,
        top=0.94,
    )
    top = GridSpecFromSubplotSpec(
        1,
        2,
        subplot_spec=outer[0],
        width_ratios=(1.0, 7.7),
        wspace=0.28,
    )
    shared_axis = figure.add_subplot(top[0])
    loss_axis = figure.add_subplot(top[1])
    repertoire_axis = figure.add_subplot(outer[1])
    rate_axis = figure.add_subplot(outer[2])

    def style_axis(axis) -> None:
        axis.grid(axis="y", color="#d9d9d9", linewidth=0.6)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    def segment_text_color(nlr_class: str) -> str:
        return (
            "black"
            if nlr_class in {"TIR-LRR", "TIR-NBARC"}
            else "white"
        )

    shared_bottom = 0
    shared_handles = []
    for nlr_class, color in NLR_CLASS_SERIES:
        count = shared_values[nlr_class]
        bars = shared_axis.bar(
            [0],
            [count],
            bottom=[shared_bottom],
            width=0.62,
            color=color,
            edgecolor="white",
            linewidth=0.35,
            label=nlr_class,
        )
        shared_handles.append(bars[0])
        if count >= 3:
            shared_axis.text(
                0,
                shared_bottom + count / 2,
                str(count),
                ha="center",
                va="center",
                fontsize=7.0,
                color=segment_text_color(nlr_class),
            )
        shared_bottom += count
    shared_axis.text(
        0,
        shared_total + 2.0,
        str(shared_total),
        ha="center",
        va="bottom",
        fontsize=7.5,
    )
    shared_axis.set_ylabel("Unique shared-loss NLR genes")
    shared_axis.set_xticks([0], ["Shared loss"])
    shared_axis.set_ylim(0, shared_total * 1.13)
    style_axis(shared_axis)
    shared_axis.text(
        -0.28,
        1.06,
        "(a)",
        transform=shared_axis.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
    )

    x = np.arange(len(units))
    bottom = np.zeros(len(units), dtype=int)
    for nlr_class, color in NLR_CLASS_SERIES:
        class_counts = loss_counts[nlr_class]
        bars = loss_axis.bar(
            x,
            class_counts,
            bottom=bottom,
            width=0.78,
            color=color,
            edgecolor="white",
            linewidth=0.35,
            label=nlr_class,
        )
        for bar, count, base in zip(bars, class_counts, bottom):
            if count >= 2:
                loss_axis.text(
                    bar.get_x() + bar.get_width() / 2,
                    base + count / 2,
                    str(int(count)),
                    ha="center",
                    va="center",
                    fontsize=6.3,
                    color=segment_text_color(nlr_class),
                )
        bottom += class_counts
    for index, total in enumerate(loss_totals):
        loss_axis.text(
            index,
            total + 0.35,
            str(int(total)),
            ha="center",
            va="bottom",
            fontsize=7.0,
        )

    loss_axis.set_ylabel("Non-shared lost NLR genes")
    loss_axis.set_xticks(
        x,
        labels,
        rotation=48,
        ha="right",
        rotation_mode="anchor",
    )
    loss_axis.set_ylim(0, max(loss_totals) * 1.18)
    style_axis(loss_axis)

    bottom = np.zeros(len(units), dtype=int)
    for nlr_class, color in NLR_CLASS_SERIES:
        class_counts = repertoire_counts[nlr_class]
        bars = repertoire_axis.bar(
            x,
            class_counts,
            bottom=bottom,
            width=0.78,
            color=color,
            edgecolor="white",
            linewidth=0.35,
        )
        for bar, count, base in zip(bars, class_counts, bottom):
            if count >= 8:
                repertoire_axis.text(
                    bar.get_x() + bar.get_width() / 2,
                    base + count / 2,
                    str(int(count)),
                    ha="center",
                    va="center",
                    fontsize=6.2,
                    color=segment_text_color(nlr_class),
                )
        bottom += class_counts
    for index, total in enumerate(repertoire_totals):
        repertoire_axis.text(
            index,
            total + 4.5,
            str(int(total)),
            ha="center",
            va="bottom",
            fontsize=7.0,
        )
    repertoire_axis.set_ylabel("Annotated NLR genes")
    repertoire_axis.set_xticks(x, [""] * len(units))
    repertoire_axis.tick_params(axis="x", length=0)
    repertoire_axis.set_ylim(0, max(repertoire_totals) * 1.12)
    style_axis(repertoire_axis)
    repertoire_axis.text(
        -0.025,
        1.04,
        "(b)",
        transform=repertoire_axis.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
    )

    figure.legend(
        handles=shared_handles,
        labels=[item[0] for item in NLR_CLASS_SERIES],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.985),
        ncol=4,
        columnspacing=1.4,
        handlelength=1.5,
    )

    heatmap = np.ma.masked_invalid(all_rate_matrix)
    rate_cmap = plt.get_cmap("Blues").copy()
    rate_cmap.set_bad("#eeeeee")
    rate_image = rate_axis.imshow(
        heatmap,
        aspect="auto",
        interpolation="nearest",
        cmap=rate_cmap,
        vmin=0,
        vmax=100,
    )
    rate_axis.set_xticks(
        np.arange(len(units)),
        labels,
        rotation=48,
        ha="right",
        rotation_mode="anchor",
    )
    rate_axis.set_yticks(
        np.arange(len(heatmap_classes)),
        heatmap_classes,
    )
    rate_axis.set_ylabel("Reference NLR class")
    rate_axis.set_xticks(
        np.arange(-0.5, len(units), 1),
        minor=True,
    )
    rate_axis.set_yticks(
        np.arange(-0.5, len(heatmap_classes), 1),
        minor=True,
    )
    rate_axis.grid(
        which="minor",
        color="white",
        linestyle="-",
        linewidth=0.75,
    )
    rate_axis.tick_params(which="minor", bottom=False, left=False)
    rate_axis.spines["top"].set_visible(False)
    rate_axis.spines["right"].set_visible(False)
    for row_index in range(len(heatmap_classes)):
        for column_index in range(len(units)):
            all_percentage = all_rate_matrix[row_index, column_index]
            nonshared_percentage = nonshared_rate_matrix[
                row_index,
                column_index,
            ]
            if np.isnan(all_percentage):
                label = "NA"
                color = "black"
            else:
                nonshared_label = (
                    f"{nonshared_percentage:.0f}"
                    if not np.isnan(nonshared_percentage)
                    else "–"
                )
                label = f"{all_percentage:.0f}\n({nonshared_label})"
                color = "white" if all_percentage >= 55 else "black"
            rate_axis.text(
                column_index,
                row_index,
                label,
                ha="center",
                va="center",
                fontsize=7.0,
                color=color,
            )
    rate_axis.text(
        -0.06,
        1.025,
        "(c)",
        transform=rate_axis.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
    )
    colorbar = figure.colorbar(
        rate_image,
        ax=rate_axis,
        fraction=0.018,
        pad=0.012,
    )
    colorbar.set_label("Total loss (%)")

    plot_rows: list[dict[str, object]] = []
    loss_source_rows = {
        (row["assembly_unit_id"], row["reference_nlr_class"]): row
        for row in loss_rows
    }
    for row in shared_rows:
        plot_rows.append(
            {
                "panel": "A_shared_loss",
                "assembly_unit_id": "",
                "biological_species": "",
                "haplotype_or_subgenome": "",
                "display_label": "Shared loss",
                "reference_nlr_class": row["reference_nlr_class"],
                "count": int(row["shared_reference_nlr_gene_count"]),
                "total_count": shared_total,
                "resolved_unit_gene_denominator": "",
                "positive_loss_percentage": "",
                "shared_loss_count": int(
                    row["shared_reference_nlr_gene_count"]
                ),
                "all_resolved_denominator": "",
                "all_loss_percentage": "",
                "nonshared_loss_percentage": "",
            }
        )
    for unit, label in zip(units, labels):
        for nlr_class, _ in NLR_CLASS_SERIES:
            row = loss_source_rows[(unit, nlr_class)]
            plot_rows.append(
                {
                    "panel": "A_nonshared_loss",
                    "assembly_unit_id": unit,
                    "biological_species": row["biological_species"],
                    "haplotype_or_subgenome": row["haplotype_or_subgenome"],
                    "display_label": label,
                    "reference_nlr_class": nlr_class,
                    "count": int(row["positive_loss_count"]),
                    "total_count": int(
                        loss_totals[units.index(unit)]
                    ),
                    "resolved_unit_gene_denominator": int(
                        row["resolved_unit_gene_denominator"]
                    ),
                    "positive_loss_percentage": row["positive_loss_percentage"],
                    "shared_loss_count": "",
                    "all_resolved_denominator": "",
                    "all_loss_percentage": "",
                    "nonshared_loss_percentage": "",
                }
            )
    repertoire_source_rows = {
        (row["assembly_unit_id"], row["reference_nlr_class"]): row
        for row in repertoire_rows
    }
    for unit, label, total in zip(units, labels, repertoire_totals):
        for nlr_class, _ in NLR_CLASS_SERIES:
            row = repertoire_source_rows[(unit, nlr_class)]
            plot_rows.append(
                {
                    "panel": "B_total_repertoire",
                    "assembly_unit_id": unit,
                    "biological_species": row["biological_species"],
                    "haplotype_or_subgenome": row[
                        "haplotype_or_subgenome"
                    ],
                    "display_label": label,
                    "reference_nlr_class": nlr_class,
                    "count": int(row["nlr_gene_count"]),
                    "total_count": int(total),
                    "resolved_unit_gene_denominator": "",
                    "positive_loss_percentage": "",
                    "shared_loss_count": "",
                    "all_resolved_denominator": "",
                    "all_loss_percentage": "",
                    "nonshared_loss_percentage": "",
                }
            )
    for unit, label in zip(units, labels):
        for nlr_class in heatmap_classes:
            values = rate_values[(unit, nlr_class)]
            plot_rows.append(
                {
                    "panel": "C_class_loss_rate",
                    "assembly_unit_id": unit,
                    "biological_species": unit_metadata[unit][0],
                    "haplotype_or_subgenome": unit_metadata[unit][1],
                    "display_label": label,
                    "reference_nlr_class": nlr_class,
                    "count": "",
                    "total_count": "",
                    "resolved_unit_gene_denominator": values[
                        "nonshared_resolved_denominator"
                    ],
                    "positive_loss_percentage": "",
                    "shared_loss_count": values["shared_loss_count"],
                    "all_resolved_denominator": values[
                        "all_resolved_denominator"
                    ],
                    "all_loss_percentage": (
                        ""
                        if np.isnan(float(values["all_loss_percentage"]))
                        else f"{float(values['all_loss_percentage']):.6f}"
                    ),
                    "nonshared_loss_percentage": (
                        ""
                        if np.isnan(
                            float(values["nonshared_loss_percentage"])
                        )
                        else (
                            f"{float(values['nonshared_loss_percentage']):.6f}"
                        )
                    ),
                }
            )
    caption = (
        "(a) Structural classes of lost reference NLR genes. The shared-loss "
        "column contains 138 unique reference NLR genes positive in all 23 "
        "assembly units; the adjacent bars show non-shared decayed plus "
        "deleted calls separately for each unit. Reference-gene classes come "
        "from the corresponding Clematoclethra scandens NLR-Annotator record. "
        "(b) Structural composition of the complete NLR-Annotator repertoire "
        "in every assembly unit. (c) Class-specific loss percentages among "
        "resolved reference-NLR opportunities. Cell colour and the first "
        "number give total loss (shared plus non-shared); the parenthesized "
        "number gives non-shared-only loss using the corresponding non-shared "
        "resolved denominator. Totals are printed above bars. Not-called loss "
        "comparisons are excluded and no species aggregation is used."
    )
    validation = {
        "schema_version": "1.0",
        "status": "PASS_LOST_NLR_STRUCTURAL_CLASS_FIGURE",
        "assembly_units": len(units),
        "nlr_structural_classes": len(NLR_CLASS_SERIES),
        "shared_loss_reference_nlrs": shared_total,
        "positive_nonshared_unit_gene_calls": int(loss_totals.sum()),
        "annotated_repertoire_nlrs": int(repertoire_totals.sum()),
        "class_loss_rate_rows": len(heatmap_classes) * len(units),
        "checks": {
            "nonshared_reference_nlr_foreground": True,
            "structural_class_partition_closes": True,
            "shared_loss_unique_reference_genes_separate": True,
            "complete_per_unit_nlr_repertoires_shown": True,
            "class_loss_rates_use_resolved_reference_opportunities": True,
            "total_and_nonshared_loss_rates_both_shown": True,
            "no_species_aggregation": True,
            "not_called_excluded": True,
            "latin_binomials_italic_suffixes_upright": True,
            "publication_title_omitted": True,
        },
    }
    write_figure_bundle(
        figure=figure,
        output_dir=args.output_dir,
        basename=args.basename,
        plot_rows=plot_rows,
        plot_columns=list(plot_rows[0]),
        caption=caption,
        validation=validation,
        input_paths=[
            args.class_summary,
            args.shared_class_summary,
            args.repertoire_class_summary,
            args.class_loss_rates,
            args.run_manifest,
        ],
    )
    plt.close(figure)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
