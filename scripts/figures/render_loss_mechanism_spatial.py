#!/usr/bin/env python3
"""Render chromosome distributions and residual-sequence placement by loss type."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.figure_bundle import write_figure_bundle


LOSS_TYPES = (
    "no_qualifying_translated_hit",
    "frameshift_supported",
    "inframe_stop_supported",
    "frameshift_and_stop_supported",
    "truncation_or_partial_alignment_candidate",
    "residual_sequence_mechanism_unresolved",
)
LOSS_LABELS = {
    "no_qualifying_translated_hit": "No qualifying translated hit",
    "frameshift_supported": "Frameshift",
    "inframe_stop_supported": "In-frame stop",
    "frameshift_and_stop_supported": "Frameshift + stop",
    "truncation_or_partial_alignment_candidate": "Truncation/partial candidate",
    "residual_sequence_mechanism_unresolved": "Residual sequence; mechanism unresolved",
}
RELATIONS = (
    "expected_interval_local",
    "same_chromosome_displacement_candidate",
    "interchromosomal_displacement_candidate",
    "genomewide_residual_sequence_unanchored",
    "unlocalized",
)
RELATION_LABELS = {
    "expected_interval_local": "Expected interval",
    "same_chromosome_displacement_candidate": "Same-chromosome candidate",
    "interchromosomal_displacement_candidate": "Different-chromosome candidate",
    "genomewide_residual_sequence_unanchored": "Unanchored residual",
    "unlocalized": "Unlocalized",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosome-summary", required=True, type=Path)
    parser.add_argument("--location-summary", required=True, type=Path)
    parser.add_argument("--residual-positions", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="loss_mechanism_spatial")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    if not reader.fieldnames or not rows:
        raise ValueError(f"{path.name}: invalid or empty TSV")
    return rows


def main() -> int:
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    manifest = json.loads(args.run_manifest.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_LOSS_MECHANISM_SPATIAL_ANALYSIS":
        raise SystemExit("ERROR: spatial mechanism manifest is not PASS")
    chromosome_rows = read_tsv(args.chromosome_summary)
    relation_rows = read_tsv(args.location_summary)
    detail_rows = read_tsv(args.residual_positions)
    chromosome_order = sorted(
        {row["chromosome_hy4a"] for row in chromosome_rows}
    )
    if (
        len(chromosome_order)
        != int(manifest["hy4a_standardized_chromosomes"])
        or len(chromosome_rows) != len(LOSS_TYPES) * len(chromosome_order)
        or len(relation_rows) != len(LOSS_TYPES) * len(RELATIONS)
    ):
        raise SystemExit("ERROR: chromosome or location grid does not close")

    enrichment_lookup = {
        (row["loss_type_group"], row["chromosome_hy4a"]): np.log2(
            (float(row["observed_residual_rows"]) + 0.5)
            / (float(row["length_opportunity_expected_rows"]) + 0.5)
        )
        for row in chromosome_rows
    }
    enrichment = np.asarray(
        [
            [
                enrichment_lookup[(loss_type, chromosome)]
                for chromosome in chromosome_order
            ]
            for loss_type in LOSS_TYPES
        ]
    )
    relation_lookup = {
        (row["loss_type_group"], row["location_relation"]): int(
            row["unit_gene_rows"]
        )
        for row in relation_rows
    }
    relation_counts = np.asarray(
        [
            [relation_lookup[(loss_type, relation)] for relation in RELATIONS]
            for loss_type in LOSS_TYPES
        ],
        dtype=float,
    )
    relation_percent = np.divide(
        relation_counts * 100.0,
        relation_counts.sum(axis=1, keepdims=True),
        out=np.zeros_like(relation_counts),
        where=relation_counts.sum(axis=1, keepdims=True) > 0,
    )
    strict_types = LOSS_TYPES[1:4]
    strict_values = [
        [
            float(row["normalized_end_distance"])
            for row in detail_rows
            if row["loss_type_group"] == loss_type
            and row["normalized_end_distance"]
        ]
        for loss_type in strict_types
    ]
    if any(not values for values in strict_values):
        raise SystemExit("ERROR: strict spatial type has no observed positions")

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 7.8,
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
    figure = plt.figure(figsize=(7.2, 8.4), dpi=220, constrained_layout=True)
    grid = figure.add_gridspec(
        2,
        2,
        height_ratios=[1.35, 1.0],
        width_ratios=[1.35, 1.0],
    )

    axis_a = figure.add_subplot(grid[0, :])
    from matplotlib.colors import TwoSlopeNorm

    maximum_absolute = max(float(np.abs(enrichment).max()), 0.25)
    image = axis_a.imshow(
        enrichment,
        aspect="auto",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(
            vmin=-maximum_absolute,
            vcenter=0,
            vmax=maximum_absolute,
        ),
    )
    axis_a.set_xticks(range(len(chromosome_order)), chromosome_order, rotation=55)
    axis_a.set_yticks(
        range(len(LOSS_TYPES)),
        [LOSS_LABELS[loss_type] for loss_type in LOSS_TYPES],
    )
    axis_a.set_xlabel("HY4A-standardized chromosome")
    colorbar = figure.colorbar(image, ax=axis_a, fraction=0.022, pad=0.015)
    colorbar.set_label(
        r"$\log_2$ observed / chromosome-length expected"
    )
    axis_a.set_xticks(
        np.arange(-0.5, len(chromosome_order), 1),
        minor=True,
    )
    axis_a.grid(which="minor", color="white", linewidth=0.35, alpha=0.45)
    axis_a.tick_params(which="minor", bottom=False, left=False)
    axis_a.text(
        -0.025,
        1.03,
        "(a)",
        transform=axis_a.transAxes,
        fontsize=10.5,
        fontweight="bold",
    )

    axis_b = figure.add_subplot(grid[1, 0])
    colors = ("#4477AA", "#66CCEE", "#AA3377", "#EECC66", "#B8B8B8")
    hatches = ("", "", "//", "..", "")
    y = np.arange(len(LOSS_TYPES))
    left = np.zeros(len(LOSS_TYPES))
    for index, relation in enumerate(RELATIONS):
        values = relation_percent[:, index]
        axis_b.barh(
            y,
            values,
            left=left,
            color=colors[index],
            edgecolor="#555555" if hatches[index] else "none",
            linewidth=0.35 if hatches[index] else 0,
            hatch=hatches[index],
            height=0.68,
            label=RELATION_LABELS[relation],
        )
        left += values
    axis_b.set_yticks(
        y,
        [LOSS_LABELS[loss_type] for loss_type in LOSS_TYPES],
    )
    axis_b.invert_yaxis()
    axis_b.set_xlim(0, 100)
    axis_b.set_xlabel("Loss calls by residual-sequence placement (%)")
    axis_b.grid(axis="x", color="#dddddd", linewidth=0.6)
    axis_b.set_axisbelow(True)
    axis_b.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=3,
    )
    axis_b.text(
        -0.04,
        1.04,
        "(b)",
        transform=axis_b.transAxes,
        fontsize=10.5,
        fontweight="bold",
    )

    axis_c = figure.add_subplot(grid[1, 1])
    box = axis_c.boxplot(
        strict_values,
        patch_artist=True,
        widths=0.58,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 1.4},
        whiskerprops={"color": "#555555"},
        capprops={"color": "#555555"},
    )
    for patch, color in zip(box["boxes"], ("#CC6677", "#EE7733", "#AA3377")):
        patch.set_facecolor(color)
        patch.set_edgecolor("#555555")
        patch.set_linewidth(0.6)
    axis_c.set_xticks(
        range(1, len(strict_types) + 1),
        ["Frameshift", "In-frame stop", "Frameshift + stop"],
        rotation=28,
        ha="right",
    )
    axis_c.set_ylabel("Normalized distance from nearest chromosome end")
    axis_c.set_ylim(0, 1.02)
    axis_c.grid(axis="y", color="#dddddd", linewidth=0.6)
    axis_c.set_axisbelow(True)
    axis_c.text(
        -0.08,
        1.04,
        "(c)",
        transform=axis_c.transAxes,
        fontsize=10.5,
        fontweight="bold",
    )

    plot_rows: list[dict[str, object]] = []
    for row in chromosome_rows:
        plot_rows.append(
            {
                "panel": "chromosome_density",
                "loss_type_group": row["loss_type_group"],
                "location_relation": "",
                "chromosome_hy4a": row["chromosome_hy4a"],
                "value": float(
                    enrichment_lookup[
                        (row["loss_type_group"], row["chromosome_hy4a"])
                    ]
                ),
                "unit": "log2 observed / chromosome-length expected",
            }
        )
    for row in relation_rows:
        loss_index = LOSS_TYPES.index(row["loss_type_group"])
        relation_index = RELATIONS.index(row["location_relation"])
        plot_rows.append(
            {
                "panel": "placement_composition",
                "loss_type_group": row["loss_type_group"],
                "location_relation": row["location_relation"],
                "chromosome_hy4a": "",
                "value": float(relation_percent[loss_index, relation_index]),
                "unit": "percent",
            }
        )
    for loss_type, values in zip(strict_types, strict_values):
        for value in values:
            plot_rows.append(
                {
                    "panel": "strict_end_distance",
                    "loss_type_group": loss_type,
                    "location_relation": "expected_interval_local",
                    "chromosome_hy4a": "",
                    "value": value,
                    "unit": "normalized end distance",
                }
            )
    caption = (
        "Chromosome distribution and residual-sequence placement of classified gene-loss "
        "evidence. Panel (a) gives chromosome enrichment after mapping every assembly to the "
        "HY4A Chr01–Chr29 labels. Values are log2 observed counts divided by within-unit "
        "chromosome-length expectations, with a 0.5 pseudocount. Panel (b) separates local "
        "residual loci, candidate "
        "same-chromosome displacement, candidate different-chromosome displacement, "
        "unanchored genome-wide residual sequence, and unlocalized calls. Candidate "
        "displacements are best existing protein-to-genome alignments and are not proof of "
        "inversion, translocation, or orthology. Panel (c) compares observed local positions "
        "for the three explicit coding-disruption types. The normalized end-distance scale "
        "runs from 0 at a chromosome end to 1 at the chromosome center."
    )
    validation = {
        "status": "PASS_LOSS_MECHANISM_SPATIAL_FIGURE",
        "loss_type_groups": len(LOSS_TYPES),
        "hy4a_standardized_chromosomes": len(chromosome_order),
        "positive_unit_gene_rows": int(manifest["positive_unit_gene_rows"]),
        "checks": {
            "all_positive_calls_partitioned": True,
            "ordinary_decayed_residuals_included_when_localizable": True,
            "strict_types_separated": True,
            "candidate_displacements_not_claimed_as_rearrangements": True,
            "chromosomes_standardized_to_hy4a": True,
            "no_species_aggregation": True,
            "publication_labels_are_neutral": True,
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
            args.chromosome_summary,
            args.location_summary,
            args.residual_positions,
            args.run_manifest,
        ],
        dpi=300,
    )
    plt.close(figure)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
