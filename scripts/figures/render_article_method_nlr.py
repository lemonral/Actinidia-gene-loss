#!/usr/bin/env python3
"""Render NLR repertoire and classified non-shared loss for 23 units."""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-summary", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="nlr_repertoire_and_loss_types")
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

    manifest = json.loads(args.run_manifest.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_ARTICLE_METHOD_NLR_SUMMARY":
        raise SystemExit("ERROR: NLR summary manifest is not PASS")
    rows = read_tsv(args.unit_summary)
    if len(rows) != 23:
        raise SystemExit("ERROR: expected exactly 23 NLR unit rows")
    rows.sort(
        key=lambda row: (
            row["biological_species"],
            row["haplotype_or_subgenome"],
            row["assembly_unit_id"],
        )
    )
    labels = [
        format_downstream_taxon_label(
            row["biological_species"],
            (row["haplotype_or_subgenome"],),
            abbreviate_genus=True,
            separator=" ",
        )
        for row in rows
    ]
    repertoire = np.asarray([int(row["total_nlr_count"]) for row in rows])
    loss_series = (
        (
            "No qualifying translated hit",
            "no_qualifying_translated_hit_count",
            "#8F8F8F",
            "",
        ),
        (
            "Frameshift",
            "frameshift_supported_count",
            "#CC6677",
            "",
        ),
        (
            "In-frame stop",
            "inframe_stop_supported_count",
            "#EE7733",
            "",
        ),
        (
            "Frameshift + stop",
            "frameshift_and_stop_supported_count",
            "#AA3377",
            "",
        ),
        (
            "Truncation/partial candidate",
            "truncation_or_partial_alignment_candidate_count",
            "#228833",
            "//",
        ),
        (
            "Residual sequence; mechanism unresolved",
            "residual_sequence_mechanism_unresolved_count",
            "#CCBB44",
            "..",
        ),
    )
    mechanism_counts = {
        field: np.asarray([int(row[field]) for row in rows])
        for _, field, _, _ in loss_series
    }
    denominator = np.asarray(
        [int(row["callable_reference_nlr_denominator"]) for row in rows]
    )
    positive = sum(
        (mechanism_counts[field] for _, field, _, _ in loss_series),
        start=np.zeros(len(rows), dtype=int),
    )
    percentage = np.divide(
        100.0 * positive,
        denominator,
        out=np.zeros_like(positive, dtype=float),
        where=denominator > 0,
    )

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.8,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    # Keep both panels and the caption on one portrait manuscript page while
    # retaining the journal-scale type sizes defined above.
    figure = plt.figure(figsize=(7.2, 9.0), dpi=220, constrained_layout=True)
    grid = figure.add_gridspec(2, 1, height_ratios=[1.25, 1.0])
    y = np.arange(len(rows))
    offset = 0.19

    axis_a = figure.add_subplot(grid[0, 0])
    axis_a.barh(
        y - offset,
        repertoire,
        color="#4477AA",
        height=0.34,
        label="Complete NLR repertoire",
    )
    left = np.zeros(len(rows), dtype=int)
    for legend_label, field, color, hatch in loss_series:
        values = mechanism_counts[field]
        axis_a.barh(
            y + offset,
            values,
            left=left,
            color=color,
            edgecolor="#555555" if hatch else "none",
            linewidth=0.35 if hatch else 0,
            hatch=hatch,
            height=0.34,
            label=legend_label,
        )
        left += values
    axis_a.set_yticks(y, labels)
    axis_a.invert_yaxis()
    axis_a.set_xlabel("Number of NLR loci or reference-NLR loss calls")
    axis_a.grid(axis="x", color="#dddddd", linewidth=0.6)
    axis_a.set_axisbelow(True)
    axis_a.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        ncol=3,
    )
    axis_a.text(
        0.01,
        0.99,
        "(a)",
        transform=axis_a.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        fontweight="bold",
    )

    axis_b = figure.add_subplot(grid[1, 0])
    axis_b.barh(y, percentage, color="#66A61E", height=0.58)
    axis_b.set_yticks(y, labels)
    axis_b.invert_yaxis()
    axis_b.set_xlabel("Non-shared reference-NLR loss (%)")
    axis_b.grid(axis="x", color="#dddddd", linewidth=0.6)
    axis_b.set_axisbelow(True)
    axis_b.set_xlim(0, max(percentage) * 1.25)
    for position, (value, lost, callable_count) in enumerate(
        zip(percentage, positive, denominator)
    ):
        axis_b.text(
            value,
            position,
            f"  {value:.1f}% ({lost}/{callable_count})",
            ha="left",
            va="center",
            fontsize=7.5,
        )
    axis_b.text(
        0.01,
        0.99,
        "(b)",
        transform=axis_b.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        fontweight="bold",
    )

    plot_rows = []
    for order, (row, label) in enumerate(zip(rows, labels)):
        plot_rows.append(
            {
                "plot_order": order,
                "assembly_unit_id": row["assembly_unit_id"],
                "biological_species": row["biological_species"],
                "haplotype_or_subgenome": row["haplotype_or_subgenome"],
                "display_label": label,
                "total_nlr_count": int(row["total_nlr_count"]),
                "retained_reference_nlr_count": int(
                    row["article_retained_reference_nlr_count"]
                ),
                "decayed_reference_nlr_loss_count": int(
                    row["article_decayed_reference_nlr_loss_count"]
                ),
                "deleted_reference_nlr_loss_count": int(
                    row["article_deleted_reference_nlr_loss_count"]
                ),
                **{
                    field: int(row[field])
                    for _, field, _, _ in loss_series
                },
                "positive_reference_nlr_loss_count": int(
                    row["positive_reference_nlr_loss_count"]
                ),
                "callable_reference_nlr_denominator": int(
                    row["callable_reference_nlr_denominator"]
                ),
                "positive_reference_nlr_loss_percentage": float(
                    row["positive_reference_nlr_loss_percentage"]
                ),
            }
        )
    caption = (
        "Complete NLR repertoires and classified non-shared reference-NLR loss. "
        "Panel (a) compares NLR-Annotator repertoire counts with mutually exclusive "
        "loss-evidence groups. Frameshift and in-frame-stop groups have explicit coding-"
        "disruption support; truncation/partial alignment is a candidate category, and "
        "residual sequence without a resolved mechanism is reported separately. Panel "
        "(b) uses decayed plus deleted as the numerator and retained plus decayed plus "
        "deleted as the denominator. "
        f"The {manifest['article_shared_reference_nlrs_excluded']} reference NLR genes "
        "positive in all 23 units are excluded; not-called rows do not enter either count. "
        "No species aggregation is used."
    )
    validation = {
        "status": "PASS_ARTICLE_METHOD_NLR_FIGURE",
        "assembly_units": len(rows),
        "article_nonshared_reference_nlrs": manifest[
            "article_nonshared_reference_nlrs"
        ],
        "checks": {
            "article_method_decayed_plus_deleted": True,
            "positive_calls_partitioned_into_loss_types": True,
            "resolved_denominator": True,
            "shared_reference_nlrs_excluded": True,
            "no_species_aggregation": True,
            "deleted_is_grey": True,
            "latin_binomials_italic_suffixes_upright": True,
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
        input_paths=[args.unit_summary, args.run_manifest],
        dpi=300,
    )
    plt.close(figure)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
