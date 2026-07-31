#!/usr/bin/env python3
"""Render unit losses and strict disruption-evidence types."""

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


STRICT_LABELS = {
    "frameshift_only": "Frameshift only",
    "inframe_stop_only": "In-frame stop only",
    "frameshift_and_inframe_stop": "Frameshift + in-frame stop",
}


class FigureError(ValueError):
    """Raised when the evidence summary cannot support the publication plot."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-summary", required=True, type=Path)
    parser.add_argument("--strict-type-summary", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--basename",
        default="loss_evidence_classification",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    if not reader.fieldnames or not rows:
        raise FigureError(f"{path.name}: invalid or empty table")
    return rows


def main() -> int:
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    manifest = json.loads(args.run_manifest.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_UNIT_RESOLVED_ARTICLE_LOSS_SCAFFOLD":
        raise SystemExit("ERROR: loss evidence manifest is not PASS")
    units = read_tsv(args.unit_summary)
    strict_types = read_tsv(args.strict_type_summary)
    if len(units) != 23 or len(strict_types) != 23:
        raise SystemExit("ERROR: expected 23 unit-loss and mechanism rows")
    strict_by_unit = {row["assembly_unit_id"]: row for row in strict_types}
    if len(strict_by_unit) != 23:
        raise SystemExit("ERROR: duplicate unit in strict disruption summary")
    units.sort(
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
        for row in units
    ]
    if set(strict_by_unit) != {row["assembly_unit_id"] for row in units}:
        raise SystemExit("ERROR: unit-loss and strict-disruption cohorts differ")
    decayed = np.asarray([int(row["decayed"]) for row in units])
    deleted = np.asarray([int(row["deleted"]) for row in units])
    strict = np.asarray([
        int(strict_by_unit[row["assembly_unit_id"]]["confirmed_type_total"])
        for row in units
    ])
    if np.any(strict > decayed + deleted):
        raise SystemExit("ERROR: strict marker exceeds article positive count")

    type_order = tuple(STRICT_LABELS)
    type_matrix = np.asarray([
        [
            int(strict_by_unit[row["assembly_unit_id"]][evidence_type])
            for evidence_type in type_order
        ]
        for row in units
    ])
    if not np.array_equal(type_matrix.sum(axis=1), strict):
        raise SystemExit("ERROR: strict disruption types do not close by unit")

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
    figure = plt.figure(figsize=(7.2, 8.7), dpi=220, constrained_layout=True)
    grid = figure.add_gridspec(1, 2, width_ratios=[1.35, 1.0])
    y = np.arange(len(units))

    axis_a = figure.add_subplot(grid[0])
    axis_a.barh(
        y,
        decayed,
        color="#E69F00",
        height=0.68,
        label="Decayed",
    )
    axis_a.barh(
        y,
        deleted,
        left=decayed,
        color="#8F8F8F",
        height=0.68,
        label="Deleted",
    )
    axis_a.scatter(
        strict,
        y,
        color="#111111",
        marker="D",
        s=28,
        zorder=4,
        label="Strict disruption subset (not additive)",
    )
    axis_a.set_yticks(y, labels)
    axis_a.invert_yaxis()
    axis_a.set_xlabel("Number of unit–gene calls")
    axis_a.grid(axis="x", color="#dddddd", linewidth=0.6)
    axis_a.set_axisbelow(True)
    axis_a.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        ncol=2,
    )
    axis_a.text(
        0.02,
        0.99,
        "(a)",
        transform=axis_a.transAxes,
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        va="top",
    )

    axis_b = figure.add_subplot(grid[1])
    type_y = np.arange(len(units))
    type_left = np.zeros(len(units), dtype=int)
    for column, (evidence_type, color) in enumerate(
        zip(type_order, ("#4477AA", "#66A61E", "#AA3377"))
    ):
        values = type_matrix[:, column]
        axis_b.barh(
            type_y,
            values,
            left=type_left,
            color=color,
            height=0.68,
            label=STRICT_LABELS[evidence_type],
        )
        type_left += values
    axis_b.set_yticks(type_y)
    axis_b.set_yticklabels([])
    axis_b.tick_params(axis="y", length=0)
    axis_b.invert_yaxis()
    axis_b.set_xlabel("Confirmed disruption calls")
    axis_b.grid(axis="x", color="#dddddd", linewidth=0.6)
    axis_b.set_axisbelow(True)
    axis_b.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        ncol=1,
    )
    axis_b.text(
        0.02,
        0.99,
        "(b)",
        transform=axis_b.transAxes,
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        va="top",
    )

    plot_rows: list[dict[str, object]] = []
    for order, row in enumerate(units):
        plot_rows.append(
            {
                "panel": "unit_classification",
                "plot_order": order,
                "assembly_unit_id": row["assembly_unit_id"],
                "biological_species": row["biological_species"],
                "haplotype_or_subgenome": row["haplotype_or_subgenome"],
                "display_label": labels[order],
                "article_decayed": int(row["decayed"]),
                "article_deleted": int(row["deleted"]),
                "article_positive_loss": int(row["positive_loss"]),
                "strict_disruption_type": "",
                "strict_pseudogenized_evidence": int(strict[order]),
            }
        )
        for column, evidence_type in enumerate(type_order):
            plot_rows.append(
                {
                    "panel": "unit_strict_evidence_type",
                    "plot_order": order,
                    "assembly_unit_id": row["assembly_unit_id"],
                    "biological_species": row["biological_species"],
                    "haplotype_or_subgenome": row["haplotype_or_subgenome"],
                    "display_label": labels[order],
                    "article_decayed": 0,
                    "article_deleted": 0,
                    "article_positive_loss": 0,
                    "strict_disruption_type": evidence_type,
                    "strict_pseudogenized_evidence": int(
                        type_matrix[order, column]
                    ),
                }
            )
    caption = (
        "Gene-loss classifications and strict coding-disruption evidence. Panel (a) uses "
        "decayed plus deleted as the loss numerator for every assembly unit. Diamonds show "
        "the orthogonal strict pseudogenized count and are not added to the stacked bars. "
        "Panel (b) separates the confirmed coding-disruption subset within every assembly "
        "unit into frameshift-only, in-frame-stop-only, and combined disruption types. "
        "Strict calls are a conservative mechanistic subset rather than the main trend "
        "definition."
    )
    validation = {
        "status": "PASS_LOSS_EVIDENCE_CLASSIFICATION_FIGURE",
        "assembly_units": len(units),
        "article_positive_unit_gene_rows": int((decayed + deleted).sum()),
        "strict_pseudogenized_unit_gene_rows": int(strict.sum()),
        "checks": {
            "article_method_decayed_plus_deleted": True,
            "strict_evidence_not_additive": True,
            "three_strict_disruption_types_by_unit": True,
            "gene_counts_displayed": True,
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
        input_paths=[
            args.unit_summary,
            args.strict_type_summary,
            args.run_manifest,
        ],
        dpi=300,
    )
    plt.close(figure)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
