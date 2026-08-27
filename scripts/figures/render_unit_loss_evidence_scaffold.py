#!/usr/bin/env python3
"""Render unit losses, decayed mechanisms, sharing, and scaffold events."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.figure_bundle import _register_arial_fonts, write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label


class FigureError(ValueError):
    """Raised when an input cannot support the publication figure."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-summary", required=True, type=Path)
    parser.add_argument("--mechanism-summary", required=True, type=Path)
    parser.add_argument("--shared-summary", required=True, type=Path)
    parser.add_argument("--scaffold-nodes", required=True, type=Path)
    parser.add_argument("--branch-summary", required=True, type=Path)
    parser.add_argument("--pattern-summary", required=True, type=Path)
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


def label_for(row: dict[str, str]) -> str:
    return format_downstream_taxon_label(
        row["biological_species"],
        (row["haplotype_or_subgenome"],),
        abbreviate_genus=True,
        separator=" ",
    )


def panel_letter(axis, letter: str) -> None:
    axis.text(
        -0.025,
        1.015,
        f"({letter})",
        transform=axis.transAxes,
        fontsize=9.5,
        fontweight="bold",
        ha="right",
        va="bottom",
    )


def main() -> int:
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    _register_arial_fonts()
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D

    manifest = json.loads(args.run_manifest.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_UNIT_RESOLVED_ARTICLE_LOSS_SCAFFOLD":
        raise SystemExit("ERROR: unit-resolved scaffold manifest is not PASS")
    units = read_tsv(args.unit_summary)
    mechanisms = read_tsv(args.mechanism_summary)
    sharing = read_tsv(args.shared_summary)
    nodes = read_tsv(args.scaffold_nodes)
    branches = read_tsv(args.branch_summary)
    patterns = read_tsv(args.pattern_summary)
    if len(units) != 23 or len(mechanisms) != 23 or len(sharing) != 23:
        raise SystemExit("ERROR: unit tables must each contain exactly 23 rows")
    units_by_id = {row["assembly_unit_id"]: row for row in units}
    mechanisms_by_id = {row["assembly_unit_id"]: row for row in mechanisms}
    sharing_by_id = {row["assembly_unit_id"]: row for row in sharing}
    if (
        len(units_by_id) != 23
        or set(units_by_id) != set(mechanisms_by_id)
        or set(units_by_id) != set(sharing_by_id)
    ):
        raise SystemExit("ERROR: unit universes differ between figure inputs")
    units.sort(
        key=lambda row: (
            row["biological_species"],
            row["haplotype_or_subgenome"],
            row["assembly_unit_id"],
        )
    )
    unit_ids = [row["assembly_unit_id"] for row in units]
    labels = [label_for(row) for row in units]
    y = np.arange(len(units))

    decayed = np.asarray([int(row["decayed"]) for row in units])
    deleted = np.asarray([int(row["deleted"]) for row in units])
    positive = decayed + deleted
    if positive.sum() != 179827:
        raise SystemExit("ERROR: historical positive counts do not close to 179,827")

    frameshift = np.asarray(
        [int(mechanisms_by_id[unit]["frameshift_only"]) for unit in unit_ids]
    )
    stop = np.asarray(
        [int(mechanisms_by_id[unit]["inframe_stop_only"]) for unit in unit_ids]
    )
    both = np.asarray(
        [
            int(mechanisms_by_id[unit]["frameshift_and_inframe_stop"])
            for unit in unit_ids
        ]
    )
    confirmed = frameshift + stop + both
    if (
        int(frameshift.sum()),
        int(stop.sum()),
        int(both.sum()),
    ) != (11559, 3258, 5071):
        raise SystemExit("ERROR: confirmed decayed mechanism totals changed")
    if np.any(confirmed > decayed):
        raise SystemExit("ERROR: confirmed decayed mechanisms exceed decayed")

    shared_decayed = np.asarray(
        [int(sharing_by_id[unit]["shared_decayed"]) for unit in unit_ids]
    )
    shared_deleted = np.asarray(
        [int(sharing_by_id[unit]["shared_deleted"]) for unit in unit_ids]
    )
    nonshared_decayed = np.asarray(
        [int(sharing_by_id[unit]["nonshared_decayed"]) for unit in unit_ids]
    )
    nonshared_deleted = np.asarray(
        [int(sharing_by_id[unit]["nonshared_deleted"]) for unit in unit_ids]
    )
    if np.any(
        shared_decayed
        + shared_deleted
        + nonshared_decayed
        + nonshared_deleted
        != positive
    ):
        raise SystemExit("ERROR: shared/non-shared components do not close")

    node_by_id = {row["node_id"]: row for row in nodes}
    branch_by_id = {row["node_id"]: row for row in branches}
    if len(node_by_id) != len(nodes) or set(node_by_id) != set(branch_by_id):
        raise SystemExit("ERROR: scaffold node and branch universes differ")
    roots = [row for row in nodes if not row["parent_node_id"]]
    terminal_nodes = [
        row for row in nodes if row["node_type"] == "unit_terminal"
    ]
    if len(roots) != 1 or len(terminal_nodes) != 23:
        raise SystemExit("ERROR: scaffold must have one root and 23 terminals")
    children: dict[str, list[str]] = defaultdict(list)
    for row in nodes:
        if row["parent_node_id"]:
            children[row["parent_node_id"]].append(row["node_id"])
    for parent in children:
        children[parent].sort(
            key=lambda node_id: int(node_by_id[node_id]["minimum_leaf_plot_order"])
        )

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.2,
            "axes.labelsize": 8.0,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "legend.fontsize": 6.5,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure = plt.figure(figsize=(7.25, 9.55), dpi=300, constrained_layout=True)
    grid = figure.add_gridspec(
        2,
        2,
        width_ratios=(1.03, 0.97),
        height_ratios=(1.0, 1.08),
    )
    axis_a = figure.add_subplot(grid[0, 0])
    axis_b = figure.add_subplot(grid[0, 1])
    axis_c = figure.add_subplot(grid[1, 0])
    axis_d = figure.add_subplot(grid[1, 1])

    axis_a.barh(y, decayed, color="#E69F00", height=0.70, label="Decayed")
    axis_a.barh(
        y,
        deleted,
        left=decayed,
        color="#8A8A8A",
        height=0.70,
        label="Deleted",
    )
    axis_a.set_yticks(y, labels)
    axis_a.invert_yaxis()
    axis_a.set_xlabel("Number of genes")
    axis_a.grid(axis="x", color="#dddddd", linewidth=0.6)
    axis_a.set_axisbelow(True)
    axis_a.set_xlim(0, positive.max() * 1.19)
    for position, value in enumerate(positive):
        axis_a.text(value, position, f"  {value:,}", va="center", fontsize=6.2)
    axis_a.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        ncol=2,
    )
    panel_letter(axis_a, "a")

    mechanism_left = np.zeros(len(units), dtype=int)
    for values, color, label in (
        (frameshift, "#4477AA", "Frameshift only"),
        (stop, "#66A61E", "Premature stop only"),
        (both, "#AA3377", "Frameshift + premature stop"),
    ):
        axis_b.barh(
            y,
            values,
            left=mechanism_left,
            color=color,
            height=0.70,
            label=label,
        )
        mechanism_left += values
    axis_b.set_yticks(y, labels)
    axis_b.set_yticklabels([])
    axis_b.tick_params(axis="y", length=0)
    axis_b.invert_yaxis()
    axis_b.set_xlabel("Decayed genes with confirmed coding disruption")
    axis_b.grid(axis="x", color="#dddddd", linewidth=0.6)
    axis_b.set_axisbelow(True)
    axis_b.set_xlim(0, max(confirmed) * 1.54)
    for position, (typed, total) in enumerate(zip(confirmed, decayed)):
        axis_b.text(
            typed,
            position,
            f"  {typed:,} ({typed / total:.1%})",
            va="center",
            fontsize=5.9,
        )
    axis_b.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        ncol=1,
    )
    panel_letter(axis_b, "b")

    shared_left = np.zeros(len(units), dtype=int)
    for values, color, label in (
        (shared_decayed, "#F3C567", "Shared decayed"),
        (shared_deleted, "#BDBDBD", "Shared deleted"),
        (nonshared_decayed, "#D97A00", "Non-shared decayed"),
        (nonshared_deleted, "#565656", "Non-shared deleted"),
    ):
        axis_c.barh(
            y,
            values,
            left=shared_left,
            color=color,
            height=0.70,
            label=label,
        )
        shared_left += values
    axis_c.set_yticks(y, labels)
    axis_c.invert_yaxis()
    axis_c.set_xlabel("Number of genes")
    axis_c.grid(axis="x", color="#dddddd", linewidth=0.6)
    axis_c.set_axisbelow(True)
    axis_c.set_xlim(0, positive.max() * 1.19)
    for position, value in enumerate(positive):
        axis_c.text(value, position, f"  {value:,}", va="center", fontsize=6.2)
    axis_c.legend(
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        ncol=2,
    )
    panel_letter(axis_c, "c")

    scaffold_unit_labels = {
        f"unit__{row['assembly_unit_id']}": label_for(row) for row in units
    }
    terminal_y = {
        row["node_id"]: float(row["minimum_leaf_plot_order"])
        for row in terminal_nodes
    }
    node_y: dict[str, float] = {}

    def calculate_y(node_id: str) -> float:
        if node_id in node_y:
            return node_y[node_id]
        if node_id in terminal_y:
            node_y[node_id] = terminal_y[node_id]
        else:
            child_positions = [calculate_y(child) for child in children[node_id]]
            node_y[node_id] = sum(child_positions) / len(child_positions)
        return node_y[node_id]

    root_id = roots[0]["node_id"]
    calculate_y(root_id)
    node_x = {row["node_id"]: int(row["depth"]) for row in nodes}
    for parent, child_ids in children.items():
        positions = [node_y[child] for child in child_ids]
        axis_d.plot(
            [node_x[parent], node_x[parent]],
            [min(positions), max(positions)],
            color="#303030",
            linewidth=1.1,
        )
        for child in child_ids:
            axis_d.plot(
                [node_x[parent], node_x[child]],
                [node_y[child], node_y[child]],
                color="#303030",
                linewidth=1.1,
            )
            event_count = int(branch_by_id[child]["loss_event_gene_count"])
            if event_count:
                color = (
                    "#1F5A7A"
                    if node_by_id[child]["node_type"] == "unit_terminal"
                    else "#8B2C4A"
                )
                axis_d.text(
                    (node_x[parent] + node_x[child]) / 2,
                    node_y[child] - 0.16,
                    f"{event_count:,}",
                    ha="center",
                    va="bottom",
                    fontsize=5.8,
                    color=color,
                )
    max_depth = max(node_x.values())
    for node_id, position in terminal_y.items():
        axis_d.text(
            max_depth + 0.18,
            position,
            scaffold_unit_labels[node_id],
            ha="left",
            va="center",
            fontsize=6.3,
        )
    root_events = int(branch_by_id[root_id]["loss_event_gene_count"])
    axis_d.text(
        node_x[root_id] - 0.05,
        node_y[root_id] - 0.35,
        f"{root_events:,}",
        ha="right",
        va="bottom",
        fontsize=6.4,
        color="#8B2C4A",
        fontweight="bold",
    )
    axis_d.legend(
        handles=[
            Line2D([0], [0], color="#1F5A7A", marker="", label="Species-specific"),
            Line2D([0], [0], color="#8B2C4A", marker="", label="Tree-node"),
        ],
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        ncol=2,
    )
    axis_d.set_xlim(-0.6, max_depth + 3.2)
    axis_d.set_ylim(len(terminal_nodes) - 0.25, -0.75)
    axis_d.axis("off")
    panel_letter(axis_d, "d")

    plot_rows: list[dict[str, object]] = []
    for order, (row, label) in enumerate(zip(units, labels)):
        unit = row["assembly_unit_id"]
        mechanism = mechanisms_by_id[unit]
        shared = sharing_by_id[unit]
        plot_rows.append(
            {
                "panel": "unit",
                "plot_order": order,
                "record_id": unit,
                "display_label": label,
                "decayed": int(row["decayed"]),
                "deleted": int(row["deleted"]),
                "frameshift_only": int(mechanism["frameshift_only"]),
                "inframe_stop_only": int(mechanism["inframe_stop_only"]),
                "frameshift_and_inframe_stop": int(
                    mechanism["frameshift_and_inframe_stop"]
                ),
                "shared_decayed": int(shared["shared_decayed"]),
                "shared_deleted": int(shared["shared_deleted"]),
                "nonshared_decayed": int(shared["nonshared_decayed"]),
                "nonshared_deleted": int(shared["nonshared_deleted"]),
                "node_type": "",
                "parent_node_id": "",
                "loss_event_gene_count": "",
            }
        )
    for order, row in enumerate(
        sorted(nodes, key=lambda item: (int(item["depth"]), item["node_id"]))
    ):
        plot_rows.append(
            {
                "panel": "scaffold_node",
                "plot_order": order,
                "record_id": row["node_id"],
                "display_label": scaffold_unit_labels.get(row["node_id"], ""),
                "decayed": "",
                "deleted": "",
                "frameshift_only": "",
                "inframe_stop_only": "",
                "frameshift_and_inframe_stop": "",
                "shared_decayed": "",
                "shared_deleted": "",
                "nonshared_decayed": "",
                "nonshared_deleted": "",
                "node_type": row["node_type"],
                "parent_node_id": row["parent_node_id"],
                "loss_event_gene_count": int(
                    branch_by_id[row["node_id"]]["loss_event_gene_count"]
                ),
            }
        )

    caption = (
        "Unit-resolved gene-loss evidence under the historical threshold rule. "
        "Panel (a) reports decayed and deleted genes for each assembly unit. "
        "Panel (b) shows only the decayed subset with confirmed frameshift, premature-"
        "stop, or combined coding disruption; labels give confirmed/total decayed. "
        "Panel (c) separates genes positive in all 23 units from other positive calls. "
        "Panel (d) places exact maximal all-positive events on a topology-only scaffold "
        "that expands multi-unit species into parallel terminals. Branch labels are "
        "event-gene counts. Genes containing any not-called unit are not forced onto a "
        "branch. The scaffold is not a newly inferred or dated 23-species phylogeny."
    )
    validation = {
        "status": "PASS_UNIT_LOSS_EVIDENCE_SCAFFOLD_FIGURE",
        "assembly_units": len(units),
        "scaffold_nodes": len(nodes),
        "article_positive_unit_gene_rows": int(positive.sum()),
        "confirmed_decayed_unit_gene_rows": int(confirmed.sum()),
        "checks": {
            "actual_gene_counts": True,
            "deleted_is_grey": True,
            "aggregate_mechanism_panel_removed": True,
            "mechanisms_intersect_decayed_only": True,
            "shared_nonshared_closure": True,
            "parallel_unit_terminals": True,
            "not_called_not_forced_to_branch": True,
            "latin_binomials_italic_suffixes_upright": True,
            "no_internal_stage_language_on_figure": True,
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
            args.mechanism_summary,
            args.shared_summary,
            args.scaffold_nodes,
            args.branch_summary,
            args.pattern_summary,
            args.run_manifest,
        ],
        dpi=300,
    )
    plt.close(figure)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
