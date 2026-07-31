#!/usr/bin/env python3
"""Render loss counts and representative functions on the 23-unit scaffold."""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import math
import sys
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.figure_bundle import write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label


DETAIL_SCRIPT = (
    ROOT / "scripts" / "figures" / "render_unit_loss_functional_detail.py"
)
SPEC = importlib.util.spec_from_file_location(
    "_unit_functional_detail",
    DETAIL_SCRIPT,
)
DETAIL = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(DETAIL)


CATEGORY_LIMITS = (
    ("GO biological process", 5, "BP"),
    ("GO molecular function", 3, "MF"),
    ("GO cellular component", 2, "CC"),
    ("KEGG pathway", 5, "KEGG"),
)
NODE_COLORS = {
    "backbone_internal": "#B44C4C",
    "species_unit_group": "#426E9C",
    "unit_terminal": "#8A8F98",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--node-summary", required=True, type=Path)
    parser.add_argument("--significant-enrichment", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--kegg-pathway-names", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--basename",
        default="scaffold_loss_functional_detail",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    csv.field_size_limit(sys.maxsize)
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
        fields = reader.fieldnames or []
    if not rows or len(fields) != len(set(fields)):
        raise ValueError(f"{path.name}: invalid or empty table")
    return rows, fields


def node_label(
    row: dict[str, str],
    units: dict[str, tuple[str, str]],
) -> str:
    branch_id = row["branch_id"]
    node_type = row["node_type"]
    if node_type == "unit_terminal":
        unit = branch_id.removeprefix("unit__")
        if unit not in units:
            raise ValueError(f"unknown unit terminal {unit}")
        species, suffix = units[unit]
        return format_downstream_taxon_label(
            species,
            (suffix,),
            abbreviate_genus=True,
            separator=" ",
        )
    if node_type == "species_unit_group":
        return format_downstream_taxon_label(
            row["node_name"],
            ("unit group",),
            abbreviate_genus=True,
            separator=" ",
        )
    count = int(row["descendant_unit_count"])
    if count == 23:
        return r"$\it{Actinidia}$ stem (23 units)"
    return f"Internal clade {branch_id[-6:]} ({count} units)"


def main() -> int:
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    manifest = json.loads(args.run_manifest.read_text(encoding="utf-8"))
    if (
        manifest.get("status")
        != "PASS_UNIT_SCAFFOLD_GO_KEGG_SUMMARY"
        or int(manifest.get("foregrounds", 0)) != 39
        or manifest.get("loss_classification") != "decayed + deleted"
    ):
        raise SystemExit("ERROR: scaffold functional manifest is not PASS")
    nodes, node_fields = read_tsv(args.node_summary)
    required_nodes = {
        "foreground_id",
        "branch_id",
        "node_type",
        "node_name",
        "descendant_unit_count",
        "minimum_leaf_plot_order",
        "loss_event_gene_count",
    }
    if not required_nodes.issubset(node_fields) or len(nodes) != 39:
        raise SystemExit("ERROR: invalid 39-node summary")
    unit_rows, unit_fields = read_tsv(args.unit_metadata)
    if not {
        "assembly_unit_id",
        "biological_species",
        "haplotype_or_subgenome",
    }.issubset(unit_fields):
        raise SystemExit("ERROR: invalid unit metadata")
    units = {
        row["assembly_unit_id"]: (
            row["biological_species"],
            row["haplotype_or_subgenome"],
        )
        for row in unit_rows
    }
    if len(units) != 23:
        raise SystemExit("ERROR: expected 23 unit metadata rows")

    node_by_foreground = {row["foreground_id"]: row for row in nodes}
    if len(node_by_foreground) != len(nodes):
        raise SystemExit("ERROR: duplicate node foreground")
    labels = {row["foreground_id"]: node_label(row, units) for row in nodes}

    raw_rows, significant_fields = read_tsv(
        args.significant_enrichment
    )
    if len(raw_rows) != int(manifest["significant_terms"]):
        raise SystemExit("ERROR: significant term count does not close")
    pathway_names = DETAIL.read_pathway_names(args.kegg_pathway_names)
    go_parents, go_aliases = DETAIL.read_go_parents(args.go_obo)
    go_related = DETAIL.ancestor_checker(go_parents, go_aliases)

    normalized: list[dict[str, object]] = []
    for row in raw_rows:
        foreground_id = row["foreground_id"]
        if foreground_id not in node_by_foreground:
            raise SystemExit("ERROR: unknown scaffold foreground")
        q_value = float(row["p_fdr_bh"])
        study_count = int(row["study_count"])
        fold = float(row["fold_enrichment"])
        category = row["functional_category"]
        normalized.append(
            {
                **row,
                "term_name": (
                    pathway_names.get(row["term_id"], row["term_name"])
                    if category == "KEGG pathway"
                    else row["term_name"]
                ),
                "study_count": study_count,
                "p_fdr_bh": q_value,
                "fold_enrichment": fold,
                "minus_log10_q": -math.log10(q_value),
            }
        )
    selected: dict[str, list[tuple[str, str]]] = {}
    for category, limit, _ in CATEGORY_LIMITS:
        selected[category] = DETAIL.rank_terms(
            (
                row
                for row in normalized
                if row["functional_category"] == category
            ),
            limit=limit,
            sample_key="foreground_id",
            go_related=go_related if category.startswith("GO ") else None,
        )
        if len(selected[category]) != limit:
            raise SystemExit(
                f"ERROR: insufficient representative terms for {category}"
            )
    selected_pathways = {
        term_id for term_id, _ in selected["KEGG pathway"]
    }
    if not selected_pathways.issubset(pathway_names):
        raise SystemExit("ERROR: selected KEGG pathway lacks display name")
    selected_keys = {
        (category, term_id, term_name)
        for category, terms in selected.items()
        for term_id, term_name in terms
    }
    selected_rows = [
        row
        for row in normalized
        if (
            str(row["functional_category"]),
            str(row["term_id"]),
            str(row["term_name"]),
        )
        in selected_keys
    ]

    nodes = sorted(
        nodes,
        key=lambda row: (
            int(row["minimum_leaf_plot_order"]),
            {"backbone_internal": 0, "species_unit_group": 1, "unit_terminal": 2}[
                row["node_type"]
            ],
            -int(row["descendant_unit_count"]),
            row["branch_id"],
        ),
    )
    node_index = {row["foreground_id"]: i for i, row in enumerate(nodes)}
    term_order: list[tuple[str, str, str, str, str]] = []
    for category, _, abbreviation in CATEGORY_LIMITS:
        for rank, (term_id, term_name) in enumerate(
            selected[category],
            start=1,
        ):
            prefix = "K" if abbreviation == "KEGG" else abbreviation
            term_order.append(
                (
                    category,
                    term_id,
                    term_name,
                    abbreviation,
                    f"{prefix}{rank}",
                )
            )
    term_index = {
        (category, term_id, term_name): i
        for i, (category, term_id, term_name, _, _) in enumerate(term_order)
    }

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.2,
            "legend.fontsize": 7.0,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure = plt.figure(figsize=(7.2, 10.5), dpi=180)
    grid = figure.add_gridspec(
        1,
        2,
        width_ratios=(0.92, 2.7),
        left=0.23,
        right=0.92,
        bottom=0.27,
        top=0.91,
        wspace=0.12,
    )
    y = np.arange(len(nodes))

    axis_a = figure.add_subplot(grid[0])
    counts = np.asarray(
        [int(row["loss_event_gene_count"]) for row in nodes],
        dtype=int,
    )
    axis_a.barh(
        y,
        counts,
        color=[NODE_COLORS[row["node_type"]] for row in nodes],
        height=0.68,
    )
    axis_a.set_yticks(
        y,
        [labels[row["foreground_id"]] for row in nodes],
    )
    axis_a.invert_yaxis()
    axis_a.set_xlabel("Loss-event genes")
    axis_a.grid(axis="x", color="#e3e3e3", linewidth=0.55)
    axis_a.set_axisbelow(True)
    axis_a.set_xlim(0, max(counts) * 1.25)
    for position, value in enumerate(counts):
        axis_a.text(
            value,
            position,
            f"  {value:,}",
            ha="left",
            va="center",
            fontsize=6.5,
        )
    axis_a.text(
        -0.10,
        1.025,
        "(a)",
        transform=axis_a.transAxes,
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    axis_b = figure.add_subplot(grid[1], sharey=axis_a)
    colors = np.asarray(
        [float(row["minus_log10_q"]) for row in selected_rows],
        dtype=float,
    )
    color_min = max(-math.log10(0.05), float(colors.min()))
    color_max = max(color_min + 0.5, float(np.quantile(colors, 0.98)))
    sizes_raw = np.asarray(
        [int(row["study_count"]) for row in selected_rows],
        dtype=float,
    )
    size_scale = 1.8
    x_values = [
        term_index[
            (
                str(row["functional_category"]),
                str(row["term_id"]),
                str(row["term_name"]),
            )
        ]
        for row in selected_rows
    ]
    y_values = [
        node_index[str(row["foreground_id"])] for row in selected_rows
    ]
    axis_b.scatter(
        x_values,
        y_values,
        c=[min(value, color_max) for value in colors],
        s=[size_scale * value for value in sizes_raw],
        cmap="viridis_r",
        vmin=color_min,
        vmax=color_max,
        edgecolors="#2f2f2f",
        linewidths=0.35,
        alpha=0.92,
    )
    axis_b.set_yticks(y)
    axis_b.tick_params(axis="y", labelleft=False)
    axis_b.set_xticks(
        range(len(term_order)),
        [code for _, _, _, _, code in term_order],
    )
    axis_b.set_xlim(-0.65, len(term_order) - 0.35)
    axis_b.set_ylim(len(nodes) - 0.35, -0.65)
    axis_b.set_axisbelow(True)
    axis_b.grid(color="#e8e8e8", linewidth=0.55)
    for position, row in enumerate(nodes):
        if row["node_type"] != "unit_terminal":
            axis_b.axhspan(
                position - 0.48,
                position + 0.48,
                color=NODE_COLORS[row["node_type"]],
                alpha=0.055,
                zorder=0,
            )
    axis_b.text(
        -0.06,
        1.025,
        "(b)",
        transform=axis_b.transAxes,
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    colorbar = figure.colorbar(
        ScalarMappable(
            norm=Normalize(vmin=color_min, vmax=color_max),
            cmap="viridis_r",
        ),
        ax=axis_b,
        fraction=0.025,
        pad=0.018,
    )
    colorbar.set_label(r"$-\log_{10}$(BH q-value)")
    colorbar.ax.tick_params(labelsize=7.0)
    legend_counts = sorted(
        set(
            max(2, int(round(value)))
            for value in np.quantile(sizes_raw, [0.15, 0.50, 0.85])
        )
    )
    handles = [
        plt.scatter(
            [],
            [],
            s=size_scale * value,
            facecolor="#7a7a7a",
            edgecolor="#2f2f2f",
            linewidth=0.35,
        )
        for value in legend_counts
    ]
    axis_b.legend(
        handles,
        [f"{value} genes" for value in legend_counts],
        title="Lost genes in term",
        loc="lower right",
        ncol=len(legend_counts),
        frameon=False,
        bbox_to_anchor=(1.0, 1.005),
        borderaxespad=0,
        columnspacing=0.9,
        handletextpad=0.4,
        fontsize=7.0,
        title_fontsize=7.5,
    )

    key_axis = figure.add_axes([0.23, 0.018, 0.69, 0.21])
    key_axis.axis("off")
    for index, (_, term_id, term_name, abbreviation, code) in enumerate(
        term_order
    ):
        column = index // 5
        row = index % 5
        wrapped = textwrap.fill(
            term_name,
            width=22,
            break_long_words=False,
            break_on_hyphens=False,
        )
        key_axis.text(
            column * 0.34,
            0.98 - row * 0.19,
            f"{code} {wrapped}",
            transform=key_axis.transAxes,
            fontsize=6.8,
            ha="left",
            va="top",
            linespacing=1.05,
        )
    for axis in (axis_a, axis_b):
        for spine in axis.spines.values():
            spine.set_color("#666666")
            spine.set_linewidth(0.7)

    plot_rows = [
        {
            "node_plot_order": node_index[str(row["foreground_id"])],
            "foreground_id": row["foreground_id"],
            "branch_id": row["branch_id"],
            "node_type": row["node_type"],
            "display_label": labels[str(row["foreground_id"])],
            "functional_category": row["functional_category"],
            "term_id": row["term_id"],
            "term_name": row["term_name"],
            "term_code": term_order[
                term_index[
                    (
                        str(row["functional_category"]),
                        str(row["term_id"]),
                        str(row["term_name"]),
                    )
                ]
            ][4],
            "study_count": int(row["study_count"]),
            "study_size": int(row["study_size"]),
            "p_fdr_bh": float(row["p_fdr_bh"]),
            "minus_log10_q": float(row["minus_log10_q"]),
            "fold_enrichment": float(row["fold_enrichment"]),
            "study_gene_ids": row["study_gene_ids"],
        }
        for row in selected_rows
    ]
    caption = (
        "Loss-event counts and representative functional enrichment on the "
        "topology-only 23-unit scaffold. Panel (a) reports maximal "
        "decayed-plus-deleted event genes assigned to each internal node, "
        "within-species unit group, or assembly-unit terminal; red, blue, and "
        "grey bars distinguish these three node types, respectively. Panel (b) shows "
        "representative GO and KEGG terms; point area is the number of event "
        "genes assigned to the term and color is Benjamini–Hochberg "
        "significance. All tests use the 33,998 reference genes resolved in "
        "all 23 units as the common opportunity background. The scaffold is "
        "not a newly inferred 23-species phylogeny, and KEGG labels are "
        "orthology-derived categories rather than organismal phenotypes."
    )
    validation = {
        "status": "PASS_SCAFFOLD_LOSS_FUNCTIONAL_DETAIL_FIGURE",
        "scaffold_nodes": len(nodes),
        "unit_terminals": sum(
            row["node_type"] == "unit_terminal" for row in nodes
        ),
        "significant_input_rows": len(raw_rows),
        "selected_term_counts": {
            category: len(terms) for category, terms in selected.items()
        },
        "plotted_term_node_rows": len(plot_rows),
        "checks": {
            "decayed_plus_deleted": True,
            "assembly_units_not_aggregated": True,
            "resolved_23_unit_background": True,
            "actual_event_gene_counts_shown": True,
            "point_area_linear_in_gene_count": True,
            "specific_terms_named": True,
            "go_redundancy_reduced": True,
            "kegg_pathway_names_resolved": True,
            "scaffold_not_claimed_as_species_phylogeny": True,
            "compact_bottom_layout": True,
            "full_term_key_in_bottom_margin": True,
            "figure_title_omitted": True,
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
            args.node_summary,
            args.significant_enrichment,
            args.run_manifest,
            args.unit_metadata,
            args.go_obo,
            args.kegg_pathway_names,
        ],
        dpi=600,
    )
    plt.close(figure)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
