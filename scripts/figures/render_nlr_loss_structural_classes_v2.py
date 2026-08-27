#!/usr/bin/env python3
"""Render extant, lost, and branch-mapped NLR structural classes."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.figure_bundle import _register_arial_fonts, write_figure_bundle
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
NLR_CLASSES = tuple(item[0] for item in NLR_CLASS_SERIES)
NLR_COLORS = dict(NLR_CLASS_SERIES)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repertoire-class-summary", required=True, type=Path)
    parser.add_argument("--class-loss-rates", required=True, type=Path)
    parser.add_argument("--reference-classes", required=True, type=Path)
    parser.add_argument("--tree-loss-events", required=True, type=Path)
    parser.add_argument("--branch-metadata", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--basename",
        default="lost_nlr_structural_classes",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    if not reader.fieldnames or not rows:
        raise ValueError(f"{path.name}: invalid or empty table")
    return rows


def short_lineage_label(lineage: str) -> str:
    if lineage == "Actinidia x zhejiangensis parental lineage A":
        return (
            r"$\mathit{A.\ zhejiangensis}$ "
            r"$\mathrm{A}$"
        )
    if lineage == "Actinidia x zhejiangensis parental lineage B":
        return (
            r"$\mathit{A.\ zhejiangensis}$ "
            r"$\mathrm{B}$"
        )
    parts = lineage.split()
    if len(parts) == 2 and parts[0] == "Actinidia":
        return rf"$\mathit{{A.\ {parts[1]}}}$"
    return lineage


def unit_label(species: str, suffix: str) -> str:
    normalized = suffix.strip().lower()
    if species.startswith("Actinidia x zhejiangensis parental lineage"):
        return format_downstream_taxon_label(
            species,
            (),
            abbreviate_genus=True,
            separator=" ",
        )
    base = short_lineage_label(species)
    if normalized in {"", "unphased", "unresolved_polyploid_unit", "actinidiabase v1"}:
        return base
    return base + rf" $\mathrm{{{suffix}}}$"


def branch_label(row: dict[str, str]) -> str:
    descendants = row["descendant_lineages"].split(";")
    if row["branch_type"] == "terminal":
        return short_lineage_label(descendants[0])
    if row["branch_id"] == "internal__actinidia_all":
        return r"$\mathit{Actinidia}$ stem"
    if len(descendants) == 2:
        return " + ".join(short_lineage_label(item) for item in descendants)
    return f"{len(descendants)}-lineage clade"


def main() -> int:
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    _register_arial_fonts()
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.gridspec import GridSpecFromSubplotSpec

    manifest = json.loads(args.run_manifest.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_LOST_NLR_STRUCTURAL_CLASSIFICATION":
        raise SystemExit("ERROR: NLR structural-class manifest is not PASS")

    repertoire_rows = read_tsv(args.repertoire_class_summary)
    rate_rows = read_tsv(args.class_loss_rates)
    reference_rows = read_tsv(args.reference_classes)
    event_rows = read_tsv(args.tree_loss_events)
    branch_rows = read_tsv(args.branch_metadata)

    expected_unit_class_rows = int(manifest["assembly_units"]) * len(NLR_CLASSES)
    if (
        len(repertoire_rows) != expected_unit_class_rows
        or len(rate_rows) != expected_unit_class_rows
    ):
        raise SystemExit("ERROR: unit-by-NLR-class grid changed")

    unit_metadata: dict[str, tuple[str, str]] = {}
    repertoire: dict[tuple[str, str], int] = {}
    repertoire_totals: dict[str, int] = {}
    for row in repertoire_rows:
        unit = row["assembly_unit_id"]
        nlr_class = row["reference_nlr_class"]
        if nlr_class not in NLR_COLORS:
            raise SystemExit(f"ERROR: unexpected NLR class {nlr_class!r}")
        metadata = (row["biological_species"], row["haplotype_or_subgenome"])
        if unit in unit_metadata and unit_metadata[unit] != metadata:
            raise SystemExit(f"ERROR: inconsistent unit metadata for {unit}")
        unit_metadata[unit] = metadata
        key = (unit, nlr_class)
        if key in repertoire:
            raise SystemExit(f"ERROR: duplicate repertoire row {key}")
        repertoire[key] = int(row["nlr_gene_count"])
        total = int(row["total_nlr_count"])
        if unit in repertoire_totals and repertoire_totals[unit] != total:
            raise SystemExit(f"ERROR: inconsistent repertoire total for {unit}")
        repertoire_totals[unit] = total

    loss: dict[tuple[str, str], int] = {}
    for row in rate_rows:
        unit = row["assembly_unit_id"]
        nlr_class = row["reference_nlr_class"]
        key = (unit, nlr_class)
        if (
            unit not in unit_metadata
            or nlr_class not in NLR_COLORS
            or key in loss
        ):
            raise SystemExit("ERROR: invalid class-loss row")
        loss[key] = int(row["all_positive_loss_count"])

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
        unit_label(
            unit_metadata[unit][0],
            unit_metadata[unit][1],
        )
        for unit in units
    ]

    repertoire_matrix = np.asarray(
        [
            [repertoire[(unit, nlr_class)] for nlr_class in NLR_CLASSES]
            for unit in units
        ],
        dtype=int,
    )
    loss_matrix = np.asarray(
        [
            [loss[(unit, nlr_class)] for nlr_class in NLR_CLASSES]
            for unit in units
        ],
        dtype=int,
    )
    if (
        int(repertoire_matrix.sum()) != int(manifest["repertoire_nlr_calls"])
        or any(
            int(total) != repertoire_totals[unit]
            for unit, total in zip(units, repertoire_matrix.sum(axis=1))
        )
    ):
        raise SystemExit("ERROR: NLR repertoire totals do not close")

    relative_burden = np.divide(
        100.0 * loss_matrix,
        repertoire_matrix + loss_matrix,
        out=np.full(loss_matrix.shape, np.nan, dtype=float),
        where=(repertoire_matrix + loss_matrix) > 0,
    )

    reference_classes: dict[str, str] = {}
    for row in reference_rows:
        gene = row["reference_nlr_id"]
        nlr_class = row["reference_nlr_class"]
        if (
            not gene
            or nlr_class not in NLR_COLORS
            or gene in reference_classes
        ):
            raise SystemExit("ERROR: invalid reference-NLR class row")
        reference_classes[gene] = nlr_class
    if len(reference_classes) != int(manifest["reference_nlrs"]):
        raise SystemExit("ERROR: reference-NLR universe changed")

    branch_metadata: dict[str, dict[str, str]] = {}
    branch_order: dict[str, int] = {}
    for row in branch_rows:
        branch = row["branch_id"]
        if (
            branch in branch_metadata
            or row["branch_type"] not in {"terminal", "internal"}
        ):
            raise SystemExit("ERROR: invalid branch metadata")
        branch_metadata[branch] = row
        branch_order[branch] = int(row["plot_order"])

    branch_class_counts: dict[str, Counter[str]] = {
        branch: Counter() for branch in branch_metadata
    }
    seen_events: set[tuple[str, str]] = set()
    for row in event_rows:
        gene = row["reference_gene_id"]
        if gene not in reference_classes:
            continue
        branch = row["branch_id"]
        key = (gene, branch)
        if branch not in branch_metadata or key in seen_events:
            raise SystemExit("ERROR: invalid or duplicate NLR tree event")
        seen_events.add(key)
        branch_class_counts[branch][reference_classes[gene]] += 1

    terminal_branches = sorted(
        (
            branch
            for branch, row in branch_metadata.items()
            if row["branch_type"] == "terminal"
        ),
        key=lambda branch: branch_order[branch],
    )
    internal_branches = sorted(
        (
            branch
            for branch, row in branch_metadata.items()
            if row["branch_type"] == "internal"
            and sum(branch_class_counts[branch].values()) > 0
        ),
        key=lambda branch: branch_order[branch],
    )
    if len(terminal_branches) != 13:
        raise SystemExit("ERROR: expected 13 biological-lineage terminals")

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.2,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    figure = plt.figure(figsize=(8.6, 10.8), dpi=240)
    outer = figure.add_gridspec(
        2,
        1,
        height_ratios=(2.05, 1.0),
        hspace=0.27,
        left=0.075,
        right=0.975,
        bottom=0.105,
        top=0.91,
    )
    top = GridSpecFromSubplotSpec(
        1,
        4,
        subplot_spec=outer[0],
        width_ratios=(2.25, 1.65, 1.40, 2.85),
        wspace=0.045,
    )
    repertoire_axis = figure.add_subplot(top[0])
    label_axis = figure.add_subplot(top[1])
    loss_axis = figure.add_subplot(top[2])
    burden_axis = figure.add_subplot(top[3])
    y = np.arange(len(units))

    def clean_bar_axis(axis) -> None:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_visible(False)
        axis.grid(axis="x", color="#dddddd", linewidth=0.55)
        axis.set_axisbelow(True)
        axis.tick_params(axis="y", left=False, labelleft=False)
        axis.set_ylim(len(units) - 0.35, -0.65)

    left = np.zeros(len(units), dtype=int)
    legend_handles = []
    for class_index, (nlr_class, color) in enumerate(NLR_CLASS_SERIES):
        values = repertoire_matrix[:, class_index]
        bars = repertoire_axis.barh(
            y,
            values,
            left=left,
            height=0.72,
            color=color,
            edgecolor="white",
            linewidth=0.25,
            label=nlr_class,
        )
        legend_handles.append(bars[0])
        left += values
    clean_bar_axis(repertoire_axis)
    repertoire_axis.invert_xaxis()
    repertoire_axis.set_xlabel("Annotated NLR genes")
    repertoire_axis.text(
        -0.02,
        1.025,
        "(a)",
        transform=repertoire_axis.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
    )
    for position, total in enumerate(repertoire_matrix.sum(axis=1)):
        repertoire_axis.text(
            total + 4,
            position,
            str(int(total)),
            ha="right",
            va="center",
            fontsize=6.8,
        )

    label_axis.set_xlim(0, 1)
    label_axis.set_ylim(len(units) - 0.35, -0.65)
    label_axis.axis("off")
    for position, label in enumerate(labels):
        label_axis.text(
            0.5,
            position,
            label,
            ha="center",
            va="center",
            fontsize=7.6,
        )

    left = np.zeros(len(units), dtype=int)
    for class_index, (nlr_class, color) in enumerate(NLR_CLASS_SERIES):
        values = loss_matrix[:, class_index]
        loss_axis.barh(
            y,
            values,
            left=left,
            height=0.72,
            color=color,
            edgecolor="white",
            linewidth=0.25,
        )
        left += values
    clean_bar_axis(loss_axis)
    loss_axis.set_xlim(0, max(loss_matrix.sum(axis=1)) * 1.13)
    loss_axis.set_xlabel("Inferred NLR losses")
    for position, total in enumerate(loss_matrix.sum(axis=1)):
        loss_axis.text(
            total + 1.5,
            position,
            str(int(total)),
            ha="left",
            va="center",
            fontsize=6.8,
        )

    burden_cmap = plt.get_cmap("YlOrRd").copy()
    burden_cmap.set_bad("#eeeeee")
    burden_image = burden_axis.imshow(
        np.ma.masked_invalid(relative_burden),
        aspect="auto",
        interpolation="nearest",
        cmap=burden_cmap,
        vmin=0,
        vmax=100,
    )
    burden_axis.set_ylim(len(units) - 0.5, -0.5)
    burden_axis.set_yticks([])
    burden_axis.set_xticks(
        np.arange(len(NLR_CLASSES)),
        NLR_CLASSES,
        rotation=54,
        ha="left",
        rotation_mode="anchor",
    )
    burden_axis.tick_params(
        axis="x",
        top=True,
        labeltop=True,
        bottom=False,
        labelbottom=False,
        pad=1.5,
    )
    burden_axis.set_xticks(np.arange(-0.5, len(NLR_CLASSES), 1), minor=True)
    burden_axis.set_yticks(np.arange(-0.5, len(units), 1), minor=True)
    burden_axis.grid(which="minor", color="white", linewidth=0.55)
    burden_axis.tick_params(which="minor", bottom=False, left=False)
    burden_axis.text(
        -0.10,
        1.025,
        "(b)",
        transform=burden_axis.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
    )
    for row_index in range(len(units)):
        for column_index in range(len(NLR_CLASSES)):
            value = relative_burden[row_index, column_index]
            if np.isnan(value):
                label = "–"
                color = "#555555"
            else:
                label = f"{value:.0f}"
                color = "white" if value >= 58 else "black"
            burden_axis.text(
                column_index,
                row_index,
                label,
                ha="center",
                va="center",
                fontsize=6.1,
                color=color,
            )
    colorbar_axis = burden_axis.inset_axes([0.24, -0.105, 0.64, 0.026])
    colorbar = figure.colorbar(
        burden_image,
        cax=colorbar_axis,
        orientation="horizontal",
    )
    colorbar.set_label(
        "Relative loss burden: loss / (annotated + loss), %",
        fontsize=7.5,
        labelpad=2,
    )
    colorbar.ax.tick_params(labelsize=6.5, pad=1)

    figure.legend(
        handles=legend_handles,
        labels=NLR_CLASSES,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.50, 0.985),
        ncol=5,
        columnspacing=1.15,
        handlelength=1.4,
    )

    bottom = GridSpecFromSubplotSpec(
        1,
        2,
        subplot_spec=outer[1],
        width_ratios=(2.05, 1.0),
        wspace=0.24,
    )
    terminal_axis = figure.add_subplot(bottom[0])
    internal_axis = figure.add_subplot(bottom[1])

    def draw_branch_composition(axis, branches: list[str]) -> None:
        totals = np.asarray(
            [sum(branch_class_counts[branch].values()) for branch in branches],
            dtype=int,
        )
        x = np.arange(len(branches))
        base = np.zeros(len(branches), dtype=float)
        for nlr_class, color in NLR_CLASS_SERIES:
            raw = np.asarray(
                [branch_class_counts[branch][nlr_class] for branch in branches],
                dtype=float,
            )
            values = np.divide(
                100.0 * raw,
                totals,
                out=np.zeros_like(raw),
                where=totals > 0,
            )
            axis.bar(
                x,
                values,
                bottom=base,
                width=0.75,
                color=color,
                edgecolor="white",
                linewidth=0.25,
            )
            base += values
        for position, total in enumerate(totals):
            axis.text(
                position,
                102.0 if total else 2.0,
                str(int(total)),
                ha="center",
                va="bottom",
                fontsize=6.8,
            )
        axis.set_ylim(0, 111)
        axis.set_xlim(-0.65, len(branches) - 0.35)
        axis.grid(axis="y", color="#dddddd", linewidth=0.55)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.set_xticks(
            x,
            [branch_label(branch_metadata[branch]) for branch in branches],
            rotation=50,
            ha="right",
            rotation_mode="anchor",
        )

    draw_branch_composition(terminal_axis, terminal_branches)
    terminal_axis.set_ylabel("NLR loss-event composition (%)")
    terminal_axis.set_xlabel("Species-specific events")
    terminal_axis.text(
        -0.055,
        1.035,
        "(c)",
        transform=terminal_axis.transAxes,
        fontsize=10.5,
        fontweight="bold",
        va="bottom",
    )

    draw_branch_composition(internal_axis, internal_branches)
    internal_axis.set_xlabel("Tree-node events with NLR losses")
    internal_axis.set_yticklabels([])
    internal_axis.tick_params(axis="y", length=0)

    plot_rows: list[dict[str, object]] = []
    for unit, display_label in zip(units, labels):
        for class_index, nlr_class in enumerate(NLR_CLASSES):
            extant = int(repertoire_matrix[units.index(unit), class_index])
            lost = int(loss_matrix[units.index(unit), class_index])
            burden = relative_burden[units.index(unit), class_index]
            plot_rows.append(
                {
                    "panel": "A_B_unit_extant_loss",
                    "assembly_unit_id": unit,
                    "biological_species": unit_metadata[unit][0],
                    "haplotype_or_subgenome": unit_metadata[unit][1],
                    "display_label": display_label,
                    "branch_id": "",
                    "branch_type": "",
                    "descendant_lineages": "",
                    "reference_nlr_class": nlr_class,
                    "annotated_nlr_count": extant,
                    "inferred_loss_count": lost,
                    "relative_loss_burden_percentage": (
                        "" if np.isnan(burden) else f"{burden:.6f}"
                    ),
                    "branch_loss_event_count": "",
                    "branch_class_composition_percentage": "",
                }
            )
    for branch in terminal_branches + internal_branches:
        total = sum(branch_class_counts[branch].values())
        for nlr_class in NLR_CLASSES:
            count = branch_class_counts[branch][nlr_class]
            percentage = 100.0 * count / total if total else None
            plot_rows.append(
                {
                    "panel": "C_tree_loss_events",
                    "assembly_unit_id": "",
                    "biological_species": "",
                    "haplotype_or_subgenome": "",
                    "display_label": branch_label(branch_metadata[branch]),
                    "branch_id": branch,
                    "branch_type": branch_metadata[branch]["branch_type"],
                    "descendant_lineages": branch_metadata[branch][
                        "descendant_lineages"
                    ],
                    "reference_nlr_class": nlr_class,
                    "annotated_nlr_count": "",
                    "inferred_loss_count": "",
                    "relative_loss_burden_percentage": "",
                    "branch_loss_event_count": count,
                    "branch_class_composition_percentage": (
                        "" if percentage is None else f"{percentage:.6f}"
                    ),
                }
            )

    caption = (
        "(a) Structural composition of complete NLR-Annotator repertoires "
        "(left) and inferred reference-NLR losses (right) for each of the 23 "
        "assembly units. Inferred loss combines shared and non-shared "
        "decayed-plus-deleted calls; totals are printed outside the bars. "
        "(b) Descriptive class-specific loss burden, calculated as inferred "
        "loss divided by inferred loss plus annotated target-assembly NLRs "
        "of the same structural class. This composition metric compares two "
        "NLR gene universes and is not an orthology-resolved evolutionary "
        "loss rate. (c) Structural composition of exact tree-mapped NLR loss "
        "events. For multi-assembly species, a species-specific event requires all "
        "assigned haplotypes or subgenomes to be decayed or deleted. Bar "
        "heights are normalized to 100%, with event totals printed above; "
        "tree-node branches without an NLR event are omitted. Structural "
        "classes are inherited from the matching Clematoclethra scandens "
        "reference NLR. Not-called comparisons are excluded."
    )
    validation = {
        "schema_version": "2.0",
        "status": "PASS_NLR_LOSS_STRUCTURAL_CLASS_FIGURE_V2",
        "assembly_units": len(units),
        "reference_nlrs": len(reference_classes),
        "annotated_repertoire_nlrs": int(repertoire_matrix.sum()),
        "unit_level_total_loss_calls": int(loss_matrix.sum()),
        "tree_nlr_loss_events": len(seen_events),
        "terminal_tree_nlr_loss_events": sum(
            sum(branch_class_counts[branch].values())
            for branch in terminal_branches
        ),
        "internal_tree_nlr_loss_events": sum(
            sum(branch_class_counts[branch].values())
            for branch in internal_branches
        ),
        "checks": {
            "shared_and_nonshared_loss_combined_per_unit": True,
            "complete_target_repertoires_shown": True,
            "relative_burden_uses_loss_plus_extant_denominator": True,
            "tree_events_use_complete_biological_lineage_loss": True,
            "multi_assembly_species_merged_for_tree_events": True,
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
            args.repertoire_class_summary,
            args.class_loss_rates,
            args.reference_classes,
            args.tree_loss_events,
            args.branch_metadata,
            args.run_manifest,
        ],
        dpi=300,
    )
    plt.close(figure)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
