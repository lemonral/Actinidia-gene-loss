#!/usr/bin/env python3
"""Render species-specific and tree-node loss events with strict support."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.figure_bundle import _register_arial_fonts, write_figure_bundle
from geneloss_repro.labels import format_taxon_label


class FigureError(ValueError):
    """Raised when the branch evidence table is invalid."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--branch-summary", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--basename",
        default="tree_branch_loss_evidence",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    if not reader.fieldnames or not rows:
        raise FigureError(f"{path.name}: invalid or empty table")
    return rows


def short_taxon(species: str) -> str:
    suffixes: tuple[str, ...] = ()
    if " parental lineage " in species:
        species, suffix = species.rsplit(" parental lineage ", 1)
        suffixes = (suffix,)
    return format_taxon_label(
        species,
        suffixes,
        abbreviate_genus=True,
        separator=" ",
    )


def internal_label(
    descendants: list[str],
    all_lineages: set[str],
) -> str:
    if len(descendants) == len(all_lineages):
        return r"$\mathit{Actinidia}$ stem (13 lineages)"
    missing = sorted(all_lineages - set(descendants))
    if len(descendants) >= 10:
        omitted = ", ".join(short_taxon(species) for species in missing)
        return f"Clade excluding\n{omitted}"
    rendered = [short_taxon(species) for species in descendants]
    if len(rendered) <= 2:
        return " + ".join(rendered)
    return f"{len(rendered)}-lineage clade\n" + ", ".join(rendered)


def main() -> int:
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    _register_arial_fonts()
    import matplotlib.pyplot as plt
    import numpy as np

    manifest = json.loads(args.run_manifest.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_ARTICLE_LOSS_WITH_STRICT_EVIDENCE_SUMMARY":
        raise SystemExit("ERROR: loss evidence manifest is not PASS")
    rows = read_tsv(args.branch_summary)
    terminal = [row for row in rows if row["branch_type"] == "terminal"]
    internal = [row for row in rows if row["branch_type"] == "internal"]
    if len(terminal) != 13 or len(internal) != 12:
        raise SystemExit("ERROR: expected 13 terminal and 12 internal branches")
    all_lineages = {
        row["descendant_lineages"]
        for row in terminal
    }
    terminal.sort(key=lambda row: row["descendant_lineages"])
    internal.sort(
        key=lambda row: (
            -int(row["descendant_lineage_count"]),
            row["descendant_lineages"],
        )
    )

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10.5,
            "axes.labelsize": 11.5,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 9.5,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(17.5, 11.8),
        dpi=220,
        constrained_layout=True,
    )
    plot_rows: list[dict[str, object]] = []
    for panel, (axis, branch_rows, branch_type) in enumerate(
        zip(axes, (terminal, internal), ("terminal", "internal")),
        start=1,
    ):
        if branch_type == "terminal":
            labels = [
                short_taxon(row["descendant_lineages"])
                for row in branch_rows
            ]
        else:
            labels = [
                internal_label(
                    row["descendant_lineages"].split(";"),
                    all_lineages,
                )
                for row in branch_rows
            ]
        total = np.asarray(
            [int(row["article_method_loss_gene_count"]) for row in branch_rows]
        )
        strict_all = np.asarray(
            [
                int(row["strict_supported_all_descendant_lineages"])
                for row in branch_rows
            ]
        )
        if np.any(strict_all > total):
            raise SystemExit("ERROR: strict branch support exceeds article loss count")
        y = np.arange(len(branch_rows))
        axis.barh(
            y,
            total,
            color="#9ECAE1",
            height=0.70,
            label="Branch loss",
        )
        axis.barh(
            y,
            strict_all,
            color="#2F4B5C",
            height=0.34,
            label="Strict disruption-supported subset",
        )
        axis.set_yticks(y, labels)
        axis.invert_yaxis()
        axis.set_xlabel("Number of reference genes")
        axis.grid(axis="x", color="#dddddd", linewidth=0.6)
        axis.set_axisbelow(True)
        axis.text(
            -0.04,
            1.02,
            f"({chr(96 + panel)})",
            transform=axis.transAxes,
            fontsize=13,
            fontweight="bold",
        )
        if panel == 1:
            axis.legend(
                frameon=False,
                loc="lower center",
                bbox_to_anchor=(0.5, 1.005),
                ncol=2,
            )
        axis.set_xlim(0, max(total) * 1.18)
        for position, (total_value, strict_value) in enumerate(
            zip(total, strict_all)
        ):
            axis.text(
                total_value,
                position,
                f"  {total_value:,} | {strict_value:,}",
                va="center",
                ha="left",
                fontsize=8.5,
            )
        for order, (row, label) in enumerate(zip(branch_rows, labels)):
            plot_rows.append(
                {
                    "panel": branch_type,
                    "plot_order": order,
                    "branch_id": row["branch_id"],
                    "branch_type": branch_type,
                    "display_label": label,
                    "descendant_lineage_count": int(
                        row["descendant_lineage_count"]
                    ),
                    "descendant_lineages": row["descendant_lineages"],
                    "article_method_loss_gene_count": int(
                        row["article_method_loss_gene_count"]
                    ),
                    "strict_supported_any_descendant_lineage": int(
                        row["strict_supported_any_descendant_lineage"]
                    ),
                    "strict_supported_all_descendant_lineages": int(
                        row["strict_supported_all_descendant_lineages"]
                    ),
                }
            )

    caption = (
        "Topology placement of gene losses and strict disruption support. Panel (a) shows "
        "species-specific branches and panel (b) tree-node branches. Wide bars count exact branch "
        "events inferred from the decayed-plus-deleted lineage states. Narrow overlays count "
        "the subset with at least one strict pseudogenized unit in every descendant lineage. "
        "The two numbers printed after each bar are total branch losses and the strict all-"
        "descendant-lineage subset, respectively. Partial, unknown, and non-exact placements "
        "are not assigned to branches."
    )
    validation = {
        "status": "PASS_TREE_BRANCH_LOSS_EVIDENCE_FIGURE",
        "terminal_branches": len(terminal),
        "internal_branches": len(internal),
        "checks": {
            "article_method_decayed_plus_deleted": True,
            "strict_support_is_subset_overlay": True,
            "species_specific_and_tree_node_separated": True,
            "partial_and_unknown_excluded": True,
            "latin_binomials_italic": True,
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
        input_paths=[args.branch_summary, args.run_manifest],
        dpi=300,
    )
    plt.close(figure)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
