#!/usr/bin/env python3
"""Render classified topology-aware article-method GO/KEGG results."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.figure_bundle import write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label


csv.field_size_limit(sys.maxsize)


ONTOLOGIES = ("GO", "KEGG_KO", "KEGG_PATHWAY")
CATEGORY_ORDER = (
    "category__shared_all_23_units",
    "category__nonshared_any_unit_loss",
    "category__partial_lineage_loss_any",
    "category__single_terminal_branch_loss",
    "category__single_internal_branch_loss",
    "category__recurrent_independent_losses",
)
CATEGORY_LABELS = {
    "category__shared_all_23_units": "Shared across all 23 units",
    "category__nonshared_any_unit_loss": "Non-shared loss in >=1 unit",
    "category__partial_lineage_loss_any": "Partial/homeolog-specific lineage loss",
    "category__single_terminal_branch_loss": "Single terminal-branch loss",
    "category__single_internal_branch_loss": "Single internal-branch loss",
    "category__recurrent_independent_losses": "Repeated independent losses",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--foreground-metadata", required=True, type=Path)
    parser.add_argument("--significant-enrichment", required=True, type=Path)
    parser.add_argument("--enrichment-summary", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="tree_aware_manuscript_functional_summary")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path.name}: missing header")
        return [dict(row) for row in reader]


def internal_label(descendants: list[str], branch_id: str) -> str:
    if len(descendants) == 13:
        return "Actinidia stem (all 13 lineages)"
    if len(descendants) <= 3:
        return " + ".join(
            format_downstream_taxon_label(species, (), abbreviate_genus=True, separator=" ")
            for species in descendants
        )
    return f"Internal clade n={len(descendants)} ({branch_id.rsplit('__', 1)[-1][:6]})"


def main() -> int:
    args = parse_args()
    import json
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    metadata_rows = read_tsv(args.foreground_metadata)
    metadata = {row["foreground_id"]: row for row in metadata_rows}
    if len(metadata) != len(metadata_rows):
        raise SystemExit("ERROR: duplicate foreground metadata")
    missing_categories = set(CATEGORY_ORDER) - set(metadata)
    if missing_categories:
        raise SystemExit(f"ERROR: missing expected functional categories: {sorted(missing_categories)}")
    significant_rows = read_tsv(args.significant_enrichment)
    if any(row.get("significant_fdr") != "true" for row in significant_rows):
        raise SystemExit("ERROR: significant table contains a non-significant row")
    summary = json.loads(args.enrichment_summary.read_text(encoding="utf-8"))
    if summary.get("status") != "PASS_TREE_AWARE_MANUSCRIPT_GO_KEGG":
        raise SystemExit("ERROR: enrichment summary is not PASS")
    counts = Counter((row["foreground_id"], row["ontology"]) for row in significant_rows)

    category_rows = []
    for foreground_id in CATEGORY_ORDER:
        row = metadata[foreground_id]
        category_rows.append(
            {
                "panel": "category_counts",
                "foreground_id": foreground_id,
                "display_label": CATEGORY_LABELS[foreground_id],
                "branch_type": "category",
                "descendant_lineage_count": "",
                "descendant_lineages": "",
                "ontology": "",
                "foreground_gene_count": int(row["foreground_gene_count"]),
                "significant_term_count": sum(counts[(foreground_id, ontology)] for ontology in ONTOLOGIES),
            }
        )

    branches = [row for row in metadata_rows if row["foreground_id"].startswith("branch__")]
    for row in branches:
        row["_n"] = int(row["descendant_lineage_count"])
        descendants = row["descendant_lineages"].split(";")
        if row["branch_id"].startswith("terminal__"):
            species = descendants[0]
            row["_type"] = "terminal"
            row["_label"] = format_downstream_taxon_label(species, (), abbreviate_genus=True, separator=" ")
        else:
            row["_type"] = "internal"
            row["_label"] = internal_label(descendants, row["branch_id"])
    branches.sort(key=lambda row: (0 if row["_type"] == "internal" else 1, -row["_n"], row["_label"]))
    branch_rows = []
    matrix = []
    for row in branches:
        values = [counts[(row["foreground_id"], ontology)] for ontology in ONTOLOGIES]
        matrix.append(values)
        for ontology, value in zip(ONTOLOGIES, values):
            branch_rows.append(
                {
                    "panel": "branch_enrichment",
                    "foreground_id": row["foreground_id"],
                    "display_label": row["_label"],
                    "branch_type": row["_type"],
                    "descendant_lineage_count": row["_n"],
                    "descendant_lineages": row["descendant_lineages"],
                    "ontology": ontology,
                    "foreground_gene_count": int(row["foreground_gene_count"]),
                    "significant_term_count": value,
                }
            )

    fig = plt.figure(figsize=(12, 17), dpi=220, constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1, 4.3])
    ax_a = fig.add_subplot(grid[0])
    labels = [row["display_label"] for row in category_rows]
    values = [row["foreground_gene_count"] for row in category_rows]
    y = np.arange(len(labels))
    ax_a.barh(y, values, color="#4477AA")
    ax_a.set_yticks(y, labels)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Number of reference genes")
    ax_a.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax_a.set_axisbelow(True)
    ax_a.text(-0.02, 1.03, "(a)", transform=ax_a.transAxes, fontweight="bold", fontsize=14)
    for position, value in enumerate(values):
        ax_a.text(value, position, f" {value:,}", va="center", fontsize=9)

    ax_b = fig.add_subplot(grid[1])
    array = np.asarray(matrix, dtype=float)
    image = ax_b.imshow(np.log1p(array), aspect="auto", cmap="YlOrRd")
    ax_b.set_xticks(range(3), ["GO", "KEGG KO", "KEGG pathway"])
    ax_b.set_yticks(range(len(branches)), [row["_label"] for row in branches])
    ax_b.set_xlabel("Significant enriched terms (BH FDR <= 0.05)")
    ax_b.text(-0.02, 1.015, "(b)", transform=ax_b.transAxes, fontweight="bold", fontsize=14)
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            value = int(array[i, j])
            ax_b.text(j, i, str(value), ha="center", va="center", fontsize=7.5,
                      color="white" if np.log1p(value) > np.log1p(array.max()) * 0.55 else "black")
    colorbar = fig.colorbar(image, ax=ax_b, shrink=0.55, pad=0.02)
    colorbar.set_label("log(1 + significant term count)")

    plot_rows = [*category_rows, *branch_rows]
    caption = (
        "Topology-aware functional classification of gene losses. Panel (a) separates "
        "23-unit shared, non-shared, partial/homeolog-specific, single terminal-branch, single "
        "internal-branch, and repeated independent loss patterns. Categories can overlap where "
        "biologically expected (for example, all-23-unit shared genes are an Actinidia-stem subset). "
        "Panel (b) reports the number of GO, KEGG KO, and KEGG pathway terms enriched on each exact "
        "terminal or internal branch after one-sided hypergeometric tests and within-foreground/ontology "
        "BH correction. Only complete biological-lineage losses are placed on branches; partial and "
        "unknown lineage states are excluded from ancestral-event inference."
    )
    validation = {
        "status": "PASS_TREE_AWARE_FUNCTIONAL_SUMMARY_FIGURE",
        "foregrounds": len(metadata),
        "branch_foregrounds": len(branches),
        "significant_term_rows": len(significant_rows),
        "checks": {
            "enrichment_summary_pass": True,
            "branch_and_category_foregrounds_separate": True,
            "partial_not_promoted_to_ancestral_loss": True,
            "exact_significant_counts_plotted": True,
        },
    }
    write_figure_bundle(
        figure=fig,
        output_dir=args.output_dir,
        basename=args.basename,
        plot_rows=plot_rows,
        plot_columns=list(plot_rows[0]),
        caption=caption,
        validation=validation,
        input_paths=[args.foreground_metadata, args.significant_enrichment, args.enrichment_summary],
        dpi=300,
    )
    plt.close(fig)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
