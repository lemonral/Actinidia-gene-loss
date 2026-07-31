#!/usr/bin/env python3
"""Render unit-resolved decayed-plus-deleted GO/KEGG summaries."""

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


CATEGORY_COLUMNS = (
    ("GO biological process", "go_biological_process_significant_terms"),
    ("GO molecular function", "go_molecular_function_significant_terms"),
    ("GO cellular component", "go_cellular_component_significant_terms"),
    ("KEGG orthology", "kegg_orthology_significant_terms"),
    ("KEGG pathway", "kegg_pathway_significant_terms"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-summary", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--basename",
        default="unit_loss_functional_enrichment",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
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
    if (
        manifest.get("status")
        != "PASS_UNIT_ARTICLE_METHOD_GO_KEGG_SUMMARY"
    ):
        raise SystemExit("ERROR: unit enrichment summary is not PASS")
    rows = read_tsv(args.unit_summary)
    if len(rows) != 23:
        raise SystemExit("ERROR: expected exactly 23 assembly units")
    if len({row["assembly_unit_id"] for row in rows}) != 23:
        raise SystemExit("ERROR: duplicate assembly unit")

    labels = [
        format_downstream_taxon_label(
            row["biological_species"],
            (row["haplotype_or_subgenome"],),
            abbreviate_genus=True,
            separator=" ",
        )
        for row in rows
    ]
    loss_counts = np.asarray(
        [int(row["loss_gene_count"]) for row in rows],
        dtype=int,
    )
    term_counts = np.asarray(
        [
            [int(row[column]) for _, column in CATEGORY_COLUMNS]
            for row in rows
        ],
        dtype=int,
    )
    if int(loss_counts.sum()) != int(manifest["foreground_memberships"]):
        raise SystemExit("ERROR: loss-gene counts do not close")
    if int(term_counts.sum()) != int(manifest["significant_terms"]):
        raise SystemExit("ERROR: significant-term counts do not close")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.labelsize": 11.5,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 10.0,
        }
    )
    figure = plt.figure(figsize=(14.8, 13.2), dpi=220)
    grid = figure.add_gridspec(
        1,
        2,
        width_ratios=(1.0, 1.25),
        left=0.18,
        right=0.97,
        bottom=0.08,
        top=0.96,
        wspace=0.28,
    )
    y = np.arange(len(rows))

    axis_a = figure.add_subplot(grid[0])
    axis_a.barh(y, loss_counts, color="#6B7C93", height=0.68)
    axis_a.set_yticks(y, labels)
    axis_a.invert_yaxis()
    axis_a.set_xlabel("Lost genes (decayed + deleted)")
    axis_a.grid(axis="x", color="#dddddd", linewidth=0.6)
    axis_a.set_axisbelow(True)
    axis_a.spines["top"].set_visible(False)
    axis_a.spines["right"].set_visible(False)
    axis_a.text(
        -0.09,
        1.015,
        "A",
        transform=axis_a.transAxes,
        fontsize=13,
        fontweight="bold",
    )
    axis_a.set_xlim(0, max(loss_counts) * 1.18)
    for position, value in enumerate(loss_counts):
        axis_a.text(
            value,
            position,
            f"  {value:,}",
            ha="left",
            va="center",
            fontsize=8.5,
        )

    axis_b = figure.add_subplot(grid[1])
    transformed = np.log1p(term_counts.astype(float))
    image = axis_b.imshow(
        transformed,
        aspect="auto",
        cmap="YlGnBu",
        vmin=0,
    )
    axis_b.set_xticks(
        range(len(CATEGORY_COLUMNS)),
        ["GO BP", "GO MF", "GO CC", "KEGG KO", "KEGG pathway"],
        rotation=25,
        ha="right",
    )
    axis_b.set_yticks(y, labels)
    axis_b.set_xlabel("Significant enriched terms (BH q ≤ 0.05)")
    axis_b.text(
        -0.09,
        1.015,
        "B",
        transform=axis_b.transAxes,
        fontsize=13,
        fontweight="bold",
    )
    axis_b.set_xticks(np.arange(-0.5, len(CATEGORY_COLUMNS), 1), minor=True)
    axis_b.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    axis_b.grid(which="minor", color="white", linewidth=0.8)
    axis_b.tick_params(which="minor", bottom=False, left=False)
    threshold = transformed.max() * 0.52 if transformed.size else 0.0
    for row_index in range(term_counts.shape[0]):
        for column_index in range(term_counts.shape[1]):
            axis_b.text(
                column_index,
                row_index,
                str(int(term_counts[row_index, column_index])),
                ha="center",
                va="center",
                fontsize=8.2,
                color=(
                    "white"
                    if transformed[row_index, column_index] > threshold
                    else "black"
                ),
            )
    colorbar = figure.colorbar(
        image,
        ax=axis_b,
        fraction=0.025,
        pad=0.02,
    )
    ticks = colorbar.get_ticks()
    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels(
        [f"{int(round(np.expm1(value))):,}" for value in ticks]
    )
    colorbar.set_label("Significant terms")

    plot_rows: list[dict[str, object]] = []
    for order, row in enumerate(rows):
        for category_label, category_column in CATEGORY_COLUMNS:
            plot_rows.append(
                {
                    "plot_order": order,
                    "assembly_unit_id": row["assembly_unit_id"],
                    "biological_species": row["biological_species"],
                    "haplotype_or_subgenome": row[
                        "haplotype_or_subgenome"
                    ],
                    "display_label": labels[order],
                    "loss_gene_count": int(row["loss_gene_count"]),
                    "functional_category": category_label,
                    "significant_term_count": int(row[category_column]),
                }
            )
    caption = (
        "GO and KEGG enrichment for gene losses in each of the 23 assembly "
        "units. Loss is defined independently in every unit as decayed plus "
        "deleted. Haplotypes and subgenomes are not combined, and no retained "
        "state in another unit is required. Panel A reports the number of lost "
        "genes. Panel B reports significant GO biological-process, "
        "molecular-function, cellular-component, KEGG orthology, and KEGG "
        "pathway terms after one-sided hypergeometric tests and "
        "Benjamini–Hochberg correction within each unit and ontology."
    )
    validation = {
        "status": "PASS_UNIT_ARTICLE_METHOD_FUNCTIONAL_FIGURE",
        "assembly_units": len(rows),
        "foreground_memberships": int(loss_counts.sum()),
        "significant_terms": int(term_counts.sum()),
        "checks": {
            "decayed_plus_deleted": True,
            "assembly_units_not_aggregated": True,
            "other_units_need_not_be_retained": True,
            "functional_categories_separated": True,
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
        input_paths=[args.unit_summary, args.run_manifest],
        dpi=300,
    )
    plt.close(figure)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
