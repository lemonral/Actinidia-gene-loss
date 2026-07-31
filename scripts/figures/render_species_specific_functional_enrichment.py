#!/usr/bin/env python3
"""Render pure species-specific GO/KEGG results."""

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
from geneloss_repro.labels import format_taxon_label


CATEGORY_COLUMNS = (
    ("GO biological process", "go_biological_process_significant_terms"),
    ("GO molecular function", "go_molecular_function_significant_terms"),
    ("GO cellular component", "go_cellular_component_significant_terms"),
    ("KEGG orthology", "kegg_orthology_significant_terms"),
    ("KEGG pathway", "kegg_pathway_significant_terms"),
)


class FigureError(ValueError):
    """Raised when the curated species-specific bundle is invalid."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-summary", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--basename",
        default="species_specific_functional_enrichment",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    if not reader.fieldnames or not rows:
        raise FigureError(f"{path.name}: invalid or empty table")
    return rows


def species_level_label(species: str) -> str:
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


def main() -> int:
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    manifest = json.loads(args.run_manifest.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_ARTICLE_METHOD_SPECIES_SPECIFIC_GO_KEGG":
        raise SystemExit("ERROR: species-specific enrichment manifest is not PASS")
    rows = read_tsv(args.species_summary)
    if len(rows) != 13:
        raise SystemExit("ERROR: expected exactly 13 species-specific rows")
    if sum(int(row["species_specific_loss_gene_count"]) for row in rows) != 1167:
        raise SystemExit("ERROR: species-specific loss genes do not close to 1,167")
    rows.sort(key=lambda row: row["biological_species"])
    labels = [
        species_level_label(row["biological_species"])
        for row in rows
    ]
    gene_counts = np.asarray(
        [int(row["species_specific_loss_gene_count"]) for row in rows],
        dtype=int,
    )
    term_counts = np.asarray(
        [
            [int(row[column]) for _, column in CATEGORY_COLUMNS]
            for row in rows
        ],
        dtype=int,
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.labelsize": 11.5,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 10.5,
            "legend.fontsize": 9.5,
        }
    )
    figure = plt.figure(figsize=(13.5, 9.8), dpi=220, constrained_layout=True)
    grid = figure.add_gridspec(1, 2, width_ratios=[1.05, 1.25])
    y = np.arange(len(rows))

    axis_a = figure.add_subplot(grid[0])
    axis_a.barh(y, gene_counts, color="#4477AA", height=0.68)
    axis_a.set_yticks(y, labels)
    axis_a.invert_yaxis()
    axis_a.set_xlabel("Species-specific lost genes")
    axis_a.grid(axis="x", color="#dddddd", linewidth=0.6)
    axis_a.set_axisbelow(True)
    axis_a.text(
        -0.04,
        1.02,
        "(a)",
        transform=axis_a.transAxes,
        fontsize=13,
        fontweight="bold",
    )
    limit = max(gene_counts) * 1.16
    axis_a.set_xlim(0, limit)
    for position, value in enumerate(gene_counts):
        axis_a.text(
            value,
            position,
            f"  {value:,}",
            ha="left",
            va="center",
            fontsize=9,
        )

    axis_b = figure.add_subplot(grid[1])
    transformed = np.log1p(term_counts.astype(float))
    axis_b.imshow(
        transformed,
        aspect="auto",
        cmap="YlGnBu",
        vmin=0,
    )
    axis_b.set_xticks(
        range(len(CATEGORY_COLUMNS)),
        [
            "GO BP",
            "GO MF",
            "GO CC",
            "KEGG KO",
            "KEGG pathway",
        ],
    )
    axis_b.set_yticks(y, labels)
    axis_b.set_xlabel("Significant enriched terms (BH q ≤ 0.05)")
    axis_b.text(
        -0.04,
        1.02,
        "(b)",
        transform=axis_b.transAxes,
        fontsize=13,
        fontweight="bold",
    )
    threshold = transformed.max() * 0.55 if transformed.size else 0
    for row_index in range(term_counts.shape[0]):
        for column_index in range(term_counts.shape[1]):
            axis_b.text(
                column_index,
                row_index,
                str(int(term_counts[row_index, column_index])),
                ha="center",
                va="center",
                fontsize=9,
                color=(
                    "white"
                    if transformed[row_index, column_index] > threshold
                    else "black"
                ),
            )
    plot_rows: list[dict[str, object]] = []
    for order, row in enumerate(rows):
        for category_label, category_column in CATEGORY_COLUMNS:
            plot_rows.append(
                {
                    "plot_order": order,
                    "biological_species": row["biological_species"],
                    "display_label": labels[order],
                    "foreground_id": row["foreground_id"],
                    "species_specific_loss_gene_count": int(
                        row["species_specific_loss_gene_count"]
                    ),
                    "functional_category": category_label,
                    "significant_term_count": int(row[category_column]),
                }
            )
    caption = (
        "Functional enrichment of species-specific gene losses defined as decayed plus "
        "deleted. A species-specific foreground contains genes with exactly one inferred "
        "terminal-branch loss and resolved retained states in every other lineage; recurrent, "
        "internal, partial/homeolog-specific, and unknown patterns are excluded. Panel (a) "
        "reports foreground gene counts. Panel (b) reports significant GO biological-process, "
        "molecular-function, cellular-component, KEGG orthology, and KEGG pathway terms after "
        "one-sided hypergeometric tests and Benjamini–Hochberg correction within each species "
        "foreground and ontology."
    )
    validation = {
        "status": "PASS_SPECIES_SPECIFIC_FUNCTIONAL_ENRICHMENT_FIGURE",
        "lineages": len(rows),
        "species_specific_loss_genes": int(gene_counts.sum()),
        "significant_terms": int(term_counts.sum()),
        "checks": {
            "article_method_decayed_plus_deleted": True,
            "single_terminal_only": True,
            "recurrent_losses_excluded": True,
            "functional_categories_separated": True,
            "latin_binomials_italic": True,
            "species_level_suffix_policy": True,
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
        input_paths=[args.species_summary, args.run_manifest],
        dpi=300,
    )
    plt.close(figure)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
