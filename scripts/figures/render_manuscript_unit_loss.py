#!/usr/bin/env python3
"""Render 23-unit article-method shared/non-shared loss composition."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.figure_bundle import write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label_from_metadata


PUBLICATION_CATEGORY_LABELS = {
    "shared_decayed": "Shared decayed",
    "shared_deleted": "Shared deleted",
    "nonshared_decayed": "Non-shared decayed",
    "nonshared_deleted": "Non-shared deleted",
}
PUBLICATION_CATEGORY_COLORS = {
    "shared_decayed": "#56B4E9",
    "shared_deleted": "#595959",
    "nonshared_decayed": "#E69F00",
    "nonshared_deleted": "#B3B3B3",
}
PUBLICATION_VALUE_COLUMN = "gene_count"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript-matrix", required=True, type=Path)
    parser.add_argument("--shared-genes", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="manuscript_method_unit_loss")
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matrix = pd.read_csv(
        args.manuscript_matrix,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        compression="infer",
        usecols=["reference_gene_id", "assembly_unit_id", "manuscript_classification"],
    )
    if len(matrix) != args.expected_units * args.expected_reference_genes:
        raise SystemExit("ERROR: manuscript matrix row count does not close")
    if matrix.duplicated(["reference_gene_id", "assembly_unit_id"]).any():
        raise SystemExit("ERROR: duplicate gene-unit rows")
    shared_table = pd.read_csv(args.shared_genes, sep="\t", dtype=str, keep_default_na=False)
    if "reference_gene_id" not in shared_table or shared_table["reference_gene_id"].duplicated().any():
        raise SystemExit("ERROR: invalid shared-gene table")
    shared = set(shared_table["reference_gene_id"])
    if not shared:
        raise SystemExit("ERROR: shared-gene table is empty")
    metadata = pd.read_csv(args.unit_metadata, sep="\t", dtype=str, keep_default_na=False)
    required_metadata = {"assembly_unit_id", "biological_species", "haplotype_or_subgenome", "include"}
    if not required_metadata.issubset(metadata.columns):
        raise SystemExit("ERROR: unit metadata is missing required columns")
    metadata = metadata.loc[metadata["include"].str.lower() == "true"].copy()
    if len(metadata) != args.expected_units or metadata["assembly_unit_id"].duplicated().any():
        raise SystemExit("ERROR: expected exact unique 23-unit metadata")
    if set(metadata["assembly_unit_id"]) != set(matrix["assembly_unit_id"]):
        raise SystemExit("ERROR: unit metadata and manuscript matrix disagree")

    matrix["shared"] = matrix["reference_gene_id"].isin(shared)
    matrix["category"] = ""
    matrix.loc[matrix["manuscript_classification"] == "retained", "category"] = "retained"
    matrix.loc[matrix["manuscript_classification"] == "not_called_loss", "category"] = "not_called"
    matrix.loc[
        matrix["shared"] & (matrix["manuscript_classification"] == "decayed"),
        "category",
    ] = "shared_decayed"
    matrix.loc[
        matrix["shared"] & (matrix["manuscript_classification"] == "deleted"),
        "category",
    ] = "shared_deleted"
    matrix.loc[
        ~matrix["shared"] & (matrix["manuscript_classification"] == "decayed"), "category"
    ] = "nonshared_decayed"
    matrix.loc[
        ~matrix["shared"] & (matrix["manuscript_classification"] == "deleted"), "category"
    ] = "nonshared_deleted"
    if (matrix["category"] == "").any():
        raise SystemExit("ERROR: an article-method row was not categorized")
    counts = matrix.groupby(["assembly_unit_id", "category"]).size().unstack(fill_value=0)
    categories = [
        "shared_decayed", "shared_deleted", "nonshared_decayed", "nonshared_deleted",
        "retained", "not_called"
    ]
    counts = counts.reindex(columns=categories, fill_value=0).reset_index()
    counts = counts.merge(metadata, on="assembly_unit_id", validate="one_to_one")
    counts["resolved_denominator"] = (
        counts["shared_decayed"] + counts["shared_deleted"]
        + counts["nonshared_decayed"] + counts["nonshared_deleted"] + counts["retained"]
    )
    if (counts["resolved_denominator"] + counts["not_called"] != args.expected_reference_genes).any():
        raise SystemExit("ERROR: per-unit reference denominator does not close")
    for category in PUBLICATION_CATEGORY_LABELS:
        counts[f"{category}_rate"] = counts[category] / counts["resolved_denominator"]
    counts["total_positive_rate"] = (
        counts["shared_decayed_rate"] + counts["shared_deleted_rate"]
        + counts["nonshared_decayed_rate"] + counts["nonshared_deleted_rate"]
    )
    counts = counts.sort_values(["biological_species", "haplotype_or_subgenome", "assembly_unit_id"])
    counts["display_label"] = [
        format_downstream_taxon_label_from_metadata(
            row,
            suffix_fields=("haplotype_or_subgenome",),
            abbreviate_genus=True,
            separator=" ",
        )
        for row in counts.to_dict("records")
    ]

    plot_rows = []
    labels = PUBLICATION_CATEGORY_LABELS
    for order, row in enumerate(counts.to_dict("records")):
        for category in labels:
            plot_rows.append(
                {
                    "assembly_unit_id": row["assembly_unit_id"],
                    "biological_species": row["biological_species"],
                    "haplotype_or_subgenome": row["haplotype_or_subgenome"],
                    "display_label": row["display_label"],
                    "plot_order": order,
                    "category": category,
                    "category_label": labels[category],
                    "gene_count": int(row[category]),
                    "resolved_denominator": int(row["resolved_denominator"]),
                    "loss_rate": float(row[f"{category}_rate"]),
                    "total_positive_loss_rate": float(row["total_positive_rate"]),
                    "not_called_count": int(row["not_called"]),
                }
            )

    fig, ax = plt.subplots(figsize=(10.5, 12.5), dpi=220, constrained_layout=True)
    colors = PUBLICATION_CATEGORY_COLORS
    y = list(range(len(counts)))
    left = [0] * len(counts)
    for category in labels:
        values = counts[category].astype(int).tolist()
        ax.barh(y, values, left=left, color=colors[category], label=labels[category], height=0.72)
        left = [a + b for a, b in zip(left, values)]
    for position, total in zip(y, left):
        ax.text(total, position, f"  {total:,}", va="center", ha="left", fontsize=8)
    ax.set_yticks(y, counts["display_label"].tolist())
    ax.invert_yaxis()
    ax.set_xlabel("Number of genes")
    ax.set_xlim(0, max(left) * 1.10)
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax.set_axisbelow(True)
    caption = (
        "Gene-loss composition across 23 assembly units. Bars report gene counts. Shared-loss "
        "genes are positive in all 23 units. Shared and non-shared calls are each divided into "
        "decayed and deleted evidence classes; both deleted classes are grey. Exact SynOrths "
        "anchors are retained and not plotted. Rows without a loss "
        "classification are also excluded. Latin binomials are italic and assembly-unit suffixes "
        "are upright."
    )
    validation = {
        "status": "PASS_MANUSCRIPT_METHOD_UNIT_LOSS_FIGURE",
        "assembly_units": len(counts),
        "reference_genes": args.expected_reference_genes,
        "shared_positive_genes": len(shared),
        "checks": {
            "exact_unit_universe": True,
            "per_unit_denominators_close": True,
            "article_classes_not_rewritten_by_refined_causes": True,
            "display_uses_gene_counts": True,
            "deleted_calls_are_grey": True,
            "publication_labels_omit_internal_threshold_wording": True,
            "downstream_labels_simplified": True,
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
        input_paths=[args.manuscript_matrix, args.shared_genes, args.unit_metadata],
        dpi=300,
    )
    plt.close(fig)
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
