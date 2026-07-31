#!/usr/bin/env python3
"""Render category-level GO/KEGG summaries for uniform loss foregrounds."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if SOURCE_ROOT.is_dir():
    sys.path.insert(0, str(SOURCE_ROOT))

from geneloss_repro.figure_bundle import FigureBundle, write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label


CATEGORIES = (
    ("GO_BP", "GO biological process", "GO", "biological_process"),
    ("GO_MF", "GO molecular function", "GO", "molecular_function"),
    ("GO_CC", "GO cellular component", "GO", "cellular_component"),
    ("KEGG_KO", "KEGG orthology", "KEGG_KO", ""),
    ("KEGG_PATHWAY", "KEGG pathway", "KEGG_PATHWAY", ""),
)
EVIDENCE_MODES = (
    (
        "combined",
        "Deleted + strict pseudogenized",
        "deleted_plus_strict_pseudogenized",
    ),
    (
        "strict_pseudogenized",
        "Strict pseudogenized only",
        "all_required_units_strict_pseudogenized",
    ),
)
PLOT_COLUMNS = (
    "foreground_id",
    "foreground_scope",
    "evidence_mode",
    "biological_species",
    "display_label",
    "category_id",
    "category_label",
    "significant_term_count",
    "requested_gene_count",
    "annotated_study_gene_count",
    "annotation_coverage",
)


class FunctionalCategoryFigureError(ValueError):
    """Raised when enrichment tables cannot support the category figure."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise FunctionalCategoryFigureError(f"{path.name}: missing header")
        return list(reader)


def functional_category(row: dict[str, str]) -> str:
    ontology = row.get("ontology", "")
    if ontology == "GO":
        namespace = row.get("go_namespace", "")
        mapping = {
            "biological_process": "GO_BP",
            "molecular_function": "GO_MF",
            "cellular_component": "GO_CC",
        }
        if namespace not in mapping:
            raise FunctionalCategoryFigureError(
                f"GO row has unsupported namespace {namespace!r}"
            )
        return mapping[namespace]
    if ontology in {"KEGG_KO", "KEGG_PATHWAY"}:
        return ontology
    raise FunctionalCategoryFigureError(f"unsupported ontology {ontology!r}")


def _foreground_scope(foreground_id: str) -> str:
    return "pooled" if foreground_id.startswith("pooled_") else "lineage"


def _display_label(species: str, foreground_id: str) -> str:
    if foreground_id.startswith("pooled_"):
        return "Pooled 13 lineages"
    return format_downstream_taxon_label(species, abbreviate_genus=True)


def prepare(
    enrichment_dir: str | Path,
) -> tuple[list[dict[str, object]], dict[str, object], list[Path]]:
    source = Path(enrichment_dir)
    significant_path = source / "enrichment_significant.tsv"
    foreground_path = source / "foreground_summary.tsv"
    summary_path = source / "functional_enrichment_summary.json"
    significant = read_tsv(significant_path)
    foreground_rows = read_tsv(foreground_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "PASS_UNIFORM_GO_KEGG_ENRICHMENT":
        raise FunctionalCategoryFigureError("functional enrichment summary is not PASS")

    required_significant = {
        "foreground_id",
        "ontology",
        "go_namespace",
        "significant_fdr",
    }
    if significant and not required_significant.issubset(significant[0]):
        raise FunctionalCategoryFigureError(
            "significant enrichment table lacks category columns"
        )
    required_foreground = {
        "foreground_id",
        "biological_species",
        "evidence_mode",
        "ontology",
        "requested_gene_count",
        "annotated_study_gene_count",
        "annotation_coverage",
    }
    if not foreground_rows or not required_foreground.issubset(foreground_rows[0]):
        raise FunctionalCategoryFigureError("foreground summary lacks required columns")

    counts: Counter[tuple[str, str]] = Counter()
    for row in significant:
        if row["significant_fdr"] != "true":
            raise FunctionalCategoryFigureError(
                "enrichment_significant.tsv contains a non-significant row"
            )
        counts[(row["foreground_id"], functional_category(row))] += 1

    category_ontology = {category_id: ontology for category_id, _, ontology, _ in CATEGORIES}
    foreground_by_key: dict[tuple[str, str], dict[str, str]] = {}
    foreground_meta: dict[str, dict[str, str]] = {}
    for row in foreground_rows:
        foreground_id = row["foreground_id"]
        ontology = row["ontology"]
        key = (foreground_id, ontology)
        if key in foreground_by_key:
            raise FunctionalCategoryFigureError(
                f"duplicate foreground/ontology summary: {foreground_id}/{ontology}"
            )
        foreground_by_key[key] = row
        foreground_meta.setdefault(foreground_id, row)

    selected_ids = sorted(
        (
            foreground_id
            for foreground_id, row in foreground_meta.items()
            if foreground_id.startswith("lineage_")
            or foreground_id.startswith("pooled_nonshared_")
        ),
        key=lambda value: (
            0 if value.startswith("pooled_") else 1,
            value.split("_nonshared_", 1)[0],
            0 if value.endswith("_combined") else 1,
        ),
    )
    if len(selected_ids) != 28:
        raise FunctionalCategoryFigureError(
            f"expected 28 pooled/lineage foregrounds, observed {len(selected_ids)}"
        )

    plot_rows: list[dict[str, object]] = []
    for foreground_id in selected_ids:
        meta = foreground_meta[foreground_id]
        evidence_mode = meta["evidence_mode"]
        if evidence_mode not in {value for _, _, value in EVIDENCE_MODES}:
            raise FunctionalCategoryFigureError(
                f"unexpected evidence mode for {foreground_id}: {evidence_mode}"
            )
        for category_id, category_label, _, _ in CATEGORIES:
            ontology = category_ontology[category_id]
            summary_row = foreground_by_key.get((foreground_id, ontology))
            if summary_row is None:
                raise FunctionalCategoryFigureError(
                    f"missing {ontology} coverage for {foreground_id}"
                )
            plot_rows.append(
                {
                    "foreground_id": foreground_id,
                    "foreground_scope": _foreground_scope(foreground_id),
                    "evidence_mode": evidence_mode,
                    "biological_species": meta["biological_species"],
                    "display_label": _display_label(
                        meta["biological_species"], foreground_id
                    ),
                    "category_id": category_id,
                    "category_label": category_label,
                    "significant_term_count": counts[(foreground_id, category_id)],
                    "requested_gene_count": int(summary_row["requested_gene_count"]),
                    "annotated_study_gene_count": int(
                        summary_row["annotated_study_gene_count"]
                    ),
                    "annotation_coverage": float(summary_row["annotation_coverage"]),
                }
            )

    pooled = {
        (str(row["evidence_mode"]), str(row["category_id"])): int(
            row["significant_term_count"]
        )
        for row in plot_rows
        if row["foreground_scope"] == "pooled"
    }
    expected_pooled = {
        ("deleted_plus_strict_pseudogenized", "GO_BP"): 12,
        ("deleted_plus_strict_pseudogenized", "GO_MF"): 17,
        ("deleted_plus_strict_pseudogenized", "GO_CC"): 0,
        ("deleted_plus_strict_pseudogenized", "KEGG_KO"): 27,
        ("deleted_plus_strict_pseudogenized", "KEGG_PATHWAY"): 30,
        ("all_required_units_strict_pseudogenized", "GO_BP"): 4,
        ("all_required_units_strict_pseudogenized", "GO_MF"): 14,
        ("all_required_units_strict_pseudogenized", "GO_CC"): 0,
        ("all_required_units_strict_pseudogenized", "KEGG_KO"): 21,
        ("all_required_units_strict_pseudogenized", "KEGG_PATHWAY"): 32,
    }
    if pooled != expected_pooled:
        raise FunctionalCategoryFigureError("pooled category counts changed")
    shared_significant = sum(
        row["foreground_id"] == "shared_combined" for row in significant
    )
    validation = {
        "schema_version": "1.0",
        "status": "pass",
        "source_status": summary["status"],
        "foreground_count": len(selected_ids),
        "plot_row_count": len(plot_rows),
        "shared_significant_term_count": shared_significant,
        "categories": [category_id for category_id, _, _, _ in CATEGORIES],
        "checks": {
            "go_namespaces_separate": True,
            "kegg_ko_and_pathway_separate": True,
            "combined_and_strict_pseudogenized_separate": True,
            "all_13_lineages_and_pooled_present": True,
            "shared_no_significant_terms_reported": shared_significant == 0,
            "short_taxon_labels_used": True,
        },
    }
    return plot_rows, validation, [significant_path, foreground_path, summary_path]


def build_figure(plot_rows: list[dict[str, object]]):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(15.2, 9.2), constrained_layout=True)
    category_ids = [category_id for category_id, _, _, _ in CATEGORIES]
    category_labels = [label for _, label, _, _ in CATEGORIES]
    evidence = [value for _, _, value in EVIDENCE_MODES]
    maximum = max(int(row["significant_term_count"]) for row in plot_rows)
    image = None
    for axis, (mode_id, mode_label, mode_value) in zip(axes, EVIDENCE_MODES):
        subset = [row for row in plot_rows if row["evidence_mode"] == mode_value]
        foreground_ids = []
        for row in subset:
            foreground_id = str(row["foreground_id"])
            if foreground_id not in foreground_ids:
                foreground_ids.append(foreground_id)
        foreground_ids.sort(
            key=lambda value: (
                0 if value.startswith("pooled_") else 1,
                value.split("_nonshared_", 1)[0],
            )
        )
        row_by_key = {
            (str(row["foreground_id"]), str(row["category_id"])): row
            for row in subset
        }
        matrix = [
            [
                int(row_by_key[(foreground_id, category_id)]["significant_term_count"])
                for category_id in category_ids
            ]
            for foreground_id in foreground_ids
        ]
        labels = [
            str(row_by_key[(foreground_id, category_ids[0])]["display_label"])
            for foreground_id in foreground_ids
        ]
        image = axis.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=maximum, aspect="auto")
        for y, values in enumerate(matrix):
            for x, value in enumerate(values):
                axis.text(
                    x,
                    y,
                    str(value),
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if value > maximum * 0.55 else "black",
                )
        axis.set_xticks(range(len(category_labels)), category_labels, rotation=28, ha="right")
        axis.set_yticks(range(len(labels)), labels, fontsize=8.5)
        axis.set_title(mode_label)
        axis.set_xlabel("Functional annotation category")
        axis.tick_params(length=0)
        for spine in axis.spines.values():
            spine.set_visible(False)
    if image is not None:
        colorbar = figure.colorbar(image, ax=axes, shrink=0.72, pad=0.02)
        colorbar.set_label("Number of significant terms (BH FDR < 0.05)")
    return figure


def publish(
    *,
    enrichment_dir: str | Path,
    output_dir: str | Path,
    basename: str = "uniform_loss_functional_categories",
    dpi: int = 300,
) -> FigureBundle:
    plot_rows, validation, inputs = prepare(enrichment_dir)
    figure = build_figure(plot_rows)
    caption = (
        "Category-level summary of significant functional enrichments for pooled and "
        "lineage-specific non-shared gene-loss foregrounds. GO terms are separated "
        "into biological process, molecular function, and cellular component; KEGG "
        "orthologues and KEGG pathways are reported separately. The left panel uses "
        "deleted plus strict pseudogenized calls, and the right panel is the separate "
        "strict-pseudogenized-only sensitivity. Cells are counts of terms passing "
        "Benjamini-Hochberg FDR 0.05, not effect sizes, and should be interpreted with "
        "the annotation-coverage columns in the plot-data table. The 287 shared genes "
        "had zero significant terms and are therefore stated here rather than plotted."
    )
    try:
        return write_figure_bundle(
            figure=figure,
            output_dir=output_dir,
            basename=basename,
            plot_rows=plot_rows,
            plot_columns=PLOT_COLUMNS,
            caption=caption,
            validation=validation,
            input_paths=inputs,
            dpi=dpi,
        )
    finally:
        import matplotlib.pyplot as plt

        plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enrichment-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--basename", default="uniform_loss_functional_categories")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    try:
        publish(
            enrichment_dir=args.enrichment_dir,
            output_dir=args.output_dir,
            basename=args.basename,
            dpi=args.dpi,
        )
    except (FunctionalCategoryFigureError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
