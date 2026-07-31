#!/usr/bin/env python3
"""Render representative GO/KEGG terms for 23 unit-resolved loss sets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import statistics
import sys
import tempfile
import textwrap
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping


ROOT = Path(__file__).resolve().parents[2]
if (ROOT / "src").is_dir():
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.figure_bundle import write_figure_bundle
from geneloss_repro.labels import format_downstream_taxon_label


CATEGORY_LIMITS = (
    ("GO biological process", 10, "A", "bp"),
    ("GO molecular function", 8, "B", "mf"),
    ("GO cellular component", 6, "C", "cc"),
    ("KEGG pathway", 10, "D", "kegg"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--significant-enrichment", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--kegg-pathway-names", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--basename",
        default="unit_loss_functional",
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


def read_pathway_names(path: Path) -> dict[str, str]:
    rows, fields = read_tsv(path)
    if fields != ["term_id", "term_name"]:
        raise ValueError(f"{path.name}: unexpected pathway-name header")
    names = {row["term_id"]: row["term_name"] for row in rows}
    if len(names) != len(rows) or any(
        not term_id.startswith("map") or not name
        for term_id, name in names.items()
    ):
        raise ValueError(f"{path.name}: invalid pathway-name mapping")
    return names


def read_go_parents(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    parents: dict[str, set[str]] = defaultdict(set)
    aliases: dict[str, str] = {}
    current_id = ""
    current_aliases: list[str] = []
    obsolete = False

    def finish_term() -> None:
        if current_id and not obsolete:
            for alias in current_aliases:
                aliases[alias] = current_id

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                finish_term()
                current_id = ""
                current_aliases = []
                obsolete = False
            elif line.startswith("["):
                finish_term()
                current_id = ""
                current_aliases = []
                obsolete = False
            elif line.startswith("id: GO:"):
                current_id = line.removeprefix("id: ").strip()
            elif line.startswith("alt_id: GO:"):
                current_aliases.append(line.removeprefix("alt_id: ").strip())
            elif line.startswith("is_a: GO:") and current_id:
                parents[current_id].add(line.split()[1])
            elif line == "is_obsolete: true":
                obsolete = True
    finish_term()
    return dict(parents), aliases


def ancestor_checker(
    parents: Mapping[str, set[str]],
    aliases: Mapping[str, str],
):
    cache: dict[str, frozenset[str]] = {}

    def ancestors(term_id: str) -> frozenset[str]:
        canonical = aliases.get(term_id, term_id)
        if canonical in cache:
            return cache[canonical]
        found: set[str] = set()
        stack = list(parents.get(canonical, set()))
        while stack:
            candidate = stack.pop()
            parent = aliases.get(candidate, candidate)
            if not parent or parent in found:
                continue
            found.add(parent)
            stack.extend(parents.get(parent, set()))
        cache[canonical] = frozenset(found)
        return cache[canonical]

    def related(first: str, second: str) -> bool:
        first = aliases.get(first, first)
        second = aliases.get(second, second)
        return first in ancestors(second) or second in ancestors(first)

    return related


def functional_category(row: Mapping[str, str]) -> str:
    if row["ontology"] == "GO":
        return {
            "biological_process": "GO biological process",
            "molecular_function": "GO molecular function",
            "cellular_component": "GO cellular component",
        }[row["go_namespace"]]
    if row["ontology"] == "KEGG_PATHWAY":
        return "KEGG pathway"
    if row["ontology"] == "KEGG_KO":
        return "KEGG orthology"
    raise ValueError(f"unsupported ontology {row['ontology']!r}")


def rank_terms(
    rows: Iterable[Mapping[str, object]],
    *,
    limit: int,
    sample_key: str = "assembly_unit_id",
    go_related=None,
    jaccard_cutoff: float = 0.70,
) -> list[tuple[str, str]]:
    grouped: dict[tuple[str, str], list[Mapping[str, object]]] = defaultdict(
        list
    )
    for row in rows:
        grouped[(str(row["term_id"]), str(row["term_name"]))].append(row)
    ranked: list[
        tuple[tuple[float, ...], tuple[str, str], frozenset[str]]
    ] = []
    for term, members in grouped.items():
        unit_count = len({str(row[sample_key]) for row in members})
        median_score = statistics.median(
            float(row["minus_log10_q"]) for row in members
        )
        median_fold = statistics.median(
            float(row["fold_enrichment"]) for row in members
        )
        total_genes = sum(int(row["study_count"]) for row in members)
        genes = frozenset(
            gene
            for row in members
            for gene in str(row["study_gene_ids"]).split(";")
            if gene
        )
        ranked.append(
            (
                (
                    -float(unit_count),
                    -median_score,
                    -median_fold,
                    -float(total_genes),
                ),
                term,
                genes,
            )
        )
    ranked.sort(key=lambda item: (item[0], item[1][0], item[1][1]))
    selected: list[tuple[tuple[str, str], frozenset[str]]] = []
    for _, term, genes in ranked:
        redundant = False
        if go_related is not None and term[0].startswith("GO:"):
            for chosen_term, chosen_genes in selected:
                union = genes | chosen_genes
                jaccard = len(genes & chosen_genes) / len(union) if union else 0
                if (
                    chosen_term[0].startswith("GO:")
                    and go_related(term[0], chosen_term[0])
                    and jaccard >= jaccard_cutoff
                ):
                    redundant = True
                    break
        if not redundant:
            selected.append((term, genes))
        if len(selected) == limit:
            break
    return [term for term, _ in selected]


def wrap_term(name: str, term_id: str, *, width: int = 22) -> str:
    wrapped = "\n".join(
        textwrap.wrap(
            name,
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )
    )
    return f"{wrapped}\n{term_id}"


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
        or int(manifest.get("assembly_units", 0)) != 23
        or manifest.get("loss_classification") != "decayed + deleted"
    ):
        raise SystemExit("ERROR: unit enrichment manifest is not exact PASS")

    raw_rows, fields = read_tsv(args.significant_enrichment)
    go_parents, go_aliases = read_go_parents(args.go_obo)
    go_related = ancestor_checker(go_parents, go_aliases)
    pathway_names = read_pathway_names(args.kegg_pathway_names)
    required = {
        "assembly_unit_id",
        "biological_species",
        "haplotype_or_subgenome",
        "ontology",
        "term_id",
        "term_name",
        "go_namespace",
        "study_count",
        "p_fdr_bh",
        "fold_enrichment",
        "significant_fdr",
        "study_gene_ids",
    }
    if not required.issubset(fields):
        raise SystemExit("ERROR: significant table lacks required fields")
    if len(raw_rows) != int(manifest["significant_terms"]):
        raise SystemExit("ERROR: significant row count does not close")

    unit_metadata: dict[str, tuple[str, str]] = {}
    normalized: list[dict[str, object]] = []
    for row in raw_rows:
        unit = row["assembly_unit_id"]
        metadata = (
            row["biological_species"],
            row["haplotype_or_subgenome"],
        )
        if unit in unit_metadata and unit_metadata[unit] != metadata:
            raise SystemExit("ERROR: inconsistent unit metadata")
        unit_metadata[unit] = metadata
        q_value = float(row["p_fdr_bh"])
        fold = float(row["fold_enrichment"])
        study_count = int(row["study_count"])
        if (
            row["significant_fdr"] != "true"
            or not 0.0 < q_value <= 0.05
            or study_count < 1
        ):
            raise SystemExit("ERROR: non-significant row in PASS table")
        normalized.append(
            {
                **row,
                "term_name": (
                    pathway_names.get(row["term_id"], row["term_name"])
                    if row["ontology"] == "KEGG_PATHWAY"
                    else row["term_name"]
                ),
                "functional_category": functional_category(row),
                "study_count": study_count,
                "p_fdr_bh": q_value,
                "fold_enrichment": fold,
                "minus_log10_q": -math.log10(q_value),
            }
        )
    if len(unit_metadata) != 23:
        raise SystemExit("ERROR: expected exactly 23 assembly units")

    unit_order = list(unit_metadata)
    display_labels = {
        unit: format_downstream_taxon_label(
            unit_metadata[unit][0],
            (unit_metadata[unit][1],),
            abbreviate_genus=True,
            separator=" ",
        )
        for unit in unit_order
    }

    filtered = [
        row
        for row in normalized
        if int(row["study_count"]) >= 2
        and float(row["fold_enrichment"]) > 1.0
    ]
    selected: dict[str, list[tuple[str, str]]] = {}
    for category, limit, _, _ in CATEGORY_LIMITS:
        selected[category] = rank_terms(
            (
                row
                for row in filtered
                if row["functional_category"] == category
            ),
            limit=limit,
            go_related=go_related if category.startswith("GO ") else None,
        )
        if len(selected[category]) != limit:
            raise SystemExit(
                f"ERROR: fewer than {limit} eligible terms for {category}"
            )
    selected_pathways = {
        term_id for term_id, _ in selected["KEGG pathway"]
    }
    if not selected_pathways.issubset(pathway_names):
        missing = ",".join(sorted(selected_pathways - pathway_names))
        raise SystemExit(
            f"ERROR: selected KEGG pathways lack display names: {missing}"
        )

    selected_keys = {
        (category, term_id, term_name)
        for category, terms in selected.items()
        for term_id, term_name in terms
    }
    plot_rows = [
        {
            "panel": next(
                panel
                for category_name, _, panel, _ in CATEGORY_LIMITS
                if category_name == row["functional_category"]
            ),
            "functional_category": row["functional_category"],
            "term_rank": selected[row["functional_category"]].index(
                (str(row["term_id"]), str(row["term_name"]))
            )
            + 1,
            "term_id": row["term_id"],
            "term_name": row["term_name"],
            "assembly_unit_id": row["assembly_unit_id"],
            "biological_species": row["biological_species"],
            "haplotype_or_subgenome": row["haplotype_or_subgenome"],
            "display_label": display_labels[str(row["assembly_unit_id"])],
            "study_count": int(row["study_count"]),
            "study_size": int(row["study_size"]),
            "background_count": int(row["background_count"]),
            "background_size": int(row["background_size"]),
            "p_fdr_bh": float(row["p_fdr_bh"]),
            "minus_log10_q": float(row["minus_log10_q"]),
            "fold_enrichment": float(row["fold_enrichment"]),
            "study_gene_ids": row["study_gene_ids"],
        }
        for row in filtered
        if (
            str(row["functional_category"]),
            str(row["term_id"]),
            str(row["term_name"]),
        )
        in selected_keys
    ]
    seen = {
        (
            str(row["functional_category"]),
            str(row["term_id"]),
            str(row["assembly_unit_id"]),
        )
        for row in plot_rows
    }
    if len(seen) != len(plot_rows):
        raise SystemExit("ERROR: duplicate selected term–unit row")

    color_values = np.asarray(
        [float(row["minus_log10_q"]) for row in plot_rows],
        dtype=float,
    )
    color_min = max(-math.log10(0.05), float(color_values.min()))
    color_max = max(
        color_min + 0.5,
        float(np.quantile(color_values, 0.98)),
    )
    size_counts = np.asarray(
        [int(row["study_count"]) for row in plot_rows],
        dtype=float,
    )
    size_scale = 0.65

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 7.6,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.5,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    legend_counts = [
        int(round(value))
        for value in np.quantile(size_counts, [0.15, 0.50, 0.85])
    ]
    legend_counts = sorted(set(max(2, value) for value in legend_counts))
    output_root = args.output_dir
    if output_root.is_symlink():
        raise SystemExit("ERROR: refusing symlink collection output")
    existed_empty = output_root.is_dir() and not any(output_root.iterdir())
    if output_root.exists() and not existed_empty:
        raise SystemExit("ERROR: collection output must be absent or empty")
    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging_root = Path(
        tempfile.mkdtemp(
            prefix=f".{output_root.name}.",
            suffix=".tmp",
            dir=output_root.parent,
        )
    )
    bundle_records: list[dict[str, object]] = []
    try:
        for category, _, panel, suffix in CATEGORY_LIMITS:
            terms = selected[category]
            term_index = {term: i for i, term in enumerate(terms)}
            unit_index = {unit: i for i, unit in enumerate(unit_order)}
            members = [
                row
                for row in plot_rows
                if row["functional_category"] == category
            ]
            x = [
                unit_index[str(row["assembly_unit_id"])] for row in members
            ]
            y = [
                term_index[(str(row["term_id"]), str(row["term_name"]))]
                for row in members
            ]
            colors = [
                min(float(row["minus_log10_q"]), color_max)
                for row in members
            ]
            sizes = [
                size_scale * int(row["study_count"]) for row in members
            ]
            height = {
                "GO biological process": 5.5,
                "GO molecular function": 6.8,
                "GO cellular component": 4.6,
                "KEGG pathway": 5.5,
            }[category]
            figure = plt.figure(figsize=(7.2, height), dpi=180)
            axis = figure.add_axes([0.30, 0.25, 0.61, 0.67])
            axis.scatter(
                x,
                y,
                c=colors,
                s=sizes,
                cmap="viridis_r",
                vmin=color_min,
                vmax=color_max,
                edgecolors="#2f2f2f",
                linewidths=0.45,
                alpha=0.93,
            )
            axis.set_xticks(
                range(len(unit_order)),
                [display_labels[unit] for unit in unit_order],
                rotation=58,
                ha="right",
                rotation_mode="anchor",
            )
            axis.set_yticks(
                range(len(terms)),
                [
                    wrap_term(name, term_id, width=38)
                    for term_id, name in terms
                ],
            )
            axis.set_xlim(-0.65, len(unit_order) - 0.35)
            axis.set_ylim(len(terms) - 0.35, -0.65)
            axis.set_axisbelow(True)
            axis.grid(color="#e5e5e5", linewidth=0.6)
            axis.text(
                -0.10,
                1.025,
                f"({panel.lower()})",
                transform=axis.transAxes,
                fontsize=10.5,
                fontweight="bold",
                ha="left",
                va="bottom",
            )
            for spine in axis.spines.values():
                spine.set_color("#666666")
                spine.set_linewidth(0.75)

            colorbar = figure.colorbar(
                ScalarMappable(
                    norm=Normalize(vmin=color_min, vmax=color_max),
                    cmap="viridis_r",
                ),
                ax=axis,
                fraction=0.025,
                pad=0.018,
            )
            colorbar.set_label(r"$-\log_{10}$(BH q-value)")
            colorbar.ax.tick_params(labelsize=7.5)
            handles = [
                plt.scatter(
                    [],
                    [],
                    s=size_scale * value,
                    facecolor="#7a7a7a",
                    edgecolor="#2f2f2f",
                    linewidth=0.45,
                )
                for value in legend_counts
            ]
            axis.legend(
                handles,
                [f"{value} genes" for value in legend_counts],
                title="Lost genes in term",
                loc="lower right",
                bbox_to_anchor=(1.0, 1.015),
                ncol=len(legend_counts),
                frameon=False,
                borderaxespad=0,
                columnspacing=1.0,
                handletextpad=0.45,
                fontsize=7.5,
                title_fontsize=8.0,
            )

            caption = (
                f"Representative {category.lower()} enrichment among "
                "decayed-plus-deleted genes in each of 23 independently "
                "analyzed assembly units. Assembly units are shown on the "
                "horizontal axis and functional terms on the vertical axis. "
                "Point area is directly proportional to the number of lost "
                "genes assigned to the term, and color shows "
                "Benjamini–Hochberg significance. Eligible terms require BH "
                "q <= 0.05, at least two foreground genes, and fold "
                "enrichment > 1. Terms are ranked objectively by significant "
                "unit count, median significance, median fold enrichment, "
                "total contributing genes, and stable term identifier."
            )
            validation = {
                "status": "PASS_UNIT_LOSS_FUNCTIONAL_DETAIL_PANEL",
                "panel": panel,
                "functional_category": category,
                "assembly_units": len(unit_order),
                "significant_input_rows": len(raw_rows),
                "eligible_input_rows": len(filtered),
                "selected_terms": len(terms),
                "plotted_term_unit_rows": len(members),
                "selection_rule": {
                    "q_max": 0.05,
                    "study_count_min": 2,
                    "fold_enrichment_min_exclusive": 1.0,
                    "go_ancestor_descendant_jaccard_redundancy_cutoff": 0.70,
                    "ranking": [
                        "significant_unit_count_desc",
                        "median_minus_log10_q_desc",
                        "median_fold_enrichment_desc",
                        "total_study_count_desc",
                        "term_id_asc",
                    ],
                },
                "checks": {
                    "decayed_plus_deleted": True,
                    "assembly_units_not_aggregated": True,
                    "actual_gene_counts_shown_by_point_area": True,
                    "point_area_linear_in_gene_count": True,
                    "specific_terms_named": True,
                    "go_redundancy_reduced": True,
                    "kegg_pathway_names_resolved": True,
                    "latin_binomials_italic_suffixes_upright": True,
                    "axes_transposed_for_readability": True,
                    "separate_category_figures": True,
                    "figure_title_omitted": True,
                },
            }
            basename = f"{args.basename}_{suffix}"
            bundle = write_figure_bundle(
                figure=figure,
                output_dir=staging_root / suffix,
                basename=basename,
                plot_rows=members,
                plot_columns=list(members[0]),
                caption=caption,
                validation=validation,
                input_paths=[
                    args.significant_enrichment,
                    args.run_manifest,
                    args.go_obo,
                    args.kegg_pathway_names,
                ],
                dpi=600,
            )
            plt.close(figure)
            bundle_records.append(
                {
                    "panel": panel,
                    "functional_category": category,
                    "directory": suffix,
                    "basename": basename,
                    "files": [
                        {
                            "path": str(path.relative_to(staging_root)),
                            "bytes": path.stat().st_size,
                            "sha256": hashlib.sha256(
                                path.read_bytes()
                            ).hexdigest(),
                        }
                        for path in (
                            bundle.png,
                            bundle.pdf,
                            bundle.plot_data,
                            bundle.caption,
                            bundle.validation,
                            bundle.manifest,
                        )
                    ],
                }
            )

        (staging_root / "README.md").write_text(
            "# Unit-resolved functional detail\n\n"
            "Four publication-scale figures keep the 23 assembly units "
            "independent and place units on the horizontal axis and terms on "
            "the vertical axis. The `bp`, `mf`, `cc`, and `kegg` directories "
            "contain GO biological-process, GO molecular-function, GO "
            "cellular-component, and KEGG pathway panels, respectively. "
            "Point area is directly proportional to lost-gene count and "
            "color is `-log10(BH q-value)`.\n",
            encoding="utf-8",
        )
        (staging_root / "collection_manifest.json").write_text(
            json.dumps(
                {
                    "schema_version": "1.0",
                    "status": "PASS_UNIT_LOSS_FUNCTIONAL_DETAIL_COLLECTION",
                    "assembly_units": 23,
                    "loss_classification": "decayed + deleted",
                    "axes": {
                        "x": "assembly_unit",
                        "y": "functional_term",
                    },
                    "bundles": bundle_records,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        if existed_empty:
            output_root.rmdir()
        os.replace(staging_root, output_root)
    except BaseException:
        shutil.rmtree(staging_root, ignore_errors=True)
        raise

    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
