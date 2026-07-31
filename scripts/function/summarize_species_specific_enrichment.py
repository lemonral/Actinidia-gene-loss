#!/usr/bin/env python3
"""Curate single-terminal article-method GO/KEGG enrichment results."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


csv.field_size_limit(sys.maxsize)

CATEGORY_ORDER = (
    "GO biological process",
    "GO molecular function",
    "GO cellular component",
    "KEGG orthology",
    "KEGG pathway",
)
GO_CATEGORY = {
    "biological_process": "GO biological process",
    "molecular_function": "GO molecular function",
    "cellular_component": "GO cellular component",
}


class SummaryError(ValueError):
    """Raised when enrichment output is inconsistent with pure terminal sets."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--foreground-summary", required=True, type=Path)
    parser.add_argument("--significant-enrichment", required=True, type=Path)
    parser.add_argument("--enrichment-manifest", required=True, type=Path)
    parser.add_argument("--foreground-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-lineages", type=int, default=13)
    parser.add_argument("--expected-genes", type=int, default=1167)
    parser.add_argument("--top-terms-per-category", type=int, default=3)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = [dict(row) for row in reader]
    if not fields or len(fields) != len(set(fields)):
        raise SummaryError(f"{path.name}: invalid TSV header")
    return rows, fields


def category(row: dict[str, str]) -> str:
    if row["ontology"] == "GO":
        if row["go_namespace"] not in GO_CATEGORY:
            raise SummaryError(
                f"unsupported GO namespace {row['go_namespace']!r}"
            )
        return GO_CATEGORY[row["go_namespace"]]
    if row["ontology"] == "KEGG_KO":
        return "KEGG orthology"
    if row["ontology"] == "KEGG_PATHWAY":
        return "KEGG pathway"
    raise SummaryError(f"unsupported ontology {row['ontology']!r}")


def write_tsv(
    path: Path,
    fields: list[str],
    rows: list[dict[str, object]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise SummaryError(f"output directory already exists: {args.output_dir}")
    if args.top_terms_per_category < 1:
        raise SummaryError("top-terms-per-category must be positive")
    enrichment_manifest = json.loads(
        args.enrichment_manifest.read_text(encoding="utf-8")
    )
    foreground_manifest = json.loads(
        args.foreground_manifest.read_text(encoding="utf-8")
    )
    if enrichment_manifest.get("status") != "PASS_TREE_AWARE_MANUSCRIPT_GO_KEGG":
        raise SummaryError("enrichment manifest is not PASS")
    if (
        foreground_manifest.get("status")
        != "PASS_ARTICLE_METHOD_SINGLE_TERMINAL_FOREGROUNDS"
    ):
        raise SummaryError("foreground manifest is not pure single-terminal PASS")

    summary_rows, summary_fields = read_tsv(args.foreground_summary)
    needed_summary = {
        "foreground_id",
        "analysis_scope",
        "branch_id",
        "descendant_lineage_count",
        "descendant_lineages",
        "ontology",
        "requested_gene_count",
        "annotated_study_gene_count",
        "annotation_coverage",
    }
    if not needed_summary.issubset(summary_fields):
        raise SummaryError("foreground summary is missing required columns")
    if len(summary_rows) != args.expected_lineages * 3:
        raise SummaryError("foreground summary is not the exact 13 x 3 grid")

    species_foreground: dict[str, str] = {}
    gene_count: dict[str, int] = {}
    coverage: dict[tuple[str, str], float] = {}
    for row in summary_rows:
        if (
            row["analysis_scope"]
            != "single_terminal_branch_species_specific_loss"
            or not row["branch_id"].startswith("terminal__")
            or row["descendant_lineage_count"] != "1"
        ):
            raise SummaryError("non-species-specific foreground present")
        species = row["descendant_lineages"]
        foreground_id = row["foreground_id"]
        species_foreground.setdefault(species, foreground_id)
        if species_foreground[species] != foreground_id:
            raise SummaryError("species maps to multiple foregrounds")
        current_count = int(row["requested_gene_count"])
        gene_count.setdefault(species, current_count)
        if gene_count[species] != current_count:
            raise SummaryError("foreground gene count differs by ontology")
        coverage[(species, row["ontology"])] = float(row["annotation_coverage"])
    if len(species_foreground) != args.expected_lineages:
        raise SummaryError("terminal lineage count does not close")
    if sum(gene_count.values()) != args.expected_genes:
        raise SummaryError("species-specific foreground genes do not close")

    significant, significant_fields = read_tsv(args.significant_enrichment)
    required_significant = {
        "foreground_id",
        "analysis_scope",
        "branch_id",
        "descendant_lineages",
        "ontology",
        "term_id",
        "term_name",
        "go_namespace",
        "study_count",
        "study_size",
        "background_count",
        "background_size",
        "p_fdr_bh",
        "fold_enrichment",
        "significant_fdr",
    }
    if not required_significant.issubset(significant_fields):
        raise SummaryError("significant enrichment table is missing required columns")
    allowed_foregrounds = set(species_foreground.values())
    counts: Counter[tuple[str, str]] = Counter()
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    normalized_significant: list[dict[str, object]] = []
    for row in significant:
        if (
            row["foreground_id"] not in allowed_foregrounds
            or row["analysis_scope"]
            != "single_terminal_branch_species_specific_loss"
            or row["significant_fdr"] != "true"
        ):
            raise SummaryError("significant table contains an invalid foreground row")
        species = row["descendant_lineages"]
        functional_category = category(row)
        counts[(species, functional_category)] += 1
        grouped[(species, functional_category)].append(row)
        normalized_significant.append(
            {
                **{field: row[field] for field in significant_fields},
                "functional_category": functional_category,
            }
        )

    public_summary: list[dict[str, object]] = []
    for species in sorted(species_foreground):
        public_summary.append(
            {
                "biological_species": species,
                "foreground_id": species_foreground[species],
                "species_specific_loss_gene_count": gene_count[species],
                "go_annotation_coverage": coverage[(species, "GO")],
                "kegg_ko_annotation_coverage": coverage[(species, "KEGG_KO")],
                "kegg_pathway_annotation_coverage": coverage[
                    (species, "KEGG_PATHWAY")
                ],
                "go_biological_process_significant_terms": counts[
                    (species, "GO biological process")
                ],
                "go_molecular_function_significant_terms": counts[
                    (species, "GO molecular function")
                ],
                "go_cellular_component_significant_terms": counts[
                    (species, "GO cellular component")
                ],
                "kegg_orthology_significant_terms": counts[
                    (species, "KEGG orthology")
                ],
                "kegg_pathway_significant_terms": counts[
                    (species, "KEGG pathway")
                ],
                "total_significant_terms": sum(
                    counts[(species, item)] for item in CATEGORY_ORDER
                ),
            }
        )

    top_terms: list[dict[str, object]] = []
    for species in sorted(species_foreground):
        for functional_category in CATEGORY_ORDER:
            rows = sorted(
                grouped.get((species, functional_category), []),
                key=lambda row: (
                    float(row["p_fdr_bh"]),
                    -float(row["fold_enrichment"]),
                    row["term_id"],
                ),
            )
            for rank, row in enumerate(
                rows[: args.top_terms_per_category],
                start=1,
            ):
                top_terms.append(
                    {
                        "biological_species": species,
                        "functional_category": functional_category,
                        "rank_within_species_category": rank,
                        "term_id": row["term_id"],
                        "term_name": row["term_name"],
                        "study_count": int(row["study_count"]),
                        "study_size": int(row["study_size"]),
                        "background_count": int(row["background_count"]),
                        "background_size": int(row["background_size"]),
                        "fold_enrichment": float(row["fold_enrichment"]),
                        "p_fdr_bh": float(row["p_fdr_bh"]),
                    }
                )

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.tmp.",
            dir=args.output_dir.parent,
        )
    )
    try:
        summary_path = staging / "species_summary.tsv"
        top_path = staging / "top_terms_by_species_category.tsv"
        significant_path = staging / "enrichment_significant.tsv"
        write_tsv(summary_path, list(public_summary[0]), public_summary)
        write_tsv(top_path, list(top_terms[0]), top_terms)
        write_tsv(
            significant_path,
            [*significant_fields, "functional_category"],
            normalized_significant,
        )
        output_paths = [summary_path, top_path, significant_path]
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_ARTICLE_METHOD_SPECIES_SPECIFIC_GO_KEGG",
            "loss_classification": "article-method decayed + deleted",
            "foreground_definition": (
                "exactly one terminal-branch loss; recurrent, internal, partial, "
                "and unknown patterns excluded"
            ),
            "test_definition": (
                "one-sided hypergeometric over-representation; BH within "
                "species foreground and ontology; q <= 0.05"
            ),
            "lineages": len(public_summary),
            "species_specific_loss_genes": sum(gene_count.values()),
            "significant_terms": len(significant),
            "top_terms_per_species_category": args.top_terms_per_category,
            "inputs": [
                {
                    "role": role,
                    "basename": path.name,
                    "sha256": sha256(path),
                }
                for role, path in (
                    ("foreground_summary", args.foreground_summary),
                    ("significant_enrichment", args.significant_enrichment),
                    ("enrichment_manifest", args.enrichment_manifest),
                    ("foreground_manifest", args.foreground_manifest),
                )
            ],
            "outputs": [
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in output_paths
            ],
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, args.output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (SummaryError, OSError, csv.Error, UnicodeError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
