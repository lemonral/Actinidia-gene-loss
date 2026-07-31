#!/usr/bin/env python3
"""Curate GO/KEGG enrichment for 23 independent assembly-unit loss sets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping


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
    """Raised when unit-level enrichment outputs are inconsistent."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--foreground-summary", required=True, type=Path)
    parser.add_argument("--significant-enrichment", required=True, type=Path)
    parser.add_argument("--enrichment-manifest", required=True, type=Path)
    parser.add_argument("--foreground-manifest", required=True, type=Path)
    parser.add_argument("--foreground-metadata", required=True, type=Path)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--top-terms-per-category", type=int, default=3)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = [dict(row) for row in reader]
    if not fields or len(fields) != len(set(fields)):
        raise SummaryError(f"{path.name}: invalid TSV header")
    return rows, fields


def category(row: Mapping[str, str]) -> str:
    ontology = row["ontology"]
    if ontology == "GO":
        namespace = row["go_namespace"]
        if namespace not in GO_CATEGORY:
            raise SummaryError(f"unsupported GO namespace {namespace!r}")
        return GO_CATEGORY[namespace]
    if ontology == "KEGG_KO":
        return "KEGG orthology"
    if ontology == "KEGG_PATHWAY":
        return "KEGG pathway"
    raise SummaryError(f"unsupported ontology {ontology!r}")


def write_tsv(
    path: Path,
    fields: list[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
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
    if (
        enrichment_manifest.get("status")
        != "PASS_UNIT_ARTICLE_METHOD_GO_KEGG"
        or enrichment_manifest.get("analysis_profile") != "assembly_unit"
    ):
        raise SummaryError("enrichment manifest is not assembly-unit PASS")
    if (
        foreground_manifest.get("status")
        != "PASS_UNIT_ARTICLE_METHOD_FOREGROUNDS"
    ):
        raise SummaryError("foreground manifest is not assembly-unit PASS")

    metadata_rows, metadata_fields = read_tsv(args.foreground_metadata)
    needed_metadata = {
        "foreground_id",
        "analysis_scope",
        "assembly_unit_id",
        "biological_species",
        "haplotype_or_subgenome",
        "foreground_gene_count",
    }
    if not needed_metadata.issubset(metadata_fields):
        raise SummaryError("foreground metadata is missing required columns")
    if len(metadata_rows) != args.expected_units:
        raise SummaryError("foreground metadata does not contain 23 units")
    metadata = {row["foreground_id"]: row for row in metadata_rows}
    if len(metadata) != len(metadata_rows):
        raise SummaryError("duplicate foreground metadata")
    if any(
        row["analysis_scope"] != "assembly_unit_article_method_loss"
        for row in metadata_rows
    ):
        raise SummaryError("non-unit foreground present")

    summary_rows, summary_fields = read_tsv(args.foreground_summary)
    needed_summary = {
        "foreground_id",
        "analysis_scope",
        "ontology",
        "requested_gene_count",
        "annotated_study_gene_count",
        "annotation_coverage",
    }
    if not needed_summary.issubset(summary_fields):
        raise SummaryError("enrichment summary is missing required columns")
    if len(summary_rows) != args.expected_units * 3:
        raise SummaryError("enrichment summary is not an exact 23 x 3 grid")
    coverage: dict[tuple[str, str], float] = {}
    for row in summary_rows:
        foreground_id = row["foreground_id"]
        if (
            foreground_id not in metadata
            or row["analysis_scope"] != "assembly_unit_article_method_loss"
            or int(row["requested_gene_count"])
            != int(metadata[foreground_id]["foreground_gene_count"])
        ):
            raise SummaryError("foreground summary and metadata disagree")
        coverage[(foreground_id, row["ontology"])] = float(
            row["annotation_coverage"]
        )

    significant, significant_fields = read_tsv(args.significant_enrichment)
    required_significant = {
        "foreground_id",
        "analysis_scope",
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
        raise SummaryError("significant enrichment table is missing columns")
    counts: Counter[tuple[str, str]] = Counter()
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    normalized_significant: list[dict[str, object]] = []
    for row in significant:
        foreground_id = row["foreground_id"]
        if (
            foreground_id not in metadata
            or row["analysis_scope"] != "assembly_unit_article_method_loss"
            or row["significant_fdr"] != "true"
        ):
            raise SummaryError("significant table contains an invalid row")
        functional_category = category(row)
        counts[(foreground_id, functional_category)] += 1
        grouped[(foreground_id, functional_category)].append(row)
        normalized_significant.append(
            {
                **{field: row[field] for field in significant_fields},
                "assembly_unit_id": metadata[foreground_id][
                    "assembly_unit_id"
                ],
                "biological_species": metadata[foreground_id][
                    "biological_species"
                ],
                "haplotype_or_subgenome": metadata[foreground_id][
                    "haplotype_or_subgenome"
                ],
                "functional_category": functional_category,
            }
        )

    public_summary: list[dict[str, object]] = []
    top_terms: list[dict[str, object]] = []
    for meta in metadata_rows:
        foreground_id = meta["foreground_id"]
        unit = meta["assembly_unit_id"]
        public_summary.append(
            {
                "assembly_unit_id": unit,
                "biological_species": meta["biological_species"],
                "haplotype_or_subgenome": meta[
                    "haplotype_or_subgenome"
                ],
                "foreground_id": foreground_id,
                "loss_gene_count": int(meta["foreground_gene_count"]),
                "go_annotation_coverage": coverage[(foreground_id, "GO")],
                "kegg_ko_annotation_coverage": coverage[
                    (foreground_id, "KEGG_KO")
                ],
                "kegg_pathway_annotation_coverage": coverage[
                    (foreground_id, "KEGG_PATHWAY")
                ],
                "go_biological_process_significant_terms": counts[
                    (foreground_id, "GO biological process")
                ],
                "go_molecular_function_significant_terms": counts[
                    (foreground_id, "GO molecular function")
                ],
                "go_cellular_component_significant_terms": counts[
                    (foreground_id, "GO cellular component")
                ],
                "kegg_orthology_significant_terms": counts[
                    (foreground_id, "KEGG orthology")
                ],
                "kegg_pathway_significant_terms": counts[
                    (foreground_id, "KEGG pathway")
                ],
                "total_significant_terms": sum(
                    counts[(foreground_id, item)]
                    for item in CATEGORY_ORDER
                ),
            }
        )
        for functional_category in CATEGORY_ORDER:
            ranked = sorted(
                grouped.get((foreground_id, functional_category), []),
                key=lambda row: (
                    float(row["p_fdr_bh"]),
                    -float(row["fold_enrichment"]),
                    row["term_id"],
                ),
            )
            for rank, row in enumerate(
                ranked[: args.top_terms_per_category],
                start=1,
            ):
                top_terms.append(
                    {
                        "assembly_unit_id": unit,
                        "biological_species": meta["biological_species"],
                        "haplotype_or_subgenome": meta[
                            "haplotype_or_subgenome"
                        ],
                        "functional_category": functional_category,
                        "rank_within_unit_category": rank,
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
        summary_path = staging / "unit_summary.tsv"
        top_path = staging / "top_terms_by_unit_category.tsv"
        significant_path = staging / "enrichment_significant.tsv.gz"
        write_tsv(summary_path, list(public_summary[0]), public_summary)
        write_tsv(
            top_path,
            list(top_terms[0]) if top_terms else [
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "functional_category",
                "rank_within_unit_category",
                "term_id",
                "term_name",
                "study_count",
                "study_size",
                "background_count",
                "background_size",
                "fold_enrichment",
                "p_fdr_bh",
            ],
            top_terms,
        )
        write_tsv(
            significant_path,
            [
                *significant_fields,
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "functional_category",
            ],
            normalized_significant,
        )
        output_paths = [summary_path, top_path, significant_path]
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_UNIT_ARTICLE_METHOD_GO_KEGG_SUMMARY",
            "loss_classification": "decayed + deleted",
            "foreground_definition": (
                "one independent foreground per assembly unit; no species, "
                "haplotype, or subgenome aggregation"
            ),
            "background_definition": (
                "matching unit retained + decayed + deleted; "
                "not_called_loss excluded"
            ),
            "test_definition": (
                "one-sided hypergeometric over-representation; BH within "
                "assembly-unit foreground and ontology; q <= 0.05"
            ),
            "assembly_units": len(public_summary),
            "foreground_memberships": sum(
                int(row["loss_gene_count"]) for row in public_summary
            ),
            "significant_terms": len(significant),
            "top_terms_per_unit_category": args.top_terms_per_category,
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
                    ("foreground_metadata", args.foreground_metadata),
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
    except (SummaryError, OSError, csv.Error, UnicodeError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
