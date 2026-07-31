#!/usr/bin/env python3
"""Curate GO/KEGG enrichment for 23-unit scaffold loss events."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping


CATEGORIES = (
    "GO biological process",
    "GO molecular function",
    "GO cellular component",
    "KEGG orthology",
    "KEGG pathway",
)


class SummaryError(ValueError):
    """Raised when scaffold enrichment outputs do not close."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--foreground-summary", required=True, type=Path)
    parser.add_argument("--significant-enrichment", required=True, type=Path)
    parser.add_argument("--enrichment-manifest", required=True, type=Path)
    parser.add_argument("--foreground-metadata", required=True, type=Path)
    parser.add_argument("--foreground-manifest", required=True, type=Path)
    parser.add_argument("--expected-foregrounds", type=int, default=39)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
        fields = reader.fieldnames or []
    if not rows or len(fields) != len(set(fields)):
        raise SummaryError(f"{path.name}: invalid or empty table")
    return rows, fields


def category(row: Mapping[str, str]) -> str:
    if row["ontology"] == "GO":
        return {
            "biological_process": "GO biological process",
            "molecular_function": "GO molecular function",
            "cellular_component": "GO cellular component",
        }[row["go_namespace"]]
    if row["ontology"] == "KEGG_KO":
        return "KEGG orthology"
    if row["ontology"] == "KEGG_PATHWAY":
        return "KEGG pathway"
    raise SummaryError(f"unsupported ontology {row['ontology']!r}")


def write_tsv(
    path: Path,
    fields: list[str],
    rows: Iterable[Mapping[str, object]],
    *,
    gz: bool = False,
) -> None:
    opener = gzip.open if gz else open
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
        raise SummaryError(f"output directory exists: {args.output_dir}")
    inputs = {
        "foreground_summary": args.foreground_summary,
        "significant_enrichment": args.significant_enrichment,
        "enrichment_manifest": args.enrichment_manifest,
        "foreground_metadata": args.foreground_metadata,
        "foreground_manifest": args.foreground_manifest,
    }
    enrichment_manifest = json.loads(
        args.enrichment_manifest.read_text(encoding="utf-8")
    )
    foreground_manifest = json.loads(
        args.foreground_manifest.read_text(encoding="utf-8")
    )
    if (
        enrichment_manifest.get("status") != "PASS_UNIT_SCAFFOLD_GO_KEGG"
        or enrichment_manifest.get("analysis_profile") != "scaffold"
        or foreground_manifest.get("status")
        != "PASS_UNIT_SCAFFOLD_FUNCTION_FOREGROUNDS"
    ):
        raise SummaryError("input manifests are not scaffold PASS")

    metadata_rows, metadata_fields = read_tsv(args.foreground_metadata)
    needed_metadata = {
        "foreground_id",
        "analysis_scope",
        "background_scope",
        "foreground_gene_count",
        "branch_id",
        "descendant_lineage_count",
        "descendant_lineages",
        "node_type",
        "node_name",
        "parent_node_id",
        "minimum_leaf_plot_order",
    }
    if not needed_metadata.issubset(metadata_fields):
        raise SummaryError("foreground metadata lacks required columns")
    if len(metadata_rows) != args.expected_foregrounds:
        raise SummaryError("unexpected scaffold foreground count")
    metadata = {row["foreground_id"]: row for row in metadata_rows}
    if len(metadata) != len(metadata_rows):
        raise SummaryError("duplicate scaffold foreground")

    summary_rows, summary_fields = read_tsv(args.foreground_summary)
    needed_summary = {
        "foreground_id",
        "ontology",
        "requested_gene_count",
        "annotation_coverage",
    }
    if not needed_summary.issubset(summary_fields):
        raise SummaryError("foreground summary lacks required columns")
    if len(summary_rows) != args.expected_foregrounds * 3:
        raise SummaryError("foreground summary is not an exact node x 3 grid")
    coverage: dict[tuple[str, str], float] = {}
    for row in summary_rows:
        foreground_id = row["foreground_id"]
        if (
            foreground_id not in metadata
            or int(row["requested_gene_count"])
            != int(metadata[foreground_id]["foreground_gene_count"])
        ):
            raise SummaryError("summary and foreground metadata disagree")
        coverage[(foreground_id, row["ontology"])] = float(
            row["annotation_coverage"]
        )

    significant, significant_fields = read_tsv(
        args.significant_enrichment
    )
    required_significant = {
        "foreground_id",
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
    if not required_significant.issubset(significant_fields):
        raise SummaryError("significant table lacks required columns")
    normalized: list[dict[str, object]] = []
    counts: Counter[tuple[str, str]] = Counter()
    for row in significant:
        foreground_id = row["foreground_id"]
        if (
            foreground_id not in metadata
            or row["significant_fdr"] != "true"
            or float(row["p_fdr_bh"]) > 0.05
            or int(row["study_count"]) < 2
            or float(row["fold_enrichment"]) <= 1
        ):
            raise SummaryError("invalid significant enrichment row")
        functional_category = category(row)
        counts[(foreground_id, functional_category)] += 1
        meta = metadata[foreground_id]
        normalized.append(
            {
                **{field: row[field] for field in significant_fields},
                "node_type": meta["node_type"],
                "node_name": meta["node_name"],
                "parent_node_id": meta["parent_node_id"],
                "minimum_leaf_plot_order": meta[
                    "minimum_leaf_plot_order"
                ],
                "functional_category": functional_category,
            }
        )

    node_summary: list[dict[str, object]] = []
    for meta in metadata_rows:
        foreground_id = meta["foreground_id"]
        node_summary.append(
            {
                "foreground_id": foreground_id,
                "branch_id": meta["branch_id"],
                "node_type": meta["node_type"],
                "node_name": meta["node_name"],
                "parent_node_id": meta["parent_node_id"],
                "descendant_unit_count": int(
                    meta["descendant_lineage_count"]
                ),
                "descendant_units": meta["descendant_lineages"],
                "minimum_leaf_plot_order": int(
                    meta["minimum_leaf_plot_order"]
                ),
                "loss_event_gene_count": int(
                    meta["foreground_gene_count"]
                ),
                "go_annotation_coverage": coverage[
                    (foreground_id, "GO")
                ],
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
                    counts[(foreground_id, item)] for item in CATEGORIES
                ),
            }
        )

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.",
            dir=args.output_dir.parent,
        )
    )
    try:
        node_path = temporary / "node_summary.tsv"
        significant_path = temporary / "enrichment_significant.tsv.gz"
        write_tsv(node_path, list(node_summary[0]), node_summary)
        write_tsv(
            significant_path,
            list(normalized[0]),
            normalized,
            gz=True,
        )
        outputs = [node_path, significant_path]
        manifest = {
            "schema_version": 1,
            "status": "PASS_UNIT_SCAFFOLD_GO_KEGG_SUMMARY",
            "foregrounds": len(metadata_rows),
            "foreground_memberships": enrichment_manifest["counts"][
                "foreground_memberships"
            ],
            "resolved_background_genes": enrichment_manifest["counts"][
                "resolved_scaffold_background_genes"
            ],
            "significant_terms": len(normalized),
            "loss_classification": "decayed + deleted",
            "foreground_definition": (
                "maximal loss-event gene set at one node or terminal of the "
                "topology-only 23-unit scaffold"
            ),
            "background_definition": (
                "reference genes resolved across all 23 assembly units"
            ),
            "inputs": [
                {
                    "role": role,
                    "basename": path.name,
                    "sha256": sha256(path),
                }
                for role, path in inputs.items()
            ],
            "outputs": [
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in outputs
            ],
        }
        (temporary / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, args.output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (SummaryError, OSError, json.JSONDecodeError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
