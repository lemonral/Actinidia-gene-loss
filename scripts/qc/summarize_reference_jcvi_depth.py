#!/usr/bin/env python3
"""Summarize bidirectional JCVI gene-depth coverage for fixed reference pairs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

if __package__:
    from .jcvi_depth import load_anchor_blocks, load_bed, summarize_depth
else:
    from jcvi_depth import load_anchor_blocks, load_bed, summarize_depth


SCRIPT_VERSION = "1.0.0"
STATUS = "PASS_CLEMATOCLETHRA_ACTINIDIA_JCVI_GENE_DEPTH_SUMMARY"
REQUIRED_COLUMNS = (
    "plot_order",
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "source_class",
    "anchors",
    "target_bed",
)
OUTPUT_COLUMNS = (
    "plot_order",
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "source_class",
    "anchor_blocks",
    "anchor_pair_rows",
    "reference_bed_genes",
    "reference_nonzero_depth_genes",
    "reference_gene_depth_coverage_percent",
    "target_bed_genes",
    "target_nonzero_depth_genes",
    "target_gene_depth_coverage_percent",
    "metric_definition",
)
CHROMOSOME_OUTPUT_COLUMNS = (
    "plot_order",
    "assembly_unit_id",
    "reference_chromosome",
    "reference_chromosome_genes",
    "reference_nonzero_depth_genes",
    "reference_gene_depth_coverage_percent",
)
METRIC_DEFINITION = (
    "100 * union([min_gene_index,max_gene_index) across raw JCVI anchor blocks) "
    "/ total BED genes"
)


class SummaryError(RuntimeError):
    """Raised when fixed JCVI sources cannot be summarized exactly."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise SummaryError(f"{label} is missing or empty: {resolved}")
    return resolved


def read_sources(path: Path) -> list[dict[str, str]]:
    path = require_file(path, "JCVI source table")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SummaryError("JCVI source table has no header")
        missing = sorted(set(REQUIRED_COLUMNS).difference(reader.fieldnames))
        if missing:
            raise SummaryError(f"JCVI source table is missing columns: {missing}")
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
            if any((value or "").strip() for value in row.values())
        ]
    if len(rows) != 23:
        raise SummaryError(f"expected 23 JCVI source rows, found {len(rows)}")
    orders = [int(row["plot_order"]) for row in rows]
    if sorted(orders) != list(range(1, 24)):
        raise SummaryError("plot_order must contain each integer from 1 through 23")
    units = [row["assembly_unit_id"] for row in rows]
    if any(not value for value in units) or len(units) != len(set(units)):
        raise SummaryError("assembly_unit_id values must be nonempty and unique")
    return sorted(rows, key=lambda row: int(row["plot_order"]))


def resolve_source(root: Path, value: str, label: str) -> Path:
    if not value:
        raise SummaryError(f"{label} source path is empty")
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = root / candidate
    return require_file(candidate, label)


def merge_half_open(intervals: list[tuple[int, int]]) -> int:
    ordered = sorted((start, end) for start, end in intervals if end > start)
    if not ordered:
        return 0
    covered = 0
    start, end = ordered[0]
    for next_start, next_end in ordered[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            covered += end - start
            start, end = next_start, next_end
    return covered + end - start


def summarize_reference_chromosomes(
    blocks: list[list[tuple[str, str]]],
    reference_bed: object,
) -> list[dict[str, int | float | str]]:
    sequence_order: list[str] = []
    totals: dict[str, int] = {}
    for sequence_id, _start, _end, _accession in reference_bed.rows:
        if sequence_id not in totals:
            sequence_order.append(sequence_id)
            totals[sequence_id] = 0
        totals[sequence_id] += 1
    intervals: dict[str, list[tuple[int, int]]] = {
        sequence_id: [] for sequence_id in sequence_order
    }
    for block in blocks:
        located = [reference_bed.order[pair[0]] for pair in block]
        sequence_ids = {sequence_id for _index, sequence_id in located}
        if len(sequence_ids) != 1:
            raise SummaryError("reference anchor block crosses BED sequence IDs")
        sequence_id = next(iter(sequence_ids))
        indices = [index for index, _sequence_id in located]
        intervals[sequence_id].append((min(indices), max(indices)))
    result: list[dict[str, int | float | str]] = []
    for sequence_id in sequence_order:
        covered = merge_half_open(intervals[sequence_id])
        total = totals[sequence_id]
        result.append(
            {
                "reference_chromosome": sequence_id,
                "reference_chromosome_genes": total,
                "reference_nonzero_depth_genes": covered,
                "reference_gene_depth_coverage_percent": covered * 100.0 / total,
            }
        )
    return result


def summarize(
    sources_path: Path,
    source_root: Path,
    reference_bed_path: Path,
    output_dir: Path,
) -> None:
    sources_path = require_file(sources_path, "JCVI source table")
    source_root = source_root.expanduser().resolve()
    if not source_root.is_dir():
        raise SummaryError(f"source root is not a directory: {source_root}")
    reference_bed_path = require_file(reference_bed_path, "reference BED")
    sources = read_sources(sources_path)
    reference_bed = load_bed(reference_bed_path)

    summary_rows: list[dict[str, object]] = []
    chromosome_rows: list[dict[str, object]] = []
    exact_inputs: list[Path] = [sources_path, reference_bed_path]
    for source in sources:
        unit = source["assembly_unit_id"]
        anchors_path = resolve_source(
            source_root, source["anchors"], f"{unit} raw JCVI anchors"
        )
        target_bed_path = resolve_source(
            source_root, source["target_bed"], f"{unit} target BED"
        )
        blocks, pair_rows = load_anchor_blocks(anchors_path)
        target_bed = load_bed(target_bed_path)
        reference = summarize_depth(blocks, reference_bed, 0)
        target = summarize_depth(blocks, target_bed, 1)
        summary_rows.append(
            {
                "plot_order": int(source["plot_order"]),
                "assembly_unit_id": unit,
                "biological_species": source["biological_species"],
                "haplotype_or_subgenome": source["haplotype_or_subgenome"],
                "source_class": source["source_class"],
                "anchor_blocks": len(blocks),
                "anchor_pair_rows": pair_rows,
                "reference_bed_genes": int(reference["total"]),
                "reference_nonzero_depth_genes": int(reference["nonzero"]),
                "reference_gene_depth_coverage_percent": (
                    f"{float(reference['coverage']):.9f}"
                ),
                "target_bed_genes": int(target["total"]),
                "target_nonzero_depth_genes": int(target["nonzero"]),
                "target_gene_depth_coverage_percent": (
                    f"{float(target['coverage']):.9f}"
                ),
                "metric_definition": METRIC_DEFINITION,
            }
        )
        chromosome_summary = summarize_reference_chromosomes(blocks, reference_bed)
        if sum(
            int(row["reference_nonzero_depth_genes"])
            for row in chromosome_summary
        ) != int(reference["nonzero"]):
            raise SummaryError(
                f"{unit}: chromosome gene-depth counts do not sum to overall coverage"
            )
        for chromosome_row in chromosome_summary:
            chromosome_rows.append(
                {
                    "plot_order": int(source["plot_order"]),
                    "assembly_unit_id": unit,
                    "reference_chromosome": chromosome_row[
                        "reference_chromosome"
                    ],
                    "reference_chromosome_genes": chromosome_row[
                        "reference_chromosome_genes"
                    ],
                    "reference_nonzero_depth_genes": chromosome_row[
                        "reference_nonzero_depth_genes"
                    ],
                    "reference_gene_depth_coverage_percent": (
                        f"{float(chromosome_row['reference_gene_depth_coverage_percent']):.9f}"
                    ),
                }
            )
        exact_inputs.extend((anchors_path, target_bed_path))

    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists():
        raise SummaryError(f"refusing existing output directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.tmp.", dir=output_dir.parent)
    )
    try:
        summary_path = temporary / "clematoclethra_actinidia_jcvi_gene_depth.tsv"
        with summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(summary_rows)
        chromosome_path = (
            temporary
            / "clematoclethra_actinidia_jcvi_gene_depth.chromosome_summary.tsv"
        )
        with chromosome_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=CHROMOSOME_OUTPUT_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(chromosome_rows)
        reference_values = [
            float(row["reference_gene_depth_coverage_percent"]) for row in summary_rows
        ]
        target_values = [
            float(row["target_gene_depth_coverage_percent"]) for row in summary_rows
        ]
        validation = {
            "schema_version": 1,
            "status": STATUS,
            "script": "scripts/qc/summarize_reference_jcvi_depth.py",
            "script_version": SCRIPT_VERSION,
            "unit_count": len(summary_rows),
            "reference_chromosome_count": len(
                {row["reference_chromosome"] for row in chromosome_rows}
            ),
            "chromosome_grid_rows": len(chromosome_rows),
            "metric_definition": METRIC_DEFINITION,
            "reference_bed_gene_count": len(reference_bed.rows),
            "minimum_reference_gene_depth_coverage_percent": min(reference_values),
            "maximum_reference_gene_depth_coverage_percent": max(reference_values),
            "minimum_target_gene_depth_coverage_percent": min(target_values),
            "maximum_target_gene_depth_coverage_percent": max(target_values),
            "inputs": [
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in exact_inputs
            ],
            "output": {
                "basename": summary_path.name,
                "bytes": summary_path.stat().st_size,
                "sha256": sha256(summary_path),
            },
            "chromosome_output": {
                "basename": chromosome_path.name,
                "bytes": chromosome_path.stat().st_size,
                "sha256": sha256(chromosome_path),
            },
        }
        (temporary / "validation.json").write_text(
            json.dumps(validation, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--reference-bed", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        summarize(
            args.sources,
            args.source_root,
            args.reference_bed,
            args.output_dir,
        )
        print(f"PASS: {args.output_dir}")
        return 0
    except (OSError, UnicodeError, ValueError, SummaryError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
