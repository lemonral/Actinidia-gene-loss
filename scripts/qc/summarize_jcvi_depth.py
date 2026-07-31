#!/usr/bin/env python3
"""Validate and summarize one raw-anchor JCVI non-zero gene-depth comparison."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

if __package__:
    from .jcvi_depth import (
        ReconstructionError,
        load_anchor_blocks,
        load_bed,
        summarize_depth,
    )
    from .prepare_jcvi_bed import read_fasta_ids, sha256
else:
    from jcvi_depth import (
        ReconstructionError,
        load_anchor_blocks,
        load_bed,
        summarize_depth,
    )
    from prepare_jcvi_bed import read_fasta_ids, sha256


SCRIPT_VERSION = "2.0.0"


class DepthSummaryError(RuntimeError):
    """Raised when a JCVI result or its denominators fail validation."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--display-name", required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--reference-protein", required=True, type=Path)
    parser.add_argument("--reference-bed", required=True, type=Path)
    parser.add_argument("--query-protein", required=True, type=Path)
    parser.add_argument("--query-bed", required=True, type=Path)
    parser.add_argument("--anchors", required=True, type=Path)
    parser.add_argument("--depthfile", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--minimum-block-size", type=int, default=4)
    parser.add_argument(
        "--allowed-reference-bed-only-ids",
        type=Path,
        help=(
            "Optional one-ID-per-line allow-list. By default the reference BED and "
            "protein FASTA identifier sets must be identical. When supplied, the "
            "BED-only set must exactly equal this file; FASTA-only IDs are never allowed."
        ),
    )
    return parser.parse_args()


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise DepthSummaryError(f"{label} is missing or empty: {resolved}")
    return resolved


def load_identifier_allowlist(path: Path) -> set[str]:
    """Load a strict one-identifier-per-line allow-list."""

    identifiers: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            value = raw.strip()
            if not value or value.startswith("#"):
                continue
            if any(character.isspace() for character in value):
                raise DepthSummaryError(
                    f"{path}:{line_number}: expected exactly one identifier"
                )
            if value in identifiers:
                raise DepthSummaryError(
                    f"{path}:{line_number}: duplicate allowed identifier {value!r}"
                )
            identifiers.add(value)
    if not identifiers:
        raise DepthSummaryError(f"Reference BED-only allow-list contains no identifiers: {path}")
    return identifiers


def validate_depthfile(
    path: Path,
    reference_ids: set[str],
    query_ids: set[str],
    expected_reference: dict[str, int | float],
    expected_query: dict[str, int | float],
) -> dict[str, object]:
    overlap = reference_ids.intersection(query_ids)
    if overlap:
        raise DepthSummaryError(
            f"Reference and query BED identifier sets overlap; examples: {sorted(overlap)[:5]}"
        )
    observed: Counter[str] = Counter()
    depths: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise DepthSummaryError(f"{path}:{line_number}: expected ID and integer depth")
            identifier, depth_text = fields
            try:
                depth = int(depth_text)
            except ValueError as exc:
                raise DepthSummaryError(f"{path}:{line_number}: non-integer depth") from exc
            if depth < 0:
                raise DepthSummaryError(f"{path}:{line_number}: negative depth")
            observed[identifier] += 1
            depths[identifier] = depth
    expected_ids = reference_ids.union(query_ids)
    if set(observed) != expected_ids:
        missing = sorted(expected_ids.difference(observed))
        extra = sorted(set(observed).difference(expected_ids))
        raise DepthSummaryError(
            f"Depthfile ID set differs from BED IDs: missing={len(missing)} {missing[:5]}; "
            f"extra={len(extra)} {extra[:5]}"
        )
    duplicates = sorted(identifier for identifier, count in observed.items() if count != 1)
    if duplicates:
        raise DepthSummaryError(f"Depthfile IDs do not occur exactly once: {duplicates[:5]}")

    def counts(ids: set[str]) -> tuple[int, int]:
        zero = sum(depths[identifier] == 0 for identifier in ids)
        nonzero = len(ids) - zero
        return zero, nonzero

    reference_zero, reference_nonzero = counts(reference_ids)
    query_zero, query_nonzero = counts(query_ids)
    if (reference_zero, reference_nonzero) != (
        int(expected_reference["zero"]),
        int(expected_reference["nonzero"]),
    ):
        raise DepthSummaryError("Reference depthfile counts disagree with raw-anchor reconstruction")
    if (query_zero, query_nonzero) != (
        int(expected_query["zero"]),
        int(expected_query["nonzero"]),
    ):
        raise DepthSummaryError("Query depthfile counts disagree with raw-anchor reconstruction")
    return {
        "rows": len(observed),
        "reference_zero": reference_zero,
        "reference_nonzero": reference_nonzero,
        "query_zero": query_zero,
        "query_nonzero": query_nonzero,
        "maximum_depth": max(depths.values()),
    }


def atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    path = path.expanduser().resolve()
    if path.exists():
        raise DepthSummaryError(f"Refusing to overwrite existing JSON: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def atomic_write_tsv(path: Path, row: dict[str, object]) -> None:
    path = path.expanduser().resolve()
    if path.exists():
        raise DepthSummaryError(f"Refusing to overwrite existing TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=tuple(row), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerow(row)
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def run(args: argparse.Namespace) -> None:
    if not args.sample or not args.display_name or not args.accession:
        raise DepthSummaryError("Sample, display name, and accession must be non-empty")
    if args.minimum_block_size < 1:
        raise DepthSummaryError("--minimum-block-size must be positive")
    paths = {
        "reference_protein": require_file(args.reference_protein, "Reference protein FASTA"),
        "reference_bed": require_file(args.reference_bed, "Reference BED"),
        "query_protein": require_file(args.query_protein, "Query protein FASTA"),
        "query_bed": require_file(args.query_bed, "Query BED"),
        "anchors": require_file(args.anchors, "Raw JCVI anchors"),
        "depthfile": require_file(args.depthfile, "JCVI depthfile"),
    }
    allowlist_path = getattr(args, "allowed_reference_bed_only_ids", None)
    allowed_reference_bed_only: set[str] = set()
    if allowlist_path is not None:
        resolved_allowlist = require_file(
            allowlist_path, "Reference BED-only identifier allow-list"
        )
        paths["reference_bed_only_allowlist"] = resolved_allowlist
        allowed_reference_bed_only = load_identifier_allowlist(resolved_allowlist)
    reference_fasta_ids = read_fasta_ids(paths["reference_protein"])
    query_fasta_ids = read_fasta_ids(paths["query_protein"])
    reference_bed = load_bed(paths["reference_bed"])
    query_bed = load_bed(paths["query_bed"])
    reference_bed_ids = set(reference_bed.order)
    query_bed_ids = set(query_bed.order)
    if query_bed_ids != query_fasta_ids:
        raise DepthSummaryError("Query BED and protein FASTA ID sets are not identical")
    reference_bed_only = sorted(reference_bed_ids.difference(reference_fasta_ids))
    reference_fasta_only = sorted(reference_fasta_ids.difference(reference_bed_ids))
    observed_reference_bed_only = set(reference_bed_only)
    if reference_fasta_only:
        raise DepthSummaryError(
            "Reference protein FASTA contains IDs absent from the BED: "
            f"count={len(reference_fasta_only)} examples={reference_fasta_only[:5]}"
        )
    if observed_reference_bed_only != allowed_reference_bed_only:
        unexpected = sorted(observed_reference_bed_only.difference(allowed_reference_bed_only))
        absent = sorted(allowed_reference_bed_only.difference(observed_reference_bed_only))
        raise DepthSummaryError(
            "Reference BED/FASTA ID contract failed: BED-only IDs must be empty by "
            "default or exactly match --allowed-reference-bed-only-ids; "
            f"unexpected_bed_only={len(unexpected)} examples={unexpected[:5]}; "
            f"allowed_but_not_bed_only={len(absent)} examples={absent[:5]}"
        )

    blocks, pair_rows = load_anchor_blocks(paths["anchors"])
    small_blocks = [len(block) for block in blocks if len(block) < args.minimum_block_size]
    if small_blocks:
        raise DepthSummaryError(
            f"{len(small_blocks)} raw anchor blocks are below minimum size {args.minimum_block_size}"
        )
    missing_reference_fasta = sorted(
        {pair[0] for block in blocks for pair in block}.difference(reference_fasta_ids)
    )
    missing_query_fasta = sorted(
        {pair[1] for block in blocks for pair in block}.difference(query_fasta_ids)
    )
    if missing_reference_fasta or missing_query_fasta:
        raise DepthSummaryError(
            f"Anchor IDs absent from FASTA: reference={missing_reference_fasta[:5]}; "
            f"query={missing_query_fasta[:5]}"
        )
    reference = summarize_depth(blocks, reference_bed, 0)
    query = summarize_depth(blocks, query_bed, 1)
    depth_validation = validate_depthfile(
        paths["depthfile"], reference_bed_ids, query_bed_ids, reference, query
    )
    metric = (
        "100 * union([min_gene_index,max_gene_index) across raw JCVI anchor blocks) "
        "/ total matching BED rows"
    )
    row: dict[str, object] = {
        "sample": args.sample,
        "display_name": args.display_name,
        "accession": args.accession,
        "anchor_blocks": len(blocks),
        "anchor_pair_rows": pair_rows,
        "minimum_block_size": min(len(block) for block in blocks),
        "reference_bed_rows": reference["total"],
        "reference_zero_depth_gene_indices": reference["zero"],
        "reference_nonzero_depth_gene_indices": reference["nonzero"],
        "reference_coverage_percent": f"{float(reference['coverage']):.6f}",
        "query_bed_rows": query["total"],
        "query_zero_depth_gene_indices": query["zero"],
        "query_nonzero_depth_gene_indices": query["nonzero"],
        "query_coverage_percent": f"{float(query['coverage']):.6f}",
        "reference_bed_only_ids": len(reference_bed_only),
        "reference_bed_fasta_id_identity": not reference_bed_only,
        "reference_bed_only_allowlist_used": allowlist_path is not None,
        "query_bed_fasta_id_identity": True,
        "depthfile_rows": depth_validation["rows"],
        "depthfile_maximum_depth": depth_validation["maximum_depth"],
        "metric_definition": metric,
    }
    payload: dict[str, object] = {
        "schema_version": 1,
        "script": "scripts/qc/summarize_jcvi_depth.py",
        "script_version": SCRIPT_VERSION,
        "summary": row,
        "reference_bed_only_ids": reference_bed_only,
        "allowed_reference_bed_only_ids": sorted(allowed_reference_bed_only),
        "depthfile_validation": depth_validation,
        "inputs": {
            role: {
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for role, path in paths.items()
        },
    }
    atomic_write_tsv(args.output_tsv, row)
    atomic_write_json(args.output_json, payload)


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except (OSError, UnicodeError, ReconstructionError, DepthSummaryError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
