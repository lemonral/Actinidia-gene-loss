#!/usr/bin/env python3
"""Audit one SynOrths pair table against its protein and coordinate inputs.

The archived historical SynOrths 1.5 pair files used by this project are
headerless, whitespace-separated tables whose target-query IDs are in column
5 and C. scandens reference IDs are in column 1.  The defaults preserve that
archived-file convention: ``--query-column 5 --reference-column 1``.  Column
meaning is not intrinsic to the words "query" and "reference": it follows the
biological inputs supplied to SynOrths.  New standardized reviewer runs pass
the target query as ``-a`` and C. scandens as ``-b`` and therefore explicitly
use ``--query-column 1 --reference-column 5``.  Always verify the producer
command and input IDs before selecting the 1-based columns.

Coverage conventions
--------------------
``query_coverage_percent`` and ``reference_coverage_percent`` use all unique
FASTA IDs as their denominators and count only anchors present in the matching
FASTA.  Separate non-empty-FASTA and coordinate coverage fields are emitted so
that empty protein records, missing identifiers, or coordinate inconsistencies
cannot silently change the primary denominator.  Exact repeated query/reference
pairs count as duplicate pair rows even if other SynOrths columns differ.

FASTA inputs may be plain text or gzip-compressed.  Pair and coordinate tables
also accept either form.  Blank lines and full-line comments beginning with
``#`` are ignored.  The JSON and one-row TSV contain no timestamps and use
stable field/list ordering, making repeated summaries deterministic.

Example
-------
python scripts/assembly_qc/summarize_synorth.py \
  --sample target_assembly \
  --pairs results/synorth/target_vs_reference.synorths.txt \
  --query-fasta data/target.primary.faa \
  --reference-fasta data/reference.primary.faa \
  --query-coords data/target.coords \
  --reference-coords data/reference.coords \
  --output-json results/synorth/target_vs_reference.synorth_summary.json \
  --output-tsv results/synorth/target_vs_reference.synorth_summary.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO


SCHEMA_VERSION = 2
SAFE_SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


class SynorthSummaryError(RuntimeError):
    """Raised when an input cannot produce an auditable SynOrths summary."""


@dataclass(frozen=True)
class FastaInventory:
    """Unique and record-level identifiers from one protein FASTA."""

    record_count: int
    ids: frozenset[str]
    nonempty_ids: frozenset[str]
    duplicate_id_records: int


@dataclass(frozen=True)
class IdTableInventory:
    """Identifiers read from one coordinate-style table."""

    row_count: int
    ids: frozenset[str]
    duplicate_id_rows: int


@dataclass(frozen=True)
class PairInventory:
    """Selected query/reference identifiers from a SynOrths pair table."""

    row_count: int
    pairs: Counter[tuple[str, str]]
    query_ids: frozenset[str]
    reference_ids: frozenset[str]


def open_text(path: Path) -> TextIO:
    """Open a plain or ``.gz`` text file and report compressed-read errors."""
    try:
        if path.suffix.lower() == ".gz":
            return gzip.open(path, "rt", encoding="utf-8", errors="replace")
        return path.open("rt", encoding="utf-8", errors="replace")
    except OSError as error:
        raise SynorthSummaryError(f"Cannot open {path}: {error}") from error


def require_input(path: Path, label: str) -> Path:
    """Resolve and validate a required non-empty regular file."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise SynorthSummaryError(f"{label} is not a regular file: {resolved}")
    try:
        if resolved.stat().st_size == 0:
            raise SynorthSummaryError(f"{label} is empty: {resolved}")
    except OSError as error:
        raise SynorthSummaryError(f"Cannot inspect {label} {resolved}: {error}") from error
    return resolved


def validate_sample(sample: str) -> str:
    """Require a non-empty, path-safe sample identifier."""
    if not SAFE_SAMPLE_RE.fullmatch(sample):
        raise SynorthSummaryError(
            "--sample must be a safe non-empty identifier: 1-128 characters, "
            "starting with an ASCII letter or digit and containing only ASCII "
            "letters, digits, '.', '_', or '-'"
        )
    return sample


def read_fasta(path: Path) -> FastaInventory:
    """Read FASTA identifiers while retaining only per-ID non-empty state."""
    identifiers: set[str] = set()
    nonempty_identifiers: set[str] = set()
    record_count = 0
    current_id: str | None = None
    current_nonempty = False

    def finish_record() -> None:
        if current_id is not None and current_nonempty:
            nonempty_identifiers.add(current_id)

    try:
        with open_text(path) as handle:
            for line_number, raw in enumerate(handle, start=1):
                stripped = raw.strip()
                if not stripped:
                    continue
                if stripped.startswith(">"):
                    finish_record()
                    header = stripped[1:].strip()
                    if not header:
                        raise SynorthSummaryError(
                            f"{path}:{line_number}: FASTA header has no identifier"
                        )
                    current_id = header.split()[0]
                    if not current_id:
                        raise SynorthSummaryError(
                            f"{path}:{line_number}: FASTA header has no identifier"
                        )
                    identifiers.add(current_id)
                    record_count += 1
                    current_nonempty = False
                    continue
                if current_id is None:
                    raise SynorthSummaryError(
                        f"{path}:{line_number}: sequence data precedes the first FASTA header"
                    )
                sequence = "".join(stripped.split())
                if sequence.replace("-", "").replace(".", ""):
                    current_nonempty = True
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise SynorthSummaryError(f"Cannot read FASTA {path}: {error}") from error

    finish_record()
    if record_count == 0:
        raise SynorthSummaryError(f"{path}: no FASTA records")
    return FastaInventory(
        record_count=record_count,
        ids=frozenset(identifiers),
        nonempty_ids=frozenset(nonempty_identifiers),
        duplicate_id_records=record_count - len(identifiers),
    )


def validate_column(column: int, option: str) -> None:
    """Require a positive 1-based table column."""
    if column < 1:
        raise SynorthSummaryError(f"{option} must be a positive 1-based column number")


def read_id_table(
    path: Path,
    *,
    id_column: int,
    has_header: bool,
    label: str,
) -> IdTableInventory:
    """Read unique IDs from a whitespace-separated coordinate table."""
    validate_column(id_column, f"--{label}-coords-id-column")
    identifiers: set[str] = set()
    row_count = 0
    header_skipped = not has_header

    try:
        with open_text(path) as handle:
            for line_number, raw in enumerate(handle, start=1):
                stripped = raw.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                fields = stripped.split()
                if len(fields) < id_column:
                    kind = "header" if not header_skipped else "row"
                    raise SynorthSummaryError(
                        f"{path}:{line_number}: {label} coordinate {kind} needs at least "
                        f"{id_column} columns for ID column {id_column}; found {len(fields)}"
                    )
                if not header_skipped:
                    header_skipped = True
                    continue
                identifier = fields[id_column - 1]
                if not identifier:
                    raise SynorthSummaryError(
                        f"{path}:{line_number}: empty {label} coordinate identifier"
                    )
                identifiers.add(identifier)
                row_count += 1
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise SynorthSummaryError(f"Cannot read coordinate table {path}: {error}") from error

    if has_header and not header_skipped:
        raise SynorthSummaryError(f"{path}: no non-comment coordinate header")
    if row_count == 0:
        raise SynorthSummaryError(f"{path}: no {label} coordinate records")
    return IdTableInventory(
        row_count=row_count,
        ids=frozenset(identifiers),
        duplicate_id_rows=row_count - len(identifiers),
    )


def read_pairs(
    path: Path,
    *,
    query_column: int,
    reference_column: int,
    has_header: bool,
) -> PairInventory:
    """Read selected 1-based columns from a whitespace-separated pair file."""
    validate_column(query_column, "--query-column")
    validate_column(reference_column, "--reference-column")
    if query_column == reference_column:
        raise SynorthSummaryError(
            "--query-column and --reference-column must select different columns"
        )

    required_columns = max(query_column, reference_column)
    pairs: Counter[tuple[str, str]] = Counter()
    row_count = 0
    header_skipped = not has_header

    try:
        with open_text(path) as handle:
            for line_number, raw in enumerate(handle, start=1):
                stripped = raw.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                fields = stripped.split()
                if len(fields) < required_columns:
                    kind = "header" if not header_skipped else "row"
                    raise SynorthSummaryError(
                        f"{path}:{line_number}: SynOrths {kind} needs at least "
                        f"{required_columns} columns for query/reference columns "
                        f"{query_column}/{reference_column}; found {len(fields)}"
                    )
                if not header_skipped:
                    header_skipped = True
                    continue
                query_id = fields[query_column - 1]
                reference_id = fields[reference_column - 1]
                if not query_id or not reference_id:
                    raise SynorthSummaryError(
                        f"{path}:{line_number}: selected query/reference identifier is empty"
                    )
                pairs[(query_id, reference_id)] += 1
                row_count += 1
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise SynorthSummaryError(f"Cannot read SynOrths pair table {path}: {error}") from error

    if has_header and not header_skipped:
        raise SynorthSummaryError(f"{path}: no non-comment SynOrths header")
    if row_count == 0:
        raise SynorthSummaryError(f"{path}: no SynOrths pair rows")
    return PairInventory(
        row_count=row_count,
        pairs=pairs,
        query_ids=frozenset(query for query, _ in pairs),
        reference_ids=frozenset(reference for _, reference in pairs),
    )


def percentage(numerator: int, denominator: int) -> float | None:
    """Return a stable six-decimal percentage, or null for a zero universe."""
    if denominator == 0:
        return None
    return round(100.0 * numerator / denominator, 6)


def sorted_ids(values: set[str] | frozenset[str]) -> list[str]:
    """Return identifiers in deterministic lexical order."""
    return sorted(values)


def validate_pair_columns(
    pairs: PairInventory,
    query_fasta: FastaInventory,
    reference_fasta: FastaInventory,
    query_column: int,
    reference_column: int,
) -> None:
    """Reject selected columns that have no support in their declared FASTA."""
    query_matches = pairs.query_ids & query_fasta.ids
    reference_matches = pairs.reference_ids & reference_fasta.ids
    if query_matches and reference_matches:
        return
    cross_query = len(pairs.query_ids & reference_fasta.ids)
    cross_reference = len(pairs.reference_ids & query_fasta.ids)
    problems: list[str] = []
    if not query_matches:
        problems.append(
            f"column {query_column} has 0/{len(pairs.query_ids)} unique IDs in the query FASTA"
        )
    if not reference_matches:
        problems.append(
            f"column {reference_column} has 0/{len(pairs.reference_ids)} unique IDs in the reference FASTA"
        )
    raise SynorthSummaryError(
        "Selected SynOrths ID columns failed FASTA validation: "
        + "; ".join(problems)
        + f". Cross-side matches are query-to-reference={cross_query} and "
        f"reference-to-query={cross_reference}. Verify --query-column, "
        "--reference-column, the SynOrths -a/-b order, and identifier normalization."
    )


def build_summary(
    *,
    sample: str,
    pair_path: Path,
    query_fasta_path: Path,
    reference_fasta_path: Path,
    query_coords_path: Path | None,
    reference_coords_path: Path | None,
    query_column: int,
    reference_column: int,
    query_coords_id_column: int,
    reference_coords_id_column: int,
    pairs_has_header: bool,
    query_coords_has_header: bool,
    reference_coords_has_header: bool,
) -> dict[str, Any]:
    """Read inputs and construct the stable JSON/TSV payload."""
    query_fasta = read_fasta(query_fasta_path)
    reference_fasta = read_fasta(reference_fasta_path)
    pairs = read_pairs(
        pair_path,
        query_column=query_column,
        reference_column=reference_column,
        has_header=pairs_has_header,
    )
    validate_pair_columns(
        pairs, query_fasta, reference_fasta, query_column, reference_column
    )

    query_coords = (
        read_id_table(
            query_coords_path,
            id_column=query_coords_id_column,
            has_header=query_coords_has_header,
            label="query",
        )
        if query_coords_path is not None
        else None
    )
    reference_coords = (
        read_id_table(
            reference_coords_path,
            id_column=reference_coords_id_column,
            has_header=reference_coords_has_header,
            label="reference",
        )
        if reference_coords_path is not None
        else None
    )

    query_in_fasta = pairs.query_ids & query_fasta.ids
    reference_in_fasta = pairs.reference_ids & reference_fasta.ids
    query_in_nonempty_fasta = pairs.query_ids & query_fasta.nonempty_ids
    reference_in_nonempty_fasta = pairs.reference_ids & reference_fasta.nonempty_ids
    query_in_coords = pairs.query_ids & query_coords.ids if query_coords else None
    reference_in_coords = pairs.reference_ids & reference_coords.ids if reference_coords else None

    metrics: dict[str, int | float | None] = {
        "query_fasta_records": query_fasta.record_count,
        "query_fasta_unique_ids": len(query_fasta.ids),
        "query_fasta_nonempty_ids": len(query_fasta.nonempty_ids),
        "query_fasta_empty_ids": len(query_fasta.ids - query_fasta.nonempty_ids),
        "query_fasta_duplicate_id_records": query_fasta.duplicate_id_records,
        "reference_fasta_records": reference_fasta.record_count,
        "reference_fasta_unique_ids": len(reference_fasta.ids),
        "reference_fasta_nonempty_ids": len(reference_fasta.nonempty_ids),
        "reference_fasta_empty_ids": len(reference_fasta.ids - reference_fasta.nonempty_ids),
        "reference_fasta_duplicate_id_records": reference_fasta.duplicate_id_records,
        "query_coordinate_rows": query_coords.row_count if query_coords else None,
        "query_coordinate_unique_ids": len(query_coords.ids) if query_coords else None,
        "query_coordinate_duplicate_id_rows": (
            query_coords.duplicate_id_rows if query_coords else None
        ),
        "reference_coordinate_rows": reference_coords.row_count if reference_coords else None,
        "reference_coordinate_unique_ids": len(reference_coords.ids) if reference_coords else None,
        "reference_coordinate_duplicate_id_rows": (
            reference_coords.duplicate_id_rows if reference_coords else None
        ),
        "pair_rows": pairs.row_count,
        "unique_pairs": len(pairs.pairs),
        "duplicate_pair_rows": pairs.row_count - len(pairs.pairs),
        "unique_query_anchors": len(pairs.query_ids),
        "unique_reference_anchors": len(pairs.reference_ids),
        "query_anchors_in_fasta": len(query_in_fasta),
        "reference_anchors_in_fasta": len(reference_in_fasta),
        "query_anchors_in_nonempty_fasta": len(query_in_nonempty_fasta),
        "reference_anchors_in_nonempty_fasta": len(reference_in_nonempty_fasta),
        "query_anchors_in_coordinates": len(query_in_coords) if query_in_coords is not None else None,
        "reference_anchors_in_coordinates": (
            len(reference_in_coords) if reference_in_coords is not None else None
        ),
        "query_coverage_percent": percentage(len(query_in_fasta), len(query_fasta.ids)),
        "reference_coverage_percent": percentage(
            len(reference_in_fasta), len(reference_fasta.ids)
        ),
        "query_nonempty_coverage_percent": percentage(
            len(query_in_nonempty_fasta), len(query_fasta.nonempty_ids)
        ),
        "reference_nonempty_coverage_percent": percentage(
            len(reference_in_nonempty_fasta), len(reference_fasta.nonempty_ids)
        ),
        "query_coordinate_coverage_percent": (
            percentage(len(query_in_coords), len(query_coords.ids)) if query_coords else None
        ),
        "reference_coordinate_coverage_percent": (
            percentage(len(reference_in_coords), len(reference_coords.ids))
            if reference_coords
            else None
        ),
    }

    absent_ids: dict[str, list[str] | None] = {
        "query_anchor_ids_absent_from_fasta": sorted_ids(pairs.query_ids - query_fasta.ids),
        "reference_anchor_ids_absent_from_fasta": sorted_ids(
            pairs.reference_ids - reference_fasta.ids
        ),
        "query_anchor_ids_absent_from_coordinates": (
            sorted_ids(pairs.query_ids - query_coords.ids) if query_coords else None
        ),
        "reference_anchor_ids_absent_from_coordinates": (
            sorted_ids(pairs.reference_ids - reference_coords.ids) if reference_coords else None
        ),
        "query_coordinate_ids_absent_from_fasta": (
            sorted_ids(query_coords.ids - query_fasta.ids) if query_coords else None
        ),
        "reference_coordinate_ids_absent_from_fasta": (
            sorted_ids(reference_coords.ids - reference_fasta.ids) if reference_coords else None
        ),
        "query_fasta_ids_absent_from_coordinates": (
            sorted_ids(query_fasta.ids - query_coords.ids) if query_coords else None
        ),
        "reference_fasta_ids_absent_from_coordinates": (
            sorted_ids(reference_fasta.ids - reference_coords.ids) if reference_coords else None
        ),
    }
    absent_counts = {
        f"{name}_count": len(values) if values is not None else None
        for name, values in absent_ids.items()
    }

    return {
        "sample": sample,
        "schema_version": SCHEMA_VERSION,
        "inputs": {
            "pairs": str(pair_path),
            "query_fasta": str(query_fasta_path),
            "reference_fasta": str(reference_fasta_path),
            "query_coords": str(query_coords_path) if query_coords_path else None,
            "reference_coords": str(reference_coords_path) if reference_coords_path else None,
            "query_column_1_based": query_column,
            "reference_column_1_based": reference_column,
            "query_coords_id_column_1_based": (
                query_coords_id_column if query_coords_path else None
            ),
            "reference_coords_id_column_1_based": (
                reference_coords_id_column if reference_coords_path else None
            ),
            "pairs_has_header": pairs_has_header,
            "query_coords_has_header": query_coords_has_header if query_coords_path else None,
            "reference_coords_has_header": (
                reference_coords_has_header if reference_coords_path else None
            ),
        },
        "metrics": {**metrics, **absent_counts},
        "absent_ids": absent_ids,
    }


def scalar_string(value: Any) -> str:
    """Format one TSV scalar without platform-dependent float rendering."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def flatten_for_tsv(payload: dict[str, Any]) -> dict[str, str]:
    """Flatten the JSON sections into one stable, self-describing TSV row."""
    row: dict[str, str] = {
        "sample": str(payload["sample"]),
        "schema_version": str(payload["schema_version"]),
    }
    for name, value in payload["inputs"].items():
        row[name] = scalar_string(value)
    for name, value in payload["metrics"].items():
        row[name] = scalar_string(value)
    for name, values in payload["absent_ids"].items():
        row[name] = (
            "" if values is None else json.dumps(values, ensure_ascii=False, separators=(",", ":"))
        )
    return row


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write sorted, timestamp-free JSON atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except OSError as error:
        temporary.unlink(missing_ok=True)
        raise SynorthSummaryError(f"Cannot write JSON output {path}: {error}") from error


def atomic_write_tsv(path: Path, payload: dict[str, Any]) -> None:
    """Write one summary row with stable insertion-ordered columns atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    row = flatten_for_tsv(payload)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=tuple(row), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerow(row)
        os.replace(temporary, path)
    except OSError as error:
        temporary.unlink(missing_ok=True)
        raise SynorthSummaryError(f"Cannot write TSV output {path}: {error}") from error


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample",
        required=True,
        help=(
            "Stable sample identifier (1-128 safe characters; starts with an "
            "ASCII letter or digit)"
        ),
    )
    parser.add_argument(
        "--pairs",
        "--synorth-pairs",
        dest="pairs",
        required=True,
        type=Path,
        help="SynOrths pair/result table (plain or .gz)",
    )
    parser.add_argument("--query-fasta", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--query-coords", type=Path)
    parser.add_argument("--reference-coords", type=Path)
    parser.add_argument(
        "--query-column", type=int, default=5, help="1-based query ID column (default: 5)"
    )
    parser.add_argument(
        "--reference-column",
        type=int,
        default=1,
        help="1-based reference ID column (default: 1)",
    )
    parser.add_argument(
        "--query-coords-id-column",
        type=int,
        default=1,
        help="1-based query coordinate ID column (default: 1)",
    )
    parser.add_argument(
        "--reference-coords-id-column",
        type=int,
        default=1,
        help="1-based reference coordinate ID column (default: 1)",
    )
    parser.add_argument(
        "--pairs-has-header", action="store_true", help="Skip the first non-comment pair row"
    )
    parser.add_argument(
        "--query-coords-has-header",
        action="store_true",
        help="Skip the first non-comment query coordinate row",
    )
    parser.add_argument(
        "--reference-coords-has-header",
        action="store_true",
        help="Skip the first non-comment reference coordinate row",
    )
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    return parser


def run(args: argparse.Namespace) -> int:
    """Validate paths, calculate the summary, and write both output formats."""
    sample = validate_sample(args.sample)
    pair_path = require_input(args.pairs, "SynOrths pair table")
    query_fasta_path = require_input(args.query_fasta, "Query protein FASTA")
    reference_fasta_path = require_input(args.reference_fasta, "Reference protein FASTA")
    query_coords_path = (
        require_input(args.query_coords, "Query coordinate table") if args.query_coords else None
    )
    reference_coords_path = (
        require_input(args.reference_coords, "Reference coordinate table")
        if args.reference_coords
        else None
    )
    output_json = args.output_json.expanduser().resolve()
    output_tsv = args.output_tsv.expanduser().resolve()
    if output_json == output_tsv:
        raise SynorthSummaryError("--output-json and --output-tsv must be different paths")
    input_paths = {
        pair_path,
        query_fasta_path,
        reference_fasta_path,
        *([query_coords_path] if query_coords_path else []),
        *([reference_coords_path] if reference_coords_path else []),
    }
    for output in (output_json, output_tsv):
        if output in input_paths:
            raise SynorthSummaryError(f"Output path would overwrite an input file: {output}")

    payload = build_summary(
        sample=sample,
        pair_path=pair_path,
        query_fasta_path=query_fasta_path,
        reference_fasta_path=reference_fasta_path,
        query_coords_path=query_coords_path,
        reference_coords_path=reference_coords_path,
        query_column=args.query_column,
        reference_column=args.reference_column,
        query_coords_id_column=args.query_coords_id_column,
        reference_coords_id_column=args.reference_coords_id_column,
        pairs_has_header=args.pairs_has_header,
        query_coords_has_header=args.query_coords_has_header,
        reference_coords_has_header=args.reference_coords_has_header,
    )
    atomic_write_json(output_json, payload)
    atomic_write_tsv(output_tsv, payload)
    print(
        f"Summarized {payload['metrics']['pair_rows']} SynOrths pair rows into "
        f"{output_json} and {output_tsv}"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except SynorthSummaryError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
