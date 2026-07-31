#!/usr/bin/env python3
"""Collect one-row SynOrths summary TSVs into a deterministic table.

Every input must be a one-row TSV produced with the same summary schema.  The
collector requires the exact same ordered header and ``schema_version`` in all
files, and requires a unique, non-empty ``sample`` value in each row.  It does
not merge partially compatible tables or silently normalize sample names.

Inputs can be passed by repeating ``--summary-tsv`` or through
``--summary-list``.  A list contains one path per non-empty line; lines whose
first non-whitespace character is ``#`` are comments, and relative paths are
resolved from the list file's directory.  The collected rows are sorted by
sample.  The output contains no timestamp and is written atomically.

Examples
--------
python scripts/assembly_qc/collect_synorth_summaries.py \
  --summary-tsv results/synorth/current/current.synorth_summary.tsv \
  --summary-tsv results/synorth/candidate/candidate.synorth_summary.tsv \
  --output results/synorth/synorth_summary.tsv

python scripts/assembly_qc/collect_synorth_summaries.py \
  --summary-list config/synorth_summary_paths.txt \
  --output results/synorth/synorth_summary.tsv
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


class SynorthCollectionError(RuntimeError):
    """Raised when summary inputs do not form one auditable table."""


@dataclass(frozen=True)
class SummaryRow:
    """One validated summary row and its source path."""

    source: Path
    header: tuple[str, ...]
    values: tuple[str, ...]

    def value(self, column: str) -> str:
        """Return one value using the already validated unique header."""
        return self.values[self.header.index(column)]


def require_nonempty_file(raw_path: Path, label: str) -> Path:
    """Resolve and require a non-empty regular file."""
    path = raw_path.expanduser().resolve()
    if not path.is_file():
        raise SynorthCollectionError(f"{label} is not a regular file: {path}")
    try:
        size = path.stat().st_size
    except OSError as error:
        raise SynorthCollectionError(f"Cannot inspect {label} {path}: {error}") from error
    if size == 0:
        raise SynorthCollectionError(f"{label} is empty: {path}")
    return path


def paths_from_list(raw_path: Path) -> tuple[Path, list[Path]]:
    """Read one-path-per-line input, resolving relative paths beside the list."""
    list_path = require_nonempty_file(raw_path, "Summary list")
    paths: list[Path] = []
    try:
        with list_path.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, start=1):
                value = raw.strip()
                if not value or value.startswith("#"):
                    continue
                candidate = Path(value).expanduser()
                if not candidate.is_absolute():
                    candidate = list_path.parent / candidate
                paths.append(candidate)
    except (OSError, UnicodeError) as error:
        raise SynorthCollectionError(f"Cannot read summary list {list_path}: {error}") from error
    if not paths:
        raise SynorthCollectionError(f"Summary list contains no input paths: {list_path}")
    return list_path, paths


def read_summary(raw_path: Path) -> SummaryRow:
    """Read exactly one data row from one summary TSV."""
    path = require_nonempty_file(raw_path, "SynOrths summary TSV")
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
    except (OSError, UnicodeError, csv.Error) as error:
        raise SynorthCollectionError(f"Cannot read summary TSV {path}: {error}") from error

    if not rows or not rows[0]:
        raise SynorthCollectionError(f"Summary TSV has no header: {path}")
    header = tuple(rows[0])
    if any(not column for column in header):
        raise SynorthCollectionError(f"Summary TSV has an empty header field: {path}")
    if len(set(header)) != len(header):
        raise SynorthCollectionError(f"Summary TSV has duplicate header fields: {path}")
    if len(rows) != 2:
        raise SynorthCollectionError(
            f"Summary TSV must contain exactly one data row; found {max(0, len(rows) - 1)}: {path}"
        )
    values = tuple(rows[1])
    if len(values) != len(header):
        raise SynorthCollectionError(
            f"Summary TSV row has {len(values)} fields but the header has {len(header)}: {path}"
        )
    for required in ("sample", "schema_version"):
        if required not in header:
            raise SynorthCollectionError(
                f"Summary TSV is missing required column {required!r}: {path}"
            )
    row = SummaryRow(source=path, header=header, values=values)
    sample = row.value("sample")
    if not sample.strip():
        raise SynorthCollectionError(f"Summary TSV has an empty sample value: {path}")
    if sample != sample.strip():
        raise SynorthCollectionError(
            f"Summary TSV sample has leading or trailing whitespace: {path}"
        )
    if not row.value("schema_version").strip():
        raise SynorthCollectionError(f"Summary TSV has an empty schema_version: {path}")
    return row


def collect(raw_paths: Iterable[Path]) -> tuple[tuple[str, ...], list[SummaryRow]]:
    """Validate a common schema and return rows sorted by unique sample."""
    rows = [read_summary(path) for path in raw_paths]
    if not rows:
        raise SynorthCollectionError("No SynOrths summary TSV inputs were provided")

    expected_header = rows[0].header
    expected_schema = rows[0].value("schema_version")
    sample_sources: dict[str, Path] = {}
    for row in rows:
        if row.header != expected_header:
            raise SynorthCollectionError(
                "Summary TSV header/schema mismatch: "
                f"{row.source} does not have the exact ordered header of {rows[0].source}"
            )
        schema = row.value("schema_version")
        if schema != expected_schema:
            raise SynorthCollectionError(
                "Summary TSV schema_version mismatch: "
                f"{row.source} has {schema!r}, expected {expected_schema!r} from {rows[0].source}"
            )
        sample = row.value("sample")
        if sample in sample_sources:
            raise SynorthCollectionError(
                f"Duplicate sample {sample!r} in {sample_sources[sample]} and {row.source}"
            )
        sample_sources[sample] = row.source

    return expected_header, sorted(rows, key=lambda row: row.value("sample"))


def atomic_write_tsv(path: Path, header: tuple[str, ...], rows: Iterable[SummaryRow]) -> None:
    """Write the collected TSV atomically without changing its schema."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for row in rows:
                writer.writerow(row.values)
        os.replace(temporary, path)
    except OSError as error:
        temporary.unlink(missing_ok=True)
        raise SynorthCollectionError(f"Cannot write collected TSV {path}: {error}") from error


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "--summary-tsv",
        action="append",
        type=Path,
        help="One summary TSV; repeat this option for multiple inputs",
    )
    inputs.add_argument(
        "--summary-list",
        type=Path,
        help="Text file containing one summary TSV path per line",
    )
    parser.add_argument("--output", required=True, type=Path, help="Collected output TSV")
    return parser


def run(args: argparse.Namespace) -> int:
    """Resolve inputs, validate their schema, and write the collected table."""
    list_path: Path | None = None
    if args.summary_list is not None:
        list_path, raw_paths = paths_from_list(args.summary_list)
    else:
        raw_paths = list(args.summary_tsv or [])

    resolved_paths = [
        require_nonempty_file(path, "SynOrths summary TSV") for path in raw_paths
    ]
    output = args.output.expanduser().resolve()
    protected = set(resolved_paths)
    if list_path is not None:
        protected.add(list_path)
    if output in protected:
        raise SynorthCollectionError(f"Output path would overwrite an input file: {output}")

    header, rows = collect(resolved_paths)
    atomic_write_tsv(output, header, rows)
    print(f"Collected {len(rows)} SynOrths summaries into {output}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except SynorthCollectionError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
