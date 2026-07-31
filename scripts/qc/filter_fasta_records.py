#!/usr/bin/env python3
"""Materialize or filter FASTA records with a complete, atomic audit trail.

The program streams records in input order and never loads sequence bodies into
memory.  Input compression is detected from the gzip magic bytes.  Output
compression defaults to gzip for a ``.gz`` filename and plain text otherwise;
``--output-compression`` can override that rule.  Gzip output is reproducible:
the embedded filename is empty and the modification time is fixed at zero.

A FASTA identifier is the first whitespace-delimited token after ``>``.  IDs
are matched exactly and case-sensitively.  Duplicate input IDs are fatal,
including duplicates among records that would otherwise be excluded.  Every
requested exclusion must be observed unless ``--allow-missing-exclude-id`` is
explicitly supplied.  When no ``--exclude-id`` is supplied, all validated
records are retained; this provides an audited plain/gzip materialization or
transcode for tools that cannot read a compressed FASTA.

The filtered FASTA and JSON audit are first written to temporary files in their
destination directories.  They are published with atomic replacements only
after the complete input and every requested exclusion have been validated.
The implementation uses one process and starts no worker threads.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import sys
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


READ_CHUNK_BYTES = 8 * 1024 * 1024
GZIP_MAGIC = b"\x1f\x8b"
GZIP_COMPRESSLEVEL = 6


class FastaFilterError(RuntimeError):
    """Raised when filtering cannot produce a validated output."""


@dataclass(frozen=True)
class FilterSummary:
    """Counts and logical checksums collected during one streaming pass."""

    input_records: int
    retained_records: int
    excluded_records: int
    input_sequence_characters: int
    retained_sequence_characters: int
    excluded_sequence_characters: int
    observed_excluded_ids: tuple[str, ...]
    excluded_record_lengths: tuple[tuple[str, int], ...]
    input_logical_sha256: str
    output_logical_sha256: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Stream a FASTA file, optionally remove exact record IDs, and "
            "write an atomic JSON audit with counts and SHA-256 checksums."
        )
    )
    parser.add_argument("--input", required=True, type=Path, help="Input FASTA (plain or gzip).")
    parser.add_argument("--output", required=True, type=Path, help="Output FASTA path.")
    parser.add_argument(
        "--exclude-id",
        action="append",
        default=[],
        help=(
            "Exact, case-sensitive FASTA ID to exclude; may be repeated. "
            "The ID is the first token in a FASTA header. If omitted, every "
            "record is retained for an audited materialization/transcode."
        ),
    )
    parser.add_argument(
        "--audit-json",
        required=True,
        type=Path,
        help="Required JSON audit output path.",
    )
    parser.add_argument(
        "--output-compression",
        choices=("auto", "plain", "gzip"),
        default="auto",
        help=(
            "Output encoding. 'auto' (default) uses gzip for a .gz output "
            "filename and plain text otherwise."
        ),
    )
    parser.add_argument(
        "--allow-missing-exclude-id",
        action="store_true",
        help=(
            "Permit requested exclusion IDs that are absent from the input. "
            "Missing IDs remain explicit in the JSON audit."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Atomically replace existing output and audit files.",
    )
    return parser.parse_args()


def absolute_path(path: Path) -> Path:
    """Return an expanded absolute path without dereferencing the final name."""
    expanded = path.expanduser()
    if expanded.is_absolute():
        return expanded
    return Path.cwd() / expanded


def validate_exclusion_ids(values: list[str]) -> tuple[str, ...]:
    """Validate CLI IDs and return them in deterministic lexical order."""
    invalid = [value for value in values if not value or any(char.isspace() for char in value)]
    if invalid:
        display = ", ".join(repr(value) for value in invalid)
        raise FastaFilterError(
            "--exclude-id values must be non-empty single FASTA tokens without whitespace: "
            + display
        )
    if len(values) != len(set(values)):
        duplicates = sorted({value for value in values if values.count(value) > 1})
        raise FastaFilterError(
            "Duplicate --exclude-id value(s): " + ", ".join(duplicates)
        )
    return tuple(sorted(values))


def validate_paths(input_path: Path, output_path: Path, audit_path: Path, force: bool) -> None:
    """Reject missing inputs, path collisions, and unsafe overwrites."""
    if not input_path.exists():
        raise FastaFilterError(f"Input FASTA does not exist: {input_path}")
    if not input_path.is_file():
        raise FastaFilterError(f"Input FASTA is not a regular file: {input_path}")

    resolved = {
        "input": input_path.resolve(strict=True),
        "output": output_path.resolve(strict=False),
        "audit": audit_path.resolve(strict=False),
    }
    if resolved["input"] == resolved["output"]:
        raise FastaFilterError("--output must not overwrite --input")
    if resolved["input"] == resolved["audit"]:
        raise FastaFilterError("--audit-json must not overwrite --input")
    if resolved["output"] == resolved["audit"]:
        raise FastaFilterError("--output and --audit-json must be different paths")

    for label, path in (("output", output_path), ("audit", audit_path)):
        if path.exists() and path.is_dir():
            raise FastaFilterError(f"The {label} path is a directory: {path}")
        if path.exists() and not force:
            raise FastaFilterError(
                f"The {label} path already exists (use --force to replace it): {path}"
            )


def detect_input_compression(path: Path) -> str:
    """Return ``gzip`` when the file starts with the gzip magic bytes."""
    try:
        with path.open("rb") as handle:
            return "gzip" if handle.read(2) == GZIP_MAGIC else "plain"
    except OSError as error:
        raise FastaFilterError(f"Cannot inspect input FASTA {path}: {error}") from error


def choose_output_compression(path: Path, requested: str) -> str:
    """Resolve the explicit or suffix-inferred output compression mode."""
    if requested != "auto":
        return requested
    return "gzip" if path.name.lower().endswith(".gz") else "plain"


@contextmanager
def open_input_text(path: Path, compression: str) -> Iterator[TextIO]:
    """Open a validated UTF-8 FASTA text stream."""
    try:
        if compression == "gzip":
            with gzip.open(
                path,
                "rt",
                encoding="utf-8",
                errors="strict",
                newline="",
            ) as handle:
                yield handle
        else:
            with path.open(
                "rt",
                encoding="utf-8",
                errors="strict",
                newline="",
            ) as handle:
                yield handle
    except (OSError, EOFError, UnicodeError, gzip.BadGzipFile) as error:
        raise FastaFilterError(f"Cannot read FASTA {path}: {error}") from error


@contextmanager
def open_output_text(path: Path, compression: str) -> Iterator[TextIO]:
    """Open a plain or reproducible-gzip UTF-8 output stream."""
    try:
        with path.open("wb") as raw_handle:
            if compression == "gzip":
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    compresslevel=GZIP_COMPRESSLEVEL,
                    fileobj=raw_handle,
                    mtime=0,
                ) as gzip_handle:
                    with io.TextIOWrapper(
                        gzip_handle,
                        encoding="utf-8",
                        errors="strict",
                        newline="",
                    ) as text_handle:
                        yield text_handle
            else:
                with io.TextIOWrapper(
                    raw_handle,
                    encoding="utf-8",
                    errors="strict",
                    newline="",
                ) as text_handle:
                    yield text_handle
    except (OSError, UnicodeError) as error:
        raise FastaFilterError(f"Cannot write temporary FASTA {path}: {error}") from error


def make_temporary_path(target: Path) -> Path:
    """Create and return an empty temporary path beside its final target."""
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor, name = tempfile.mkstemp(
            prefix=f".{target.name}.",
            suffix=".tmp",
            dir=target.parent,
        )
    except OSError as error:
        raise FastaFilterError(
            f"Cannot create a temporary file beside {target}: {error}"
        ) from error
    os.close(descriptor)
    return Path(name)


def sha256_file(path: Path) -> tuple[int, str]:
    """Return physical file size and SHA-256 from a bounded-memory pass."""
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            while True:
                chunk = handle.read(READ_CHUNK_BYTES)
                if not chunk:
                    break
                size += len(chunk)
                digest.update(chunk)
    except OSError as error:
        raise FastaFilterError(f"Cannot checksum {path}: {error}") from error
    return size, digest.hexdigest()


def input_stat_signature(path: Path) -> tuple[int, int, int, int]:
    """Return fields used to detect a file changed during filtering."""
    try:
        status = path.stat()
    except OSError as error:
        raise FastaFilterError(f"Cannot stat input FASTA {path}: {error}") from error
    return status.st_dev, status.st_ino, status.st_size, status.st_mtime_ns


def stream_filter(
    input_handle: TextIO,
    output_handle: TextIO,
    excluded_ids: set[str],
    source_name: str,
) -> FilterSummary:
    """Validate and filter FASTA records while preserving retained text."""
    seen_ids: set[str] = set()
    observed_excluded_ids: set[str] = set()
    excluded_record_lengths: dict[str, int] = {}
    input_logical_digest = hashlib.sha256()
    output_logical_digest = hashlib.sha256()

    input_records = 0
    retained_records = 0
    excluded_records = 0
    input_sequence_characters = 0
    retained_sequence_characters = 0
    excluded_sequence_characters = 0

    current_id: str | None = None
    current_is_excluded = False
    current_length = 0

    def finish_record() -> None:
        if current_id is not None and current_is_excluded:
            excluded_record_lengths[current_id] = current_length

    def write_retained(line: str) -> None:
        output_handle.write(line)
        output_logical_digest.update(line.encode("utf-8"))

    for line_number, line in enumerate(input_handle, start=1):
        encoded_line = line.encode("utf-8")
        input_logical_digest.update(encoded_line)

        if line.startswith(">"):
            finish_record()
            header_tokens = line[1:].split()
            if not header_tokens:
                raise FastaFilterError(
                    f"{source_name}: empty FASTA identifier at line {line_number}"
                )
            record_id = header_tokens[0]
            if record_id in seen_ids:
                raise FastaFilterError(
                    f"{source_name}: duplicate FASTA identifier {record_id!r} "
                    f"at line {line_number}"
                )
            seen_ids.add(record_id)

            current_id = record_id
            current_is_excluded = record_id in excluded_ids
            current_length = 0
            input_records += 1
            if current_is_excluded:
                observed_excluded_ids.add(record_id)
                excluded_records += 1
            else:
                retained_records += 1
                write_retained(line)
            continue

        if current_id is None:
            if not line.strip():
                # In zero-exclusion materialization mode, preserve harmless
                # leading whitespace so the logical decompressed input and
                # output text checksums remain identical.  Filtering mode
                # omits it because it belongs to no retained FASTA record.
                if not excluded_ids:
                    write_retained(line)
                continue
            raise FastaFilterError(
                f"{source_name}: sequence data before the first FASTA header "
                f"at line {line_number}"
            )

        sequence_characters = len("".join(line.split()))
        current_length += sequence_characters
        input_sequence_characters += sequence_characters
        if current_is_excluded:
            excluded_sequence_characters += sequence_characters
        else:
            retained_sequence_characters += sequence_characters
            write_retained(line)

    finish_record()
    if input_records == 0:
        raise FastaFilterError(f"{source_name}: FASTA contains no records")
    if retained_records + excluded_records != input_records:
        raise AssertionError("FASTA record accounting is inconsistent")
    if retained_sequence_characters + excluded_sequence_characters != input_sequence_characters:
        raise AssertionError("FASTA sequence-length accounting is inconsistent")

    return FilterSummary(
        input_records=input_records,
        retained_records=retained_records,
        excluded_records=excluded_records,
        input_sequence_characters=input_sequence_characters,
        retained_sequence_characters=retained_sequence_characters,
        excluded_sequence_characters=excluded_sequence_characters,
        observed_excluded_ids=tuple(sorted(observed_excluded_ids)),
        excluded_record_lengths=tuple(sorted(excluded_record_lengths.items())),
        input_logical_sha256=input_logical_digest.hexdigest(),
        output_logical_sha256=output_logical_digest.hexdigest(),
    )


def build_audit(
    *,
    input_path: Path,
    output_path: Path,
    input_compression: str,
    output_compression: str,
    input_size: int,
    input_sha256: str,
    output_size: int,
    output_sha256: str,
    requested_ids: tuple[str, ...],
    missing_ids: tuple[str, ...],
    allow_missing: bool,
    summary: FilterSummary,
) -> dict[str, object]:
    """Construct the stable, English-only JSON audit document."""
    return {
        "schema_version": "1.0",
        "status": "completed",
        "operation": (
            "exclude_exact_fasta_record_ids" if requested_ids else "materialize_fasta"
        ),
        "rules": {
            "fasta_identifier": "first whitespace-delimited token after '>'",
            "identifier_matching": "exact and case-sensitive",
            "duplicate_identifier_policy": (
                "fatal for all input records, including excluded records"
            ),
            "sequence_character_count": "all non-whitespace characters on sequence lines",
            "record_order": "retained records preserve input order and text",
            "input_compression_detection": "gzip magic bytes",
            "gzip_output": (
                f"compresslevel={GZIP_COMPRESSLEVEL}; empty embedded filename; mtime=0"
            ),
        },
        "input": {
            "path": str(input_path),
            "compression": input_compression,
            "size_bytes": input_size,
            "sha256": input_sha256,
            "logical_fasta_text_sha256": summary.input_logical_sha256,
        },
        "output": {
            "path": str(output_path),
            "compression": output_compression,
            "size_bytes": output_size,
            "sha256": output_sha256,
            "logical_fasta_text_sha256": summary.output_logical_sha256,
        },
        "exclusions": {
            "requested_ids": list(requested_ids),
            "observed_ids": list(summary.observed_excluded_ids),
            "missing_ids": list(missing_ids),
            "allow_missing_requested_ids": allow_missing,
            "excluded_record_lengths": [
                {"id": record_id, "sequence_characters": length}
                for record_id, length in summary.excluded_record_lengths
            ],
        },
        "counts": {
            "input_records": summary.input_records,
            "retained_records": summary.retained_records,
            "excluded_records": summary.excluded_records,
            "input_sequence_characters": summary.input_sequence_characters,
            "retained_sequence_characters": summary.retained_sequence_characters,
            "excluded_sequence_characters": summary.excluded_sequence_characters,
        },
    }


def write_audit(path: Path, audit: dict[str, object]) -> None:
    """Write a deterministic JSON document to a temporary path."""
    try:
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(audit, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
    except OSError as error:
        raise FastaFilterError(f"Cannot write JSON audit {path}: {error}") from error


def run(args: argparse.Namespace) -> dict[str, object]:
    """Execute filtering and atomically publish the FASTA and audit."""
    input_path = absolute_path(args.input)
    output_path = absolute_path(args.output)
    audit_path = absolute_path(args.audit_json)
    requested_ids = validate_exclusion_ids(args.exclude_id)
    validate_paths(input_path, output_path, audit_path, args.force)

    input_compression = detect_input_compression(input_path)
    output_compression = choose_output_compression(output_path, args.output_compression)
    input_signature = input_stat_signature(input_path)
    input_size, input_sha256 = sha256_file(input_path)

    output_temp: Path | None = None
    audit_temp: Path | None = None
    try:
        output_temp = make_temporary_path(output_path)
        with open_input_text(input_path, input_compression) as input_handle:
            with open_output_text(output_temp, output_compression) as output_handle:
                summary = stream_filter(
                    input_handle,
                    output_handle,
                    set(requested_ids),
                    str(input_path),
                )

        if input_stat_signature(input_path) != input_signature:
            raise FastaFilterError(
                f"Input FASTA changed while it was being filtered: {input_path}"
            )

        missing_ids = tuple(sorted(set(requested_ids) - set(summary.observed_excluded_ids)))
        if missing_ids and not args.allow_missing_exclude_id:
            raise FastaFilterError(
                "Requested exclusion ID(s) were not found; no output was published: "
                + ", ".join(missing_ids)
                + ". Use --allow-missing-exclude-id only when this is intentional."
            )

        output_size, output_sha256 = sha256_file(output_temp)
        audit = build_audit(
            input_path=input_path,
            output_path=output_path,
            input_compression=input_compression,
            output_compression=output_compression,
            input_size=input_size,
            input_sha256=input_sha256,
            output_size=output_size,
            output_sha256=output_sha256,
            requested_ids=requested_ids,
            missing_ids=missing_ids,
            allow_missing=args.allow_missing_exclude_id,
            summary=summary,
        )

        audit_temp = make_temporary_path(audit_path)
        write_audit(audit_temp, audit)

        # Recheck immediately before publication to avoid silently replacing a
        # file created during the streaming pass when --force was not given.
        if not args.force:
            if output_path.exists():
                raise FastaFilterError(f"Output appeared during the run: {output_path}")
            if audit_path.exists():
                raise FastaFilterError(f"Audit output appeared during the run: {audit_path}")

        os.replace(output_temp, output_path)
        output_temp = None
        os.replace(audit_temp, audit_path)
        audit_temp = None
        return audit
    finally:
        for temporary_path in (output_temp, audit_temp):
            if temporary_path is not None:
                try:
                    temporary_path.unlink(missing_ok=True)
                except OSError:
                    pass


def main() -> int:
    args = parse_args()
    try:
        audit = run(args)
    except FastaFilterError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    except OSError as error:
        print(f"ERROR: filesystem operation failed: {error}", file=sys.stderr)
        return 2

    counts = audit["counts"]
    assert isinstance(counts, dict)
    print(
        "Completed FASTA filtering: "
        f"retained={counts['retained_records']}, "
        f"excluded={counts['excluded_records']}, "
        f"output={absolute_path(args.output)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
