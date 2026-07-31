#!/usr/bin/env python3
"""Verify downloaded external assets against the maintained TSV manifest.

The verifier uses only the Python standard library and reads each present file
in fixed-size binary chunks.  Paths are built from ``local_subdirectory`` and
``local_filename`` beneath ``--data-root``; paths that escape that root are
rejected before any file is opened.

Exit codes
----------
0
    Every listed file is a regular file with the expected byte size and
    SHA-256 digest.
1
    At least one file is missing, unreadable, not a regular file, or differs
    from its expected size or digest.  The complete status TSV is still
    written.
2
    The manifest, command-line input, or output path is invalid.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, TextIO


REQUIRED_COLUMNS = (
    "dataset",
    "accession",
    "asset_type",
    "local_subdirectory",
    "local_filename",
    "size_bytes",
    "sha256",
    "source_url",
)

OUTPUT_COLUMNS = (
    "manifest_row",
    "dataset",
    "accession",
    "asset_type",
    "relative_path",
    "expected_size_bytes",
    "actual_size_bytes",
    "expected_sha256",
    "actual_sha256",
    "status",
    "source_url",
)

SHA256_PATTERN = re.compile(r"[0-9a-fA-F]{64}")
READ_CHUNK_BYTES = 8 * 1024 * 1024


class VerificationInputError(RuntimeError):
    """Raised when verification cannot start from a valid manifest."""


@dataclass(frozen=True)
class ManifestEntry:
    """One validated manifest row and its resolved local path."""

    manifest_row: int
    dataset: str
    accession: str
    asset_type: str
    relative_path: str
    path: Path
    expected_size: int
    expected_sha256: str
    source_url: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Verify external download sizes and SHA-256 digests and write a "
            "deterministic TSV audit report."
        )
    )
    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="Maintained tab-separated download manifest.",
    )
    parser.add_argument(
        "--data-root",
        required=True,
        type=Path,
        help="Directory containing the manifest's local subdirectories.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output status TSV path, or '-' for standard output.",
    )
    return parser.parse_args()


def require_value(row: dict[str, str | None], column: str, row_number: int) -> str:
    """Return one stripped, non-empty value or raise a clear input error."""
    value = row.get(column)
    if value is None or not value.strip():
        raise VerificationInputError(
            f"Manifest row {row_number} has an empty required value: {column}"
        )
    return value.strip()


def resolve_under_root(data_root: Path, subdirectory: str, filename: str) -> tuple[Path, str]:
    """Resolve an asset path and prove that it remains beneath data_root."""
    subdirectory_path = Path(subdirectory)
    filename_path = Path(filename)
    if subdirectory_path.is_absolute() or filename_path.is_absolute():
        raise VerificationInputError(
            "local_subdirectory and local_filename must both be relative paths"
        )

    candidate = (data_root / subdirectory_path / filename_path).resolve(strict=False)
    try:
        relative = candidate.relative_to(data_root)
    except ValueError as error:
        raise VerificationInputError(
            f"Resolved asset path escapes --data-root: {subdirectory}/{filename}"
        ) from error
    if relative == Path("."):
        raise VerificationInputError("An asset path cannot resolve to --data-root itself")
    return candidate, relative.as_posix()


def load_manifest(manifest: Path, data_root: Path) -> list[ManifestEntry]:
    """Load and validate every maintained-schema manifest row."""
    try:
        handle = manifest.open("r", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise VerificationInputError(f"Cannot open manifest {manifest}: {error}") from error

    root = data_root.expanduser().resolve(strict=False)
    entries: list[ManifestEntry] = []
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise VerificationInputError(f"Manifest has no header: {manifest}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise VerificationInputError("Manifest contains duplicate column names")
        missing_columns = [name for name in REQUIRED_COLUMNS if name not in reader.fieldnames]
        if missing_columns:
            raise VerificationInputError(
                "Manifest is missing required columns: " + ", ".join(missing_columns)
            )

        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise VerificationInputError(
                    f"Manifest row {row_number} contains more fields than the header"
                )

            subdirectory = require_value(row, "local_subdirectory", row_number)
            filename = require_value(row, "local_filename", row_number)
            size_text = require_value(row, "size_bytes", row_number)
            digest = require_value(row, "sha256", row_number).lower()

            try:
                expected_size = int(size_text)
            except ValueError as error:
                raise VerificationInputError(
                    f"Manifest row {row_number} has a non-integer size_bytes: {size_text}"
                ) from error
            if expected_size < 0:
                raise VerificationInputError(
                    f"Manifest row {row_number} has a negative size_bytes: {size_text}"
                )
            if SHA256_PATTERN.fullmatch(digest) is None:
                raise VerificationInputError(
                    f"Manifest row {row_number} sha256 must contain exactly 64 hexadecimal characters"
                )

            path, relative_path = resolve_under_root(root, subdirectory, filename)
            entries.append(
                ManifestEntry(
                    manifest_row=row_number,
                    dataset=(row.get("dataset") or "").strip(),
                    accession=(row.get("accession") or "").strip(),
                    asset_type=(row.get("asset_type") or "").strip(),
                    relative_path=relative_path,
                    path=path,
                    expected_size=expected_size,
                    expected_sha256=digest,
                    source_url=(row.get("source_url") or "").strip(),
                )
            )

    if not entries:
        raise VerificationInputError(f"Manifest contains no data rows: {manifest}")
    return entries


def stream_size_and_sha256(handle: BinaryIO) -> tuple[int, str]:
    """Return the exact byte count and SHA-256 digest from one pass."""
    digest = hashlib.sha256()
    size = 0
    while True:
        chunk = handle.read(READ_CHUNK_BYTES)
        if not chunk:
            break
        size += len(chunk)
        digest.update(chunk)
    return size, digest.hexdigest()


def verify_entry(entry: ManifestEntry) -> dict[str, str | int]:
    """Verify one entry and return one output row."""
    actual_size: int | str = ""
    actual_sha256 = ""

    if not entry.path.exists():
        status = "missing"
    elif not entry.path.is_file():
        status = "not_regular_file"
    else:
        try:
            with entry.path.open("rb") as handle:
                actual_size, actual_sha256 = stream_size_and_sha256(handle)
        except OSError:
            status = "unreadable"
        else:
            size_matches = actual_size == entry.expected_size
            hash_matches = actual_sha256 == entry.expected_sha256
            if size_matches and hash_matches:
                status = "ok"
            elif not size_matches and not hash_matches:
                status = "size_and_sha256_mismatch"
            elif not size_matches:
                status = "size_mismatch"
            else:
                status = "sha256_mismatch"

    return {
        "manifest_row": entry.manifest_row,
        "dataset": entry.dataset,
        "accession": entry.accession,
        "asset_type": entry.asset_type,
        "relative_path": entry.relative_path,
        "expected_size_bytes": entry.expected_size,
        "actual_size_bytes": actual_size,
        "expected_sha256": entry.expected_sha256,
        "actual_sha256": actual_sha256,
        "status": status,
        "source_url": entry.source_url,
    }


def open_output(output: str) -> tuple[TextIO, bool]:
    """Open an output stream and report whether the caller must close it."""
    if output == "-":
        return sys.stdout, False

    path = Path(output).expanduser()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        return path.open("w", encoding="utf-8", newline=""), True
    except OSError as error:
        raise VerificationInputError(f"Cannot open output {path}: {error}") from error


def main() -> int:
    args = parse_args()
    try:
        entries = load_manifest(args.manifest.expanduser(), args.data_root)
        rows = [verify_entry(entry) for entry in entries]
        output_handle, must_close = open_output(args.output)
        try:
            writer = csv.DictWriter(
                output_handle,
                fieldnames=OUTPUT_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        finally:
            if must_close:
                output_handle.close()
    except VerificationInputError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    return 0 if all(row["status"] == "ok" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
