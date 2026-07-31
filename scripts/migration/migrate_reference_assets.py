#!/usr/bin/env python3
"""Migrate a checksum-frozen legacy reference bundle into the data store.

The committed mapping contains only paths relative to a user-supplied legacy
root. Small files are copied and large files are soft-linked. Every source is
verified against its declared SHA-256 before any destination is accepted. The
runtime report and optional resolved manifest belong in the external data
store because they contain machine-specific absolute paths.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import sys
import tempfile
from typing import Iterable


REQUIRED_COLUMNS = (
    "reference_id",
    "biological_species",
    "role",
    "legacy_relative_path",
    "data_relative_path",
    "expected_sha256",
    "status",
    "notes",
)
SHA256_LENGTH = 64


class ReferenceMigrationError(RuntimeError):
    """Raised when a reference asset cannot be migrated without ambiguity."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_relative_path(value: str, field: str) -> Path:
    pure = PurePosixPath(value)
    if not value or pure.is_absolute() or ".." in pure.parts or "." in pure.parts:
        raise ReferenceMigrationError(f"{field} must be a normalized relative POSIX path: {value!r}")
    return Path(*pure.parts)


def read_mapping(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        missing = sorted(set(REQUIRED_COLUMNS).difference(fields))
        if missing:
            raise ReferenceMigrationError(f"{path}: missing columns {missing}")
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
            if row and any((value or "").strip() for value in row.values())
        ]
    if not rows:
        raise ReferenceMigrationError(f"{path}: mapping contains no assets")
    seen_keys: set[tuple[str, str]] = set()
    seen_destinations: set[str] = set()
    for row_number, row in enumerate(rows, 2):
        if not row["reference_id"] or not row["biological_species"] or not row["role"]:
            raise ReferenceMigrationError(f"{path}:{row_number}: identity fields cannot be empty")
        key = (row["reference_id"], row["role"])
        if key in seen_keys:
            raise ReferenceMigrationError(f"{path}:{row_number}: duplicate reference/role {key}")
        seen_keys.add(key)
        safe_relative_path(row["legacy_relative_path"], "legacy_relative_path")
        destination = safe_relative_path(row["data_relative_path"], "data_relative_path").as_posix()
        if destination in seen_destinations:
            raise ReferenceMigrationError(f"{path}:{row_number}: duplicate destination {destination}")
        seen_destinations.add(destination)
        checksum = row["expected_sha256"].lower()
        if len(checksum) != SHA256_LENGTH or any(char not in "0123456789abcdef" for char in checksum):
            raise ReferenceMigrationError(f"{path}:{row_number}: invalid expected_sha256")
        row["expected_sha256"] = checksum
    return rows


def ensure_below(root: Path, relative: Path, field: str) -> Path:
    root = root.expanduser().resolve()
    # Keep this check lexical.  Resolving an existing destination symlink would
    # replace the data-store path with its legacy target and make an otherwise
    # valid, idempotent rerun appear to escape ``data_root``.
    candidate = root / relative
    if os.path.commonpath((str(root), str(candidate))) != str(root):
        raise ReferenceMigrationError(f"{field} escapes its declared root: {relative}")
    return candidate


def migrate_one(source: Path, destination: Path, mode: str, expected_sha256: str, dry_run: bool) -> None:
    if destination.exists() or destination.is_symlink():
        if mode == "symlink":
            if not destination.is_symlink() or destination.resolve(strict=True) != source.resolve(strict=True):
                raise ReferenceMigrationError(f"existing destination is not the expected link: {destination}")
        elif destination.is_symlink() or not destination.is_file():
            raise ReferenceMigrationError(f"existing destination is not the expected copy: {destination}")
        if sha256(destination) != expected_sha256:
            raise ReferenceMigrationError(f"existing destination checksum mismatch: {destination}")
        return
    if dry_run:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    if mode == "copy":
        temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
        try:
            shutil.copy2(source, temporary)
            if sha256(temporary) != expected_sha256:
                raise ReferenceMigrationError(f"temporary copy checksum mismatch: {temporary}")
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
    else:
        temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}")
        try:
            temporary.symlink_to(source.resolve(strict=True))
            os.replace(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def write_resolved_manifest(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    fields = (
        "reference_id",
        "biological_species",
        "role",
        "path",
        "relative_path",
        "mode",
        "size_bytes",
        "sha256",
        "status",
        "notes",
    )
    output = []
    output.append("\t".join(fields))
    for row in rows:
        values = [str(row[field]).replace("\t", " ").replace("\n", " ") for field in fields]
        output.append("\t".join(values))
    atomic_write(path, "\n".join(output) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--legacy-root", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--resolved-manifest", type=Path)
    parser.add_argument("--copy-max-mib", type=float, default=25.0)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def run(args: argparse.Namespace) -> dict[str, object]:
    if args.copy_max_mib < 0:
        raise ReferenceMigrationError("--copy-max-mib must be nonnegative")
    mapping_path = args.mapping.expanduser().resolve()
    legacy_root = args.legacy_root.expanduser().resolve()
    data_root = args.data_root.expanduser().resolve()
    if not legacy_root.is_dir():
        raise ReferenceMigrationError(f"legacy root is not a directory: {legacy_root}")
    data_root.mkdir(parents=True, exist_ok=True)
    rows = read_mapping(mapping_path)
    copy_limit = int(args.copy_max_mib * 1024 * 1024)
    results: list[dict[str, object]] = []
    for row in rows:
        source = ensure_below(
            legacy_root,
            safe_relative_path(row["legacy_relative_path"], "legacy_relative_path"),
            "legacy_relative_path",
        )
        if not source.is_file():
            raise ReferenceMigrationError(f"source is missing or not a file: {source}")
        observed = sha256(source)
        if observed != row["expected_sha256"]:
            raise ReferenceMigrationError(
                f"source checksum mismatch for {row['reference_id']} {row['role']}: "
                f"expected {row['expected_sha256']}, observed {observed}"
            )
        relative_destination = safe_relative_path(row["data_relative_path"], "data_relative_path")
        destination = ensure_below(data_root, relative_destination, "data_relative_path")
        mode = "copy" if source.stat().st_size <= copy_limit else "symlink"
        migrate_one(source, destination, mode, observed, args.dry_run)
        results.append(
            {
                "reference_id": row["reference_id"],
                "biological_species": row["biological_species"],
                "role": row["role"],
                "source_path": str(source),
                "path": str(destination),
                "relative_path": relative_destination.as_posix(),
                "mode": mode,
                "size_bytes": source.stat().st_size,
                "sha256": observed,
                "status": row["status"],
                "notes": row["notes"],
            }
        )
    payload: dict[str, object] = {
        "schema_version": 1,
        "status": "dry_run" if args.dry_run else "complete",
        "mapping_path": str(mapping_path),
        "legacy_root": str(legacy_root),
        "data_root": str(data_root),
        "copy_max_bytes": copy_limit,
        "reference_count": len({str(row["reference_id"]) for row in results}),
        "asset_count": len(results),
        "assets": results,
    }
    if not args.dry_run:
        atomic_write(args.report.expanduser().resolve(), json.dumps(payload, indent=2, sort_keys=True) + "\n")
        if args.resolved_manifest is not None:
            write_resolved_manifest(args.resolved_manifest.expanduser().resolve(), results)
    return payload


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        payload = run(args)
    except (OSError, csv.Error, ReferenceMigrationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {key: payload[key] for key in ("status", "reference_count", "asset_count")},
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
