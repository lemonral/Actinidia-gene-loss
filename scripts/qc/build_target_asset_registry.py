#!/usr/bin/env python3
"""Create one exact genome/GFF/protein registry for chromosome assignment."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys
from pathlib import Path


COLUMNS = (
    "assembly_unit_id",
    "target_scope_id",
    "asset_role",
    "file_name",
    "bytes",
    "sha256",
    "status",
)


class RegistryError(RuntimeError):
    pass


def stable_binding(path: Path) -> dict[str, str | int]:
    source = path.expanduser().resolve()
    if not source.is_file() or source.is_symlink():
        raise RegistryError(f"Asset must be a regular non-symlink file: {source}")
    before = source.stat()
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    after = source.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise RegistryError(f"Asset changed while hashing: {source.name}")
    if after.st_size <= 0:
        raise RegistryError(f"Asset is empty: {source.name}")
    return {"file_name": source.name, "bytes": after.st_size, "sha256": digest.hexdigest()}


def build_rows(
    *,
    assembly_unit_id: str,
    target_scope_id: str,
    genome: Path,
    gff: Path,
    protein: Path,
) -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    for role, path in (("genome", genome), ("gff", gff), ("protein", protein)):
        rows.append(
            {
                "assembly_unit_id": assembly_unit_id,
                "target_scope_id": target_scope_id,
                "asset_role": role,
                **stable_binding(path),
                "status": "verified",
            }
        )
    if len({row["file_name"] for row in rows}) != 3:
        raise RegistryError("Asset roles must use three distinct basenames")
    if len({row["sha256"] for row in rows}) != 3:
        raise RegistryError("Asset roles must not contain copy-identical files")
    return rows


def run(args: argparse.Namespace) -> Path:
    output = args.output.expanduser().resolve()
    if output.exists():
        raise RegistryError(f"Refusing to overwrite existing registry: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = build_rows(
        assembly_unit_id=args.assembly_unit_id,
        target_scope_id=args.target_scope_id,
        genome=args.genome,
        gff=args.gff,
        protein=args.protein,
    )
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.link(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    return output


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assembly-unit-id", required=True)
    p.add_argument("--target-scope-id", required=True)
    p.add_argument("--genome", required=True, type=Path)
    p.add_argument("--gff", required=True, type=Path)
    p.add_argument("--protein", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    return p


def main() -> int:
    try:
        output = run(parser().parse_args())
        print(f"PASS\t{output}")
        return 0
    except (OSError, RegistryError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
