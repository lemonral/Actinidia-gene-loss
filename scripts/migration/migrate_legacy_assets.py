#!/usr/bin/env python3
"""Migrate legacy analysis assets into a named external data layout.

Large files are linked and small files are copied.  The old absolute source
paths are read only from a private runtime manifest; they never need to be
committed to the repository.  A checksum report records what was migrated.
"""

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
from typing import Dict, Iterable, Mapping


REQUIRED_MAP = {
    "legacy_sample",
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "accession",
    "source_url",
}
REQUIRED_LEGACY = {"sample", "genome", "gff", "protein"}
ROLES = ("genome", "gff", "protein")


class MigrationError(RuntimeError):
    """Raised when a migration cannot be proven complete and unambiguous."""


def read_tsv(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = set(required) - fields
        if missing:
            raise MigrationError(f"{path}: missing columns {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def unique_index(rows: Iterable[Mapping[str, str]], key: str, label: str) -> Dict[str, Mapping[str, str]]:
    result: Dict[str, Mapping[str, str]] = {}
    for row in rows:
        value = row[key]
        if not value:
            raise MigrationError(f"{label}: empty {key}")
        if value in result:
            raise MigrationError(f"{label}: duplicate {key} {value!r}")
        result[value] = row
    return result


def sha256(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_name_suffix(source: Path) -> str:
    suffixes = "".join(source.suffixes)
    return suffixes if suffixes else ".dat"


def migrate_one(source: Path, destination: Path, copy_limit: int, dry_run: bool) -> str:
    if not source.is_file():
        raise MigrationError(f"missing source file: {source}")
    mode = "copy" if source.stat().st_size <= copy_limit else "symlink"

    if destination.exists() or destination.is_symlink():
        if destination.is_symlink() and mode == "symlink":
            if destination.resolve(strict=True) != source.resolve(strict=True):
                raise MigrationError(f"existing link points elsewhere: {destination}")
            return mode
        if destination.is_file() and mode == "copy":
            if destination.stat().st_size != source.stat().st_size:
                raise MigrationError(f"existing copy has wrong size: {destination}")
            return mode
        raise MigrationError(f"existing destination has wrong type: {destination}")

    if dry_run:
        return mode
    destination.parent.mkdir(parents=True, exist_ok=True)
    if mode == "copy":
        shutil.copy2(source, destination)
    else:
        destination.symlink_to(source.resolve(strict=True))
    return mode


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--legacy-manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument(
        "--resolved-manifest",
        type=Path,
        help=(
            "Optional compatibility TSV for basic-stats/BUSCO consumers. "
            "This runtime file belongs in the external data store, not Git."
        ),
    )
    parser.add_argument(
        "--copy-max-mib",
        type=float,
        default=25.0,
        help="Copy files at or below this size; soft-link larger files (default: 25 MiB).",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def atomic_write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise MigrationError("cannot write an empty resolved legacy manifest")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(rows[0]),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    args = parse_args()
    try:
        mapping_rows = read_tsv(args.mapping, REQUIRED_MAP)
        legacy_rows = read_tsv(args.legacy_manifest, REQUIRED_LEGACY)
        mapping = unique_index(mapping_rows, "legacy_sample", str(args.mapping))
        legacy = unique_index(legacy_rows, "sample", str(args.legacy_manifest))
        if set(mapping) != set(legacy):
            raise MigrationError(
                "mapping and legacy sample sets differ: "
                f"mapping_only={sorted(set(mapping) - set(legacy))}; "
                f"legacy_only={sorted(set(legacy) - set(mapping))}"
            )
        unit_ids = [row["assembly_unit_id"] for row in mapping_rows]
        if len(unit_ids) != len(set(unit_ids)):
            raise MigrationError("mapping contains duplicate assembly_unit_id values")
        if args.copy_max_mib < 0:
            raise MigrationError("--copy-max-mib must be non-negative")

        copy_limit = int(args.copy_max_mib * 1024 * 1024)
        results = []
        resolved_units = []
        for sample in sorted(mapping):
            declared = mapping[sample]
            old = legacy[sample]
            unit_paths: dict[str, str] = {}
            for role in ROLES:
                source = Path(old[role]).expanduser()
                suffix = safe_name_suffix(source)
                relative = Path("legacy_linked") / declared["assembly_unit_id"] / f"{role}{suffix}"
                destination = (args.data_root / relative).resolve(strict=False)
                expected_root = args.data_root.resolve(strict=False)
                if os.path.commonpath((str(expected_root), str(destination))) != str(expected_root):
                    raise MigrationError(f"destination escapes data root: {relative}")
                mode = migrate_one(source, destination, copy_limit, args.dry_run)
                size = source.stat().st_size
                unit_paths[role] = str(destination)
                results.append(
                    {
                        "legacy_sample": sample,
                        "assembly_unit_id": declared["assembly_unit_id"],
                        "biological_species": declared["biological_species"],
                        "haplotype_or_subgenome": declared["haplotype_or_subgenome"],
                        "accession": declared["accession"],
                        "source_url": declared["source_url"],
                        "role": role,
                        "relative_path": relative.as_posix(),
                        "mode": mode,
                        "size_bytes": size,
                        "sha256": sha256(source),
                    }
                )
            resolved_units.append(
                {
                    "assembly_unit_id": declared["assembly_unit_id"],
                    "biological_species": declared["biological_species"],
                    "individual_id": declared.get("individual_id", "legacy_study_individual"),
                    "haplotype_or_subgenome": declared["haplotype_or_subgenome"],
                    "ploidy": declared.get("ploidy", "not_recorded"),
                    "accession": declared["accession"],
                    "legacy_status": declared.get("legacy_status", "legacy_input"),
                    "sample": declared["assembly_unit_id"],
                    "species": declared["biological_species"],
                    "target_haplotype": declared["assembly_unit_id"],
                    "current_or_alternative": "legacy_manuscript_reproduction",
                    "include_downstream": "false",
                    "genome": unit_paths["genome"],
                    "gff": unit_paths["gff"],
                    "protein": unit_paths["protein"],
                    "source_url": declared["source_url"],
                    "notes": declared.get("notes", ""),
                }
            )

        payload = {
            "schema_version": 1,
            "status": "dry_run" if args.dry_run else "complete",
            "mapping": str(args.mapping),
            "legacy_manifest": str(args.legacy_manifest),
            "data_root": str(args.data_root),
            "copy_max_bytes": copy_limit,
            "sample_count": len(mapping),
            "asset_count": len(results),
            "assets": results,
        }
        if not args.dry_run:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            temporary = args.report.with_suffix(args.report.suffix + ".tmp")
            with temporary.open("w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2, sort_keys=True)
                handle.write("\n")
            temporary.replace(args.report)
            if args.resolved_manifest is not None:
                atomic_write_tsv(args.resolved_manifest.expanduser().resolve(), resolved_units)
        print(json.dumps({key: payload[key] for key in ("status", "sample_count", "asset_count")}, sort_keys=True))
        return 0
    except (MigrationError, OSError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
