#!/usr/bin/env python3
"""Convert a validated protein-coordinate pair to deterministic JCVI BED6.

Input coordinates are headerless, tab-separated rows in the SynOrths project
contract::

    gene_id  seqid  start_1_based  end_1_based_inclusive  strand

Output is sorted JCVI BED6::

    seqid  start_0_based  end_0_based_exclusive  gene_id  0  strand

The program requires exact equality between protein FASTA IDs and coordinate
IDs.  It refuses duplicate IDs, invalid coordinates, unexpected strands,
checksum mismatches, or an existing output file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import tempfile
from pathlib import Path


SCRIPT_VERSION = "1.0.0"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class BedPreparationError(RuntimeError):
    """Raised when input integrity or BED conversion is unsafe."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein", required=True, type=Path)
    parser.add_argument("--coords", required=True, type=Path)
    parser.add_argument("--output-bed", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--expected-protein-sha256", required=True)
    parser.add_argument("--expected-coords-sha256", required=True)
    return parser.parse_args()


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise BedPreparationError(f"{label} is missing or empty: {resolved}")
    return resolved


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_digest(path: Path, expected: str, label: str) -> str:
    expected = expected.strip().lower()
    if not SHA256_RE.fullmatch(expected):
        raise BedPreparationError(f"{label} expected SHA-256 is malformed: {expected!r}")
    observed = sha256(path)
    if observed != expected:
        raise BedPreparationError(
            f"{label} SHA-256 mismatch: observed {observed}; expected {expected}; file {path}"
        )
    return observed


def read_fasta_ids(path: Path) -> set[str]:
    identifiers: set[str] = set()
    records = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            records += 1
            identifier = raw[1:].strip().split(None, 1)[0] if raw[1:].strip() else ""
            if not identifier:
                raise BedPreparationError(f"{path}:{line_number}: empty FASTA identifier")
            if identifier in identifiers:
                raise BedPreparationError(f"{path}:{line_number}: duplicate FASTA ID {identifier!r}")
            identifiers.add(identifier)
    if records == 0:
        raise BedPreparationError(f"Protein FASTA contains no records: {path}")
    return identifiers


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value))


def read_coords(path: Path) -> tuple[list[tuple[str, int, int, str, str]], set[str]]:
    rows: list[tuple[str, int, int, str, str]] = []
    identifiers: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                raise BedPreparationError(f"{path}:{line_number}: blank coordinate row")
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise BedPreparationError(
                    f"{path}:{line_number}: expected exactly five coordinate fields, found {len(fields)}"
                )
            identifier, seqid, start_text, end_text, strand = fields
            if not identifier or not seqid:
                raise BedPreparationError(f"{path}:{line_number}: empty gene ID or sequence ID")
            if identifier in identifiers:
                raise BedPreparationError(f"{path}:{line_number}: duplicate coordinate ID {identifier!r}")
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as exc:
                raise BedPreparationError(f"{path}:{line_number}: non-integer coordinate") from exc
            if start < 1 or end < start:
                raise BedPreparationError(f"{path}:{line_number}: invalid one-based inclusive interval")
            if strand not in {"+", "-"}:
                raise BedPreparationError(f"{path}:{line_number}: strand must be '+' or '-'")
            identifiers.add(identifier)
            rows.append((seqid, start - 1, end, identifier, strand))
    if not rows:
        raise BedPreparationError(f"Coordinate table contains no rows: {path}")
    rows.sort(key=lambda row: (natural_key(row[0]), row[1], row[3]))
    return rows, identifiers


def write_bed(path: Path, rows: list[tuple[str, int, int, str, str]]) -> None:
    if path.exists():
        raise BedPreparationError(f"Refusing to overwrite existing BED: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            for seqid, start, end, identifier, strand in rows:
                handle.write(f"{seqid}\t{start}\t{end}\t{identifier}\t0\t{strand}\n")
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def write_json(path: Path, payload: dict[str, object]) -> None:
    if path.exists():
        raise BedPreparationError(f"Refusing to overwrite existing summary JSON: {path}")
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


def run(args: argparse.Namespace) -> None:
    protein = require_file(args.protein, "Protein FASTA")
    coords = require_file(args.coords, "Coordinate table")
    output_bed = args.output_bed.expanduser().resolve()
    summary_json = args.summary_json.expanduser().resolve()
    if output_bed == summary_json:
        raise BedPreparationError("BED and summary JSON paths must differ")
    protein_digest = require_digest(
        protein, args.expected_protein_sha256, "Protein FASTA"
    )
    coords_digest = require_digest(coords, args.expected_coords_sha256, "Coordinate table")
    protein_ids = read_fasta_ids(protein)
    rows, coordinate_ids = read_coords(coords)
    only_protein = sorted(protein_ids.difference(coordinate_ids))
    only_coords = sorted(coordinate_ids.difference(protein_ids))
    if only_protein or only_coords:
        raise BedPreparationError(
            "Protein/coordinate ID sets differ: "
            f"protein_only={len(only_protein)} examples={only_protein[:5]}; "
            f"coordinate_only={len(only_coords)} examples={only_coords[:5]}"
        )
    write_bed(output_bed, rows)
    payload = {
        "schema_version": 1,
        "script": "scripts/assembly_qc/prepare_jcvi_bed.py",
        "script_version": SCRIPT_VERSION,
        "protein": {
            "path": str(protein),
            "records": len(protein_ids),
            "size_bytes": protein.stat().st_size,
            "sha256": protein_digest,
        },
        "coordinates": {
            "path": str(coords),
            "rows": len(rows),
            "size_bytes": coords.stat().st_size,
            "sha256": coords_digest,
            "input_coordinate_system": "one-based inclusive",
        },
        "bed": {
            "path": str(output_bed),
            "rows": len(rows),
            "size_bytes": output_bed.stat().st_size,
            "sha256": sha256(output_bed),
            "coordinate_system": "zero-based half-open BED6",
            "sort_key": "natural(seqid),start,gene_id",
        },
        "id_set_identity": True,
    }
    write_json(summary_json, payload)


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except (OSError, UnicodeError, BedPreparationError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
