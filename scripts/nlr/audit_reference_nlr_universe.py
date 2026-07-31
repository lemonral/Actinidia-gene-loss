#!/usr/bin/env python3
"""Audit the frozen C. scandens sequence universe used for NLR loss rates.

The historical loss-sequence extractor first required a reference ID to be
present in the protein FASTA and then extracted the same ID from the nucleotide
FASTA.  This audit reconstructs that callable universe without running
NLR-Annotator.  It also verifies that the maintained primary-CDS FASTA is an
exact sequence-level representation of the callable historical nucleotide
records.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import platform
import shlex
import shutil
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class FastaRecord:
    length: int
    sequence_sha256: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-coords", type=Path, required=True)
    parser.add_argument("--reference-proteins", type=Path, required=True)
    parser.add_argument("--reference-nucleotides", type=Path, required=True)
    parser.add_argument("--reference-cds", type=Path, required=True)
    parser.add_argument("--expected-coordinate-records", type=int, default=35558)
    parser.add_argument("--expected-callable-records", type=int, default=35547)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_coordinate_ids(path: Path) -> list[str]:
    ids: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            values = line.rstrip("\r\n").split("\t")
            if len(values) < 5:
                raise ValueError(f"Expected at least five coordinate columns at {path}:{line_no}")
            record_id = values[0].strip()
            if not record_id:
                raise ValueError(f"Empty coordinate ID at {path}:{line_no}")
            if record_id in seen:
                raise ValueError(f"Duplicate coordinate ID {record_id!r} in {path}")
            seen.add(record_id)
            ids.append(record_id)
    if not ids:
        raise ValueError(f"No coordinate records found in {path}")
    return ids


def read_fasta_catalog(path: Path) -> dict[str, FastaRecord]:
    catalog: dict[str, FastaRecord] = {}
    current_id: str | None = None
    current_length = 0
    current_digest = hashlib.sha256()

    def finish_record() -> None:
        nonlocal current_id, current_length, current_digest
        if current_id is None:
            return
        if current_length == 0:
            raise ValueError(f"FASTA record {current_id!r} has an empty sequence in {path}")
        catalog[current_id] = FastaRecord(current_length, current_digest.hexdigest())

    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_no}")
                current_id = header.split()[0]
                if current_id in catalog:
                    raise ValueError(f"Duplicate FASTA ID {current_id!r} in {path}")
                current_length = 0
                current_digest = hashlib.sha256()
                continue
            if current_id is None:
                if line.strip():
                    raise ValueError(f"Sequence data precedes the first FASTA header at {path}:{line_no}")
                continue
            sequence = "".join(line.split()).upper()
            if sequence:
                encoded = sequence.encode("ascii")
                current_digest.update(encoded)
                current_length += len(encoded)
    finish_record()
    if not catalog:
        raise ValueError(f"No FASTA records found in {path}")
    return catalog


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"ERROR: Output directory already exists: {args.output_dir}")
    if args.expected_coordinate_records < 1 or args.expected_callable_records < 1:
        raise SystemExit("ERROR: Expected record counts must be positive")
    inputs = {
        "reference_coords": args.reference_coords,
        "reference_proteins": args.reference_proteins,
        "reference_nucleotides": args.reference_nucleotides,
        "reference_cds": args.reference_cds,
    }
    for path in inputs.values():
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: Missing or empty input: {path}")

    try:
        coordinate_ids = read_coordinate_ids(args.reference_coords)
        proteins = read_fasta_catalog(args.reference_proteins)
        nucleotides = read_fasta_catalog(args.reference_nucleotides)
        cds = read_fasta_catalog(args.reference_cds)
        coordinate_set = set(coordinate_ids)

        if len(coordinate_ids) != args.expected_coordinate_records:
            raise ValueError(
                f"Coordinate record count is {len(coordinate_ids)}, "
                f"expected {args.expected_coordinate_records}"
            )
        if set(proteins) - coordinate_set:
            raise ValueError(
                f"Protein FASTA has {len(set(proteins) - coordinate_set)} IDs outside the coordinate universe"
            )
        if set(cds) != set(proteins):
            raise ValueError(
                "Primary-CDS and historical protein ID sets differ: "
                f"only_in_cds={len(set(cds) - set(proteins))}, "
                f"only_in_proteins={len(set(proteins) - set(cds))}"
            )
        missing_nucleotide = coordinate_set - set(nucleotides)
        if missing_nucleotide:
            raise ValueError(
                f"Historical nucleotide FASTA lacks {len(missing_nucleotide)} coordinate IDs"
            )

        callable_ids = coordinate_set & set(proteins) & set(cds) & set(nucleotides)
        if len(callable_ids) != args.expected_callable_records:
            raise ValueError(
                f"Callable reference record count is {len(callable_ids)}, "
                f"expected {args.expected_callable_records}"
            )
        sequence_mismatches = [
            record_id for record_id in callable_ids
            if cds[record_id] != nucleotides[record_id]
        ]
        if sequence_mismatches:
            preview = ", ".join(sorted(sequence_mismatches)[:5])
            raise ValueError(
                f"Primary-CDS and historical nucleotide sequences differ for "
                f"{len(sequence_mismatches)} callable IDs: {preview}"
            )
    except (OSError, UnicodeError, UnicodeEncodeError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}")

    excluded_ids = coordinate_set - callable_ids
    out_parent = args.output_dir.parent
    out_parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.tmp.", dir=out_parent))
    try:
        write_tsv(
            tmp_dir / "callable_reference_ids.tsv",
            ["reference_gene", "nucleotide_length", "sequence_sha256"],
            (
                {
                    "reference_gene": record_id,
                    "nucleotide_length": cds[record_id].length,
                    "sequence_sha256": cds[record_id].sequence_sha256,
                }
                for record_id in sorted(callable_ids)
            ),
        )
        write_tsv(
            tmp_dir / "excluded_reference_coordinate_ids.tsv",
            ["reference_gene", "exclusion_reasons"],
            (
                {
                    "reference_gene": record_id,
                    "exclusion_reasons": ";".join(
                        reason for reason, present in [
                            ("missing_historical_protein", record_id in proteins),
                            ("missing_primary_cds", record_id in cds),
                            ("missing_historical_nucleotide", record_id in nucleotides),
                        ] if not present
                    ),
                }
                for record_id in sorted(excluded_ids)
            ),
        )
        record_counts = {
            "reference_coords": len(coordinate_ids),
            "reference_proteins": len(proteins),
            "reference_nucleotides": len(nucleotides),
            "reference_cds": len(cds),
        }
        write_tsv(
            tmp_dir / "input_checksums.tsv",
            ["role", "path", "sha256", "records", "unique_ids"],
            (
                {
                    "role": role,
                    "path": str(path),
                    "sha256": sha256_file(path),
                    "records": record_counts[role],
                    "unique_ids": record_counts[role],
                }
                for role, path in inputs.items()
            ),
        )
        metadata_rows = [
            ("timestamp_utc", datetime.now(timezone.utc).isoformat()),
            ("command", " ".join(shlex.quote(value) for value in sys.argv)),
            ("python_version", platform.python_version()),
            ("coordinate_reference_ids", len(coordinate_set)),
            ("callable_reference_ids", len(callable_ids)),
            ("excluded_coordinate_ids", len(excluded_ids)),
            ("historical_nucleotide_ids", len(nucleotides)),
            ("historical_nucleotide_ids_outside_coordinates", len(set(nucleotides) - coordinate_set)),
            ("sequence_mismatches_between_callable_cds_and_historical_nucleotides", 0),
            ("callable_rule", "coordinate ID present in historical protein, historical nucleotide, and primary-CDS FASTAs"),
            ("thread_contract", "single Python process; no NLR-Annotator or Java process launched"),
        ]
        write_tsv(
            tmp_dir / "run_metadata.tsv",
            ["key", "value"],
            ({"key": key, "value": value} for key, value in metadata_rows),
        )
        (tmp_dir / "README.txt").write_text(
            "This audit defines the frozen reference sequence universe eligible for the "
            "reference-centric NLR loss denominator. It does not contain an NLR count; "
            "NLR-Annotator must still be run on the audited primary-CDS FASTA.\n",
            encoding="utf-8",
        )
        os.replace(tmp_dir, args.output_dir)
    except Exception as exc:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise SystemExit(f"ERROR: {exc}")

    print(
        f"Audited {len(coordinate_set)} coordinate IDs; {len(callable_ids)} are callable; "
        f"results: {args.output_dir}"
    )


if __name__ == "__main__":
    main()
