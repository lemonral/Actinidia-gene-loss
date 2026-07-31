#!/usr/bin/env python3
"""Materialize one checksum-bound plain-FASTA bundle for NLR-Annotator.

The source manifest may mix gzip-compressed and plain FASTA files.  Every
source is streamed into a regular plain FASTA under a new atomic output
directory.  The generated ``nlr_annotator_inputs.tsv`` is accepted directly
by ``run_nlr_annotator_batch.py`` with the bundle directory as ``--input-root``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterable


SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
SAFE_BASENAME = re.compile(r"^[A-Za-z0-9_.-]+\.fa$")
SOURCE_FIELDS = (
    "sample_id",
    "species",
    "ploidy",
    "analysis_role",
    "input_scope",
    "source_fasta",
    "output_basename",
    "expected_fasta_records",
)
RUNNER_FIELDS = (
    "sample_id",
    "species",
    "ploidy",
    "analysis_role",
    "input_scope",
    "relative_fasta",
    "expected_fasta_records",
)


@dataclass(frozen=True)
class SourceRow:
    sample_id: str
    species: str
    ploidy: str
    analysis_role: str
    input_scope: str
    source_relative: str
    source: Path
    output_basename: str
    expected_records: int | None


@dataclass(frozen=True)
class FastaAudit:
    records: int
    total_bases: int
    sha256: str
    bytes: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-targets", type=int, required=True)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_sources(manifest: Path, input_root: Path) -> list[SourceRow]:
    with manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SOURCE_FIELDS:
            raise ValueError(
                f"Unexpected source-manifest columns: {reader.fieldnames}; "
                f"expected {list(SOURCE_FIELDS)}"
            )
        rows: list[SourceRow] = []
        seen_ids: set[str] = set()
        seen_outputs: set[str] = set()
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise ValueError(f"Manifest line {line_number} has extra fields")
            row = {key: (value or "").strip() for key, value in raw.items()}
            sample = row["sample_id"]
            if SAFE_ID.fullmatch(sample) is None or sample in seen_ids:
                raise ValueError(f"Unsafe or duplicate sample_id at line {line_number}: {sample!r}")
            seen_ids.add(sample)
            output = row["output_basename"]
            if SAFE_BASENAME.fullmatch(output) is None or output in seen_outputs:
                raise ValueError(
                    f"Unsafe or duplicate output_basename at line {line_number}: {output!r}"
                )
            seen_outputs.add(output)
            relative = Path(row["source_fasta"])
            if not row["source_fasta"] or relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"source_fasta must be safe and relative at line {line_number}")
            role = row["analysis_role"]
            scope = row["input_scope"]
            if role == "reference_callable":
                if scope != "reference_transcript_cds":
                    raise ValueError("The reference row must use reference_transcript_cds")
            elif role == "target_repertoire":
                if scope != "whole_genome":
                    raise ValueError("Target rows must use whole_genome")
            else:
                raise ValueError(f"Unsupported analysis_role at line {line_number}: {role!r}")
            expected: int | None = None
            if row["expected_fasta_records"]:
                try:
                    expected = int(row["expected_fasta_records"])
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid expected_fasta_records at line {line_number}"
                    ) from exc
                if expected < 1:
                    raise ValueError("expected_fasta_records must be positive")
            rows.append(
                SourceRow(
                    sample_id=sample,
                    species=row["species"],
                    ploidy=row["ploidy"],
                    analysis_role=role,
                    input_scope=scope,
                    source_relative=row["source_fasta"],
                    source=input_root / relative,
                    output_basename=output,
                    expected_records=expected,
                )
            )
    if not rows:
        raise ValueError("Source manifest is empty")
    if sum(row.analysis_role == "reference_callable" for row in rows) != 1:
        raise ValueError("Source manifest must contain exactly one reference_callable row")
    return rows


def is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def audit_plain_stream(source: BinaryIO, destination: BinaryIO) -> FastaAudit:
    digest = hashlib.sha256()
    records = 0
    total_bases = 0
    total_bytes = 0
    seen: set[bytes] = set()
    current_id: bytes | None = None
    current_bases = 0
    for line_number, line in enumerate(source, start=1):
        destination.write(line)
        digest.update(line)
        total_bytes += len(line)
        if line.startswith(b">"):
            if current_id is not None and current_bases == 0:
                raise ValueError(f"Empty FASTA record {current_id!r}")
            header = line[1:].strip()
            if not header:
                raise ValueError(f"Empty FASTA header at line {line_number}")
            current_id = header.split()[0]
            if current_id in seen:
                raise ValueError(f"Duplicate FASTA ID {current_id!r}")
            seen.add(current_id)
            records += 1
            current_bases = 0
        else:
            if current_id is None:
                if line.strip():
                    raise ValueError(f"Sequence data precedes first header at line {line_number}")
                continue
            bases = len(b"".join(line.split()))
            current_bases += bases
            total_bases += bases
    if current_id is not None and current_bases == 0:
        raise ValueError(f"Empty FASTA record {current_id!r}")
    if records == 0:
        raise ValueError("No FASTA records found")
    return FastaAudit(records, total_bases, digest.hexdigest(), total_bytes)


def materialize(row: SourceRow, destination: Path) -> tuple[FastaAudit, str, int]:
    if not row.source.is_file():
        raise ValueError(f"Source FASTA is missing for {row.sample_id}: {row.source_relative}")
    source_hash = sha256_file(row.source)
    source_bytes = row.source.stat().st_size
    opener = gzip.open if is_gzip(row.source) else open
    with opener(row.source, "rb") as source, destination.open("wb") as output:
        audit = audit_plain_stream(source, output)
    if row.expected_records is not None and audit.records != row.expected_records:
        raise ValueError(
            f"{row.sample_id} has {audit.records} FASTA records; expected {row.expected_records}"
        )
    return audit, source_hash, source_bytes


def write_tsv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, object]]) -> None:
    field_list = list(fields)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=field_list, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in field_list})


def main() -> int:
    args = parse_args()
    if args.expected_targets < 1:
        print("ERROR: --expected-targets must be positive", file=sys.stderr)
        return 2
    if args.output_dir.exists():
        print(f"ERROR: output directory already exists: {args.output_dir}", file=sys.stderr)
        return 2
    try:
        rows = read_sources(args.manifest, args.input_root)
        targets = sum(row.analysis_role == "target_repertoire" for row in rows)
        if targets != args.expected_targets:
            raise ValueError(f"Observed {targets} target rows; expected {args.expected_targets}")
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(
            tempfile.mkdtemp(prefix=f".{args.output_dir.name}.tmp.", dir=args.output_dir.parent)
        )
        try:
            outputs: list[dict[str, object]] = []
            sources: list[dict[str, object]] = []
            runner_rows: list[dict[str, object]] = []
            for row in rows:
                destination = staging / row.output_basename
                audit, source_hash, source_bytes = materialize(row, destination)
                sources.append(
                    {
                        "sample_id": row.sample_id,
                        "relative_path": row.source_relative,
                        "bytes": source_bytes,
                        "sha256": source_hash,
                        "compression": "gzip" if is_gzip(row.source) else "plain",
                    }
                )
                outputs.append(
                    {
                        "sample_id": row.sample_id,
                        "basename": row.output_basename,
                        "bytes": audit.bytes,
                        "sha256": audit.sha256,
                        "fasta_records": audit.records,
                        "total_bases": audit.total_bases,
                    }
                )
                runner_rows.append(
                    {
                        "sample_id": row.sample_id,
                        "species": row.species,
                        "ploidy": row.ploidy,
                        "analysis_role": row.analysis_role,
                        "input_scope": row.input_scope,
                        "relative_fasta": row.output_basename,
                        "expected_fasta_records": audit.records,
                    }
                )
            runner_manifest = staging / "nlr_annotator_inputs.tsv"
            write_tsv(runner_manifest, RUNNER_FIELDS, runner_rows)
            outputs.append(
                {
                    "sample_id": "bundle",
                    "basename": runner_manifest.name,
                    "bytes": runner_manifest.stat().st_size,
                    "sha256": sha256_file(runner_manifest),
                }
            )
            manifest = {
                "schema_version": 1,
                "status": "PASS",
                "workflow": "plain_fasta_nlr_input_bundle",
                "source_manifest": {
                    "basename": args.manifest.name,
                    "bytes": args.manifest.stat().st_size,
                    "sha256": sha256_file(args.manifest),
                },
                "reference_count": 1,
                "target_count": targets,
                "sources": sources,
                "outputs": outputs,
            }
            (staging / "run_manifest.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            write_tsv(
                staging / "checksums.tsv",
                ("file", "bytes", "sha256"),
                (
                    {
                        "file": item["basename"],
                        "bytes": item["bytes"],
                        "sha256": item["sha256"],
                    }
                    for item in outputs
                ),
            )
            os.replace(staging, args.output_dir)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
    except (OSError, EOFError, ValueError, gzip.BadGzipFile) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}\ttargets={args.expected_targets}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
