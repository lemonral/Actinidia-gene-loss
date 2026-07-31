#!/usr/bin/env python3
"""Audit analysis FASTA scope against declared full or candidate FASTA files.

The audit is intentionally single-process and fail-closed.  FASTA sequence
bodies are streamed, normalized to uppercase after whitespace removal, and
reduced to record ID, length, and SHA-256.  They are never retained in memory.
Repeated deposit candidates are indexed once and served from an in-process
cache, which is important for combined legacy polyploid FASTA files shared by
several assembly-unit rows.

Matching proceeds in two ordered phases:

1. exact, case-sensitive FASTA ID plus exact normalized sequence;
2. a unique length-plus-sequence-SHA-256 match for renamed records.

A repeated FASTA ID or an ambiguous sequence-hash-only match is fatal.  A
same-ID/different-sequence conflict is fatal for all evidence classes except a
name-only unverified candidate, where it is recorded explicitly as evidence
that the candidate cannot safely fill the analysis-input scope.  An official
same-accession full deposit
also requires an NCBI assembly report whose sequence IDs and lengths reconcile
exactly to the deposit FASTA.  Only that evidence class may translate an
assembly-report ``unplaced-scaffold`` role into an official unplaced result.
Combined-haplotype and name-only legacy candidates remain candidate evidence:
their extra records are always labelled ``candidate_only_scope_unknown``.

Outputs are written to a sibling staging directory and installed with one
atomic rename only after every manifest row and output has passed.  The output
directory must not already exist.  The program starts no worker threads.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, TextIO


SCRIPT_VERSION = "1.2.0"
READ_CHUNK_BYTES = 8 * 1024 * 1024
GZIP_MAGIC = b"\x1f\x8b"

PAIR_CLASSES = (
    "official_same_accession_full_deposit",
    "same_genome_matched_source",
    "combined_haplotype_legacy_candidate",
    "unverified_name_candidate",
    "unavailable",
)

MANIFEST_COLUMNS = (
    "assembly_unit_id",
    "biological_species",
    "pair_class",
    "analysis_accession",
    "deposit_accession",
    "analysis_fasta",
    "deposit_fasta",
    "assembly_report",
    "analysis_scope",
    "deposit_scope",
    "evidence_note",
)

SUMMARY_COLUMNS = (
    "assembly_unit_id",
    "biological_species",
    "pair_class",
    "scope_gap_authority",
    "analysis_accession",
    "deposit_accession",
    "analysis_scope",
    "deposit_scope",
    "analysis_fasta",
    "deposit_fasta",
    "assembly_report",
    "analysis_record_count",
    "analysis_total_bp",
    "deposit_record_count",
    "deposit_total_bp",
    "exact_id_sequence_match_count",
    "exact_id_sequence_match_bp",
    "sequence_hash_only_match_count",
    "sequence_hash_only_match_bp",
    "same_id_sequence_mismatch_count",
    "same_id_sequence_mismatch_analysis_bp",
    "same_id_sequence_mismatch_deposit_bp",
    "matched_analysis_record_count",
    "matched_analysis_bp",
    "analysis_bp_reconciled_percent",
    "analysis_only_record_count",
    "analysis_only_bp",
    "deposit_only_record_count",
    "deposit_only_bp",
    "deposit_only_bp_percent_of_deposit",
    "official_unplaced_record_count",
    "official_unplaced_bp",
    "official_unplaced_bp_percent_of_deposit",
    "report_assembled_molecule_record_count",
    "report_assembled_molecule_bp",
    "report_unlocalized_scaffold_record_count",
    "report_unlocalized_scaffold_bp",
    "report_unplaced_scaffold_record_count",
    "report_unplaced_scaffold_bp",
    "report_other_role_record_count",
    "report_other_role_bp",
    "record_reconciliation_status",
    "deposit_only_interpretation",
    "evidence_note",
)

MATCH_COLUMNS = (
    "assembly_unit_id",
    "pair_class",
    "analysis_record_id",
    "deposit_record_id",
    "match_method",
    "sequence_bp",
    "sequence_sha256",
)

ANALYSIS_ONLY_COLUMNS = (
    "assembly_unit_id",
    "pair_class",
    "analysis_record_id",
    "sequence_bp",
    "sequence_sha256",
    "reason",
)

DEPOSIT_ONLY_COLUMNS = (
    "assembly_unit_id",
    "pair_class",
    "deposit_record_id",
    "sequence_bp",
    "sequence_sha256",
    "assembly_report_role",
    "scope_interpretation",
)

CONFLICT_COLUMNS = (
    "assembly_unit_id",
    "pair_class",
    "record_id",
    "analysis_sequence_bp",
    "analysis_sequence_sha256",
    "deposit_sequence_bp",
    "deposit_sequence_sha256",
    "interpretation",
)

INPUT_COLUMNS = (
    "input_type",
    "consumer_assembly_units",
    "resolved_path",
    "physical_size_bytes",
    "physical_sha256",
    "record_count",
    "logical_sequence_bp",
    "logical_index_sha256",
)


class ScopeAuditError(RuntimeError):
    """Raised when a trustworthy scope audit cannot be published."""


@dataclass(frozen=True)
class ManifestRow:
    """One validated manifest declaration."""

    assembly_unit_id: str
    biological_species: str
    pair_class: str
    analysis_accession: str
    deposit_accession: str
    analysis_fasta: Path
    deposit_fasta: Path | None
    assembly_report: Path | None
    analysis_scope: str
    deposit_scope: str
    evidence_note: str


@dataclass(frozen=True)
class FastaRecord:
    """One indexed FASTA record without its sequence body."""

    record_id: str
    length: int
    sequence_sha256: str


@dataclass(frozen=True)
class FastaIndex:
    """A deterministic, bounded-memory logical FASTA index."""

    path: Path
    physical_size: int
    physical_sha256: str
    logical_index_sha256: str
    total_bp: int
    records: tuple[FastaRecord, ...]


@dataclass(frozen=True)
class AssemblyReportRecord:
    """One NCBI assembly-report sequence row."""

    sequence_name: str
    role: str
    length: int
    genbank_accession: str
    refseq_accession: str

    @property
    def identifiers(self) -> tuple[str, ...]:
        """Return official identifiers that may occur as a FASTA first token."""
        values = (
            self.sequence_name,
            self.genbank_accession,
            self.refseq_accession,
        )
        return tuple(value for value in values if value and value.lower() != "na")


@dataclass(frozen=True)
class AssemblyReportIndex:
    """A validated NCBI assembly-report index."""

    path: Path
    physical_size: int
    physical_sha256: str
    logical_index_sha256: str
    records: tuple[AssemblyReportRecord, ...]


@dataclass(frozen=True)
class RecordMatch:
    """One analysis-to-deposit record reconciliation."""

    analysis: FastaRecord
    deposit: FastaRecord
    method: str


@dataclass(frozen=True)
class RecordConflict:
    """One same-ID/different-sequence conflict in a name-only candidate."""

    analysis: FastaRecord
    deposit: FastaRecord


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="Exact-schema paired assembly-scope TSV manifest.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="New output directory; it must not already exist.",
    )
    parser.add_argument(
        "--expected-assembly-unit-count",
        type=int,
        default=None,
        help=(
            "Optional exact number of unique assembly-unit rows. By default the "
            "cohort size is derived from the non-empty manifest."
        ),
    )
    return parser.parse_args()


def absolute_lexical(path: Path, base: Path | None = None) -> Path:
    """Return an absolute path without dereferencing the final component."""
    expanded = path.expanduser()
    if expanded.is_absolute():
        return expanded
    anchor = Path.cwd() if base is None else base
    return anchor / expanded


def require_text(row: dict[str, str | None], column: str, row_number: int) -> str:
    value = row.get(column)
    if value is None or not value.strip():
        raise ScopeAuditError(
            f"Manifest row {row_number} has an empty required value: {column}"
        )
    return value.strip()


def optional_text(row: dict[str, str | None], column: str) -> str:
    value = row.get(column)
    return "" if value is None else value.strip()


def resolve_manifest_path(value: str, manifest_parent: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest_parent / path
    return path


def require_regular_file(path: Path, description: str) -> None:
    if not path.exists():
        raise ScopeAuditError(f"{description} does not exist: {path}")
    if not path.is_file():
        raise ScopeAuditError(f"{description} is not a regular file: {path}")


def load_manifest(
    path: Path, expected_count: int | None = None
) -> list[ManifestRow]:
    """Load the complete assembly-unit pairing contract.

    The manifest itself defines the cohort unless ``expected_count`` is supplied
    as a fail-closed external assertion.
    """
    if expected_count is not None and expected_count <= 0:
        raise ScopeAuditError("--expected-assembly-unit-count must be positive")
    manifest = absolute_lexical(path)
    require_regular_file(manifest, "Manifest")
    parent = manifest.parent
    try:
        handle = manifest.open("r", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise ScopeAuditError(f"Cannot open manifest {manifest}: {error}") from error

    rows: list[ManifestRow] = []
    seen: set[str] = set()
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(MANIFEST_COLUMNS):
            observed = "\t".join(reader.fieldnames or [])
            expected = "\t".join(MANIFEST_COLUMNS)
            raise ScopeAuditError(
                "Manifest header must match exactly. "
                f"Expected {expected!r}; observed {observed!r}"
            )

        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise ScopeAuditError(
                    f"Manifest row {row_number} contains more fields than the header"
                )
            assembly_unit_id = require_text(row, "assembly_unit_id", row_number)
            if assembly_unit_id in seen:
                raise ScopeAuditError(
                    f"Duplicate assembly_unit_id in manifest: {assembly_unit_id}"
                )
            if any(character.isspace() for character in assembly_unit_id):
                raise ScopeAuditError(
                    "Manifest assembly_unit_id contains whitespace: "
                    f"{assembly_unit_id!r}"
                )
            seen.add(assembly_unit_id)

            pair_class = require_text(row, "pair_class", row_number)
            if pair_class not in PAIR_CLASSES:
                raise ScopeAuditError(
                    f"Manifest row {row_number} has unsupported pair_class {pair_class!r}; "
                    "expected one of " + ", ".join(PAIR_CLASSES)
                )

            analysis_value = require_text(row, "analysis_fasta", row_number)
            analysis_path = resolve_manifest_path(analysis_value, parent)
            require_regular_file(
                analysis_path, f"Analysis FASTA for {assembly_unit_id}"
            )

            deposit_value = optional_text(row, "deposit_fasta")
            report_value = optional_text(row, "assembly_report")
            deposit_accession = optional_text(row, "deposit_accession")
            deposit_path = (
                resolve_manifest_path(deposit_value, parent) if deposit_value else None
            )
            report_path = (
                resolve_manifest_path(report_value, parent) if report_value else None
            )

            if pair_class == "unavailable":
                if deposit_path is not None or report_path is not None or deposit_accession:
                    raise ScopeAuditError(
                        f"Unavailable row {assembly_unit_id} must leave deposit_accession, "
                        "deposit_fasta, and assembly_report empty"
                    )
                if require_text(row, "deposit_scope", row_number) != "unavailable":
                    raise ScopeAuditError(
                        f"Unavailable row {assembly_unit_id} must use "
                        "deposit_scope=unavailable"
                    )
            else:
                if deposit_path is None:
                    raise ScopeAuditError(
                        f"Manifest row {assembly_unit_id} requires a deposit_fasta"
                    )
                require_regular_file(
                    deposit_path, f"Deposit FASTA for {assembly_unit_id}"
                )
                if not deposit_accession:
                    raise ScopeAuditError(
                        f"Manifest row {assembly_unit_id} requires deposit_accession"
                    )
                if report_path is not None:
                    require_regular_file(
                        report_path, f"Assembly report for {assembly_unit_id}"
                    )

            analysis_accession = require_text(row, "analysis_accession", row_number)
            if pair_class in {
                "official_same_accession_full_deposit",
                "same_genome_matched_source",
            } and analysis_accession != deposit_accession:
                raise ScopeAuditError(
                    f"{pair_class} row {assembly_unit_id} must declare identical analysis "
                    "and deposit accessions"
                )

            if pair_class == "official_same_accession_full_deposit":
                if report_path is None:
                    raise ScopeAuditError(
                        f"Official full-deposit row {assembly_unit_id} requires "
                        "assembly_report"
                    )
                assert deposit_path is not None
                if analysis_path.resolve() == deposit_path.resolve():
                    raise ScopeAuditError(
                        f"Official full-deposit row {assembly_unit_id} points analysis and "
                        "deposit to the same file"
                    )

            rows.append(
                ManifestRow(
                    assembly_unit_id=assembly_unit_id,
                    biological_species=require_text(
                        row, "biological_species", row_number
                    ),
                    pair_class=pair_class,
                    analysis_accession=analysis_accession,
                    deposit_accession=deposit_accession,
                    analysis_fasta=analysis_path,
                    deposit_fasta=deposit_path,
                    assembly_report=report_path,
                    analysis_scope=require_text(row, "analysis_scope", row_number),
                    deposit_scope=require_text(row, "deposit_scope", row_number),
                    evidence_note=require_text(row, "evidence_note", row_number),
                )
            )

    if not rows:
        raise ScopeAuditError("Manifest contains no assembly-unit rows")
    if expected_count is not None and len(rows) != expected_count:
        raise ScopeAuditError(
            f"Manifest must contain exactly {expected_count} assembly-unit rows; "
            f"observed {len(rows)}"
        )
    return rows


def sha256_file(path: Path) -> tuple[int, str]:
    """Return exact physical bytes and SHA-256 from one bounded-memory pass."""
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
        raise ScopeAuditError(f"Cannot checksum {path}: {error}") from error
    return size, digest.hexdigest()


def detect_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == GZIP_MAGIC
    except OSError as error:
        raise ScopeAuditError(f"Cannot inspect FASTA compression {path}: {error}") from error


@contextmanager
def open_fasta_text(path: Path) -> Iterator[TextIO]:
    """Open plain/gzip FASTA according to magic bytes, not its filename."""
    try:
        if detect_gzip(path):
            with gzip.open(
                path, "rt", encoding="utf-8", errors="strict", newline=""
            ) as handle:
                yield handle
        else:
            with path.open(
                "rt", encoding="utf-8", errors="strict", newline=""
            ) as handle:
                yield handle
    except (OSError, EOFError, UnicodeError, gzip.BadGzipFile) as error:
        raise ScopeAuditError(f"Cannot read FASTA {path}: {error}") from error


def index_fasta(path: Path) -> FastaIndex:
    """Stream and index one FASTA without retaining sequence bodies."""
    resolved = path.resolve(strict=True)
    physical_size, physical_sha256 = sha256_file(resolved)
    records: list[FastaRecord] = []
    seen: set[str] = set()
    current_id: str | None = None
    current_length = 0
    current_digest = hashlib.sha256()

    def finish_record() -> None:
        nonlocal current_id, current_length, current_digest
        if current_id is None:
            return
        if current_length == 0:
            raise ScopeAuditError(
                f"FASTA record {current_id!r} has an empty sequence in {resolved}"
            )
        records.append(
            FastaRecord(
                record_id=current_id,
                length=current_length,
                sequence_sha256=current_digest.hexdigest(),
            )
        )

    with open_fasta_text(resolved) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                finish_record()
                header = stripped[1:].strip()
                if not header:
                    raise ScopeAuditError(
                        f"Blank FASTA header in {resolved} at line {line_number}"
                    )
                record_id = header.split()[0]
                if record_id in seen:
                    raise ScopeAuditError(
                        f"Duplicate FASTA identifier {record_id!r} in {resolved}"
                    )
                seen.add(record_id)
                current_id = record_id
                current_length = 0
                current_digest = hashlib.sha256()
                continue
            if current_id is None:
                raise ScopeAuditError(
                    f"Sequence text precedes the first FASTA header in {resolved} "
                    f"at line {line_number}"
                )
            normalized = "".join(stripped.split()).upper()
            if not normalized:
                continue
            try:
                encoded = normalized.encode("ascii")
            except UnicodeEncodeError as error:
                raise ScopeAuditError(
                    f"Non-ASCII sequence text in {resolved} at line {line_number}"
                ) from error
            current_length += len(encoded)
            current_digest.update(encoded)

    if current_id is None:
        raise ScopeAuditError(f"No FASTA records found: {resolved}")
    finish_record()

    logical = hashlib.sha256()
    for record in records:
        logical.update(
            (
                f"{record.record_id}\t{record.length}\t"
                f"{record.sequence_sha256}\n"
            ).encode("utf-8")
        )
    return FastaIndex(
        path=resolved,
        physical_size=physical_size,
        physical_sha256=physical_sha256,
        logical_index_sha256=logical.hexdigest(),
        total_bp=sum(record.length for record in records),
        records=tuple(records),
    )


class FastaIndexCache:
    """Cache unique resolved FASTA paths and record their consumers."""

    def __init__(self) -> None:
        self._indexes: dict[Path, FastaIndex] = {}
        self.consumers: dict[Path, set[str]] = defaultdict(set)
        self.input_types: dict[Path, set[str]] = defaultdict(set)
        self.hits = 0
        self.misses = 0

    def get(
        self, path: Path, assembly_unit_id: str, input_type: str
    ) -> FastaIndex:
        resolved = path.resolve(strict=True)
        self.consumers[resolved].add(assembly_unit_id)
        self.input_types[resolved].add(input_type)
        if resolved in self._indexes:
            self.hits += 1
            return self._indexes[resolved]
        self.misses += 1
        index = index_fasta(resolved)
        self._indexes[resolved] = index
        return index

    @property
    def indexes(self) -> dict[Path, FastaIndex]:
        return dict(self._indexes)


def parse_assembly_report(path: Path) -> AssemblyReportIndex:
    """Parse one NCBI assembly report and require unique sequence rows."""
    resolved = path.resolve(strict=True)
    physical_size, physical_sha256 = sha256_file(resolved)
    header: list[str] | None = None
    records: list[AssemblyReportRecord] = []
    seen: set[str] = set()
    try:
        handle = resolved.open("r", encoding="utf-8-sig", errors="strict")
    except OSError as error:
        raise ScopeAuditError(f"Cannot open assembly report {resolved}: {error}") from error
    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if header is None:
                if line.startswith("# Sequence-Name\t"):
                    header = line[2:].split("\t")
                continue
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != len(header):
                raise ScopeAuditError(
                    f"Assembly report {resolved} line {line_number} has "
                    f"{len(fields)} fields; expected {len(header)}"
                )
            row = dict(zip(header, fields))
            sequence_name = row.get("Sequence-Name", "").strip()
            role = row.get("Sequence-Role", "").strip()
            length_text = row.get("Sequence-Length", "").strip()
            if not sequence_name or not role or not length_text:
                raise ScopeAuditError(
                    f"Assembly report {resolved} line {line_number} has an empty "
                    "Sequence-Name, Sequence-Role, or Sequence-Length"
                )
            if sequence_name in seen:
                raise ScopeAuditError(
                    f"Duplicate Sequence-Name {sequence_name!r} in {resolved}"
                )
            try:
                length = int(length_text)
            except ValueError as error:
                raise ScopeAuditError(
                    f"Assembly report {resolved} line {line_number} has invalid "
                    f"Sequence-Length {length_text!r}"
                ) from error
            if length <= 0:
                raise ScopeAuditError(
                    f"Assembly report {resolved} line {line_number} has non-positive "
                    f"Sequence-Length {length}"
                )
            seen.add(sequence_name)
            records.append(
                AssemblyReportRecord(
                    sequence_name=sequence_name,
                    role=role,
                    length=length,
                    genbank_accession=row.get("GenBank-Accn", "").strip(),
                    refseq_accession=row.get("RefSeq-Accn", "").strip(),
                )
            )
    if header is None:
        raise ScopeAuditError(
            f"Assembly report lacks the required Sequence-Name header: {resolved}"
        )
    required = {"Sequence-Name", "Sequence-Role", "Sequence-Length"}
    missing = sorted(required.difference(header))
    if missing:
        raise ScopeAuditError(
            f"Assembly report {resolved} is missing columns: " + ", ".join(missing)
        )
    if not records:
        raise ScopeAuditError(f"Assembly report contains no sequence rows: {resolved}")

    logical = hashlib.sha256()
    for record in records:
        logical.update(
            (
                f"{record.sequence_name}\t{record.role}\t{record.length}\t"
                f"{record.genbank_accession}\t{record.refseq_accession}\n"
            ).encode("utf-8")
        )
    return AssemblyReportIndex(
        path=resolved,
        physical_size=physical_size,
        physical_sha256=physical_sha256,
        logical_index_sha256=logical.hexdigest(),
        records=tuple(records),
    )


def verify_official_report(
    deposit: FastaIndex, report: AssemblyReportIndex, assembly_unit_id: str
) -> dict[str, AssemblyReportRecord]:
    """Require one-to-one report/deposit identity through official identifiers."""
    aliases: dict[str, AssemblyReportRecord] = {}
    for report_record in report.records:
        for identifier in report_record.identifiers:
            previous = aliases.get(identifier)
            if previous is not None and previous != report_record:
                raise ScopeAuditError(
                    f"Assembly report identifier {identifier!r} maps to multiple "
                    f"sequence rows for {assembly_unit_id}"
                )
            aliases[identifier] = report_record

    report_by_deposit_id: dict[str, AssemblyReportRecord] = {}
    matched_sequence_names: set[str] = set()
    deposit_only: list[str] = []
    length_mismatches: list[str] = []
    for deposit_record in deposit.records:
        report_record = aliases.get(deposit_record.record_id)
        if report_record is None:
            deposit_only.append(deposit_record.record_id)
            continue
        if report_record.sequence_name in matched_sequence_names:
            raise ScopeAuditError(
                f"Multiple deposit FASTA records map to assembly-report row "
                f"{report_record.sequence_name!r} for {assembly_unit_id}"
            )
        matched_sequence_names.add(report_record.sequence_name)
        report_by_deposit_id[deposit_record.record_id] = report_record
        if deposit_record.length != report_record.length:
            length_mismatches.append(deposit_record.record_id)

    report_only = sorted(
        record.sequence_name
        for record in report.records
        if record.sequence_name not in matched_sequence_names
    )
    if deposit_only or report_only:
        raise ScopeAuditError(
            f"Official assembly-report ID set does not match deposit FASTA for "
            f"{assembly_unit_id}: deposit_only={deposit_only[:5]}, "
            f"report_only={report_only[:5]}"
        )
    if length_mismatches:
        raise ScopeAuditError(
            f"Official assembly-report lengths do not match deposit FASTA for "
            f"{assembly_unit_id}: {length_mismatches[:5]}"
        )
    return report_by_deposit_id


def reconcile_records(
    analysis: FastaIndex,
    deposit: FastaIndex,
    assembly_unit_id: str,
    *,
    allow_same_id_sequence_mismatch: bool = False,
) -> tuple[
    list[RecordMatch],
    list[FastaRecord],
    list[FastaRecord],
    list[RecordConflict],
]:
    """Match exact IDs first, then unique sequence hashes."""
    deposit_by_id = {record.record_id: record for record in deposit.records}
    matched_deposit_ids: set[str] = set()
    matches: list[RecordMatch] = []
    remaining_analysis: list[FastaRecord] = []
    conflicts: list[RecordConflict] = []

    for analysis_record in analysis.records:
        deposit_record = deposit_by_id.get(analysis_record.record_id)
        if deposit_record is None:
            remaining_analysis.append(analysis_record)
            continue
        if (
            analysis_record.length != deposit_record.length
            or analysis_record.sequence_sha256 != deposit_record.sequence_sha256
        ):
            if allow_same_id_sequence_mismatch:
                conflicts.append(
                    RecordConflict(analysis=analysis_record, deposit=deposit_record)
                )
                # The name conflict invalidates an exact-ID match but does not
                # preclude a unique sequence-hash match to another record.  This
                # is important for detecting a systematic chromosome renaming
                # problem without trusting the candidate's names.
                remaining_analysis.append(analysis_record)
                continue
            raise ScopeAuditError(
                f"Exact FASTA ID {analysis_record.record_id!r} has a different "
                f"sequence in analysis and deposit files for {assembly_unit_id}"
            )
        matched_deposit_ids.add(deposit_record.record_id)
        matches.append(
            RecordMatch(
                analysis=analysis_record,
                deposit=deposit_record,
                method="exact_id_and_sequence",
            )
        )

    deposit_by_hash: dict[tuple[int, str], list[FastaRecord]] = defaultdict(list)
    for deposit_record in deposit.records:
        if deposit_record.record_id not in matched_deposit_ids:
            deposit_by_hash[
                (deposit_record.length, deposit_record.sequence_sha256)
            ].append(deposit_record)

    analysis_only: list[FastaRecord] = []
    for analysis_record in remaining_analysis:
        key = (analysis_record.length, analysis_record.sequence_sha256)
        candidates = [
            record
            for record in deposit_by_hash.get(key, [])
            if record.record_id not in matched_deposit_ids
        ]
        if len(candidates) > 1:
            raise ScopeAuditError(
                f"Ambiguous sequence-hash-only match for analysis record "
                f"{analysis_record.record_id!r} in {assembly_unit_id}: "
                + ", ".join(record.record_id for record in candidates[:8])
            )
        if not candidates:
            analysis_only.append(analysis_record)
            continue
        deposit_record = candidates[0]
        matched_deposit_ids.add(deposit_record.record_id)
        matches.append(
            RecordMatch(
                analysis=analysis_record,
                deposit=deposit_record,
                method="unique_sequence_hash",
            )
        )

    deposit_only = [
        record
        for record in deposit.records
        if record.record_id not in matched_deposit_ids
    ]
    return matches, analysis_only, deposit_only, conflicts


def percent(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return ""
    return f"{100.0 * numerator / denominator:.6f}"


def scope_gap_authority(pair_class: str) -> str:
    return {
        "official_same_accession_full_deposit": "official_assembly_report",
        "same_genome_matched_source": "same_genome_reconciliation",
        "combined_haplotype_legacy_candidate": "candidate_only_no_scope_role_inference",
        "unverified_name_candidate": "candidate_only_no_scope_role_inference",
        "unavailable": "unavailable",
    }[pair_class]


def deposit_only_interpretation(pair_class: str) -> str:
    if pair_class == "official_same_accession_full_deposit":
        return "official_assembly_report_roles"
    if pair_class == "same_genome_matched_source":
        return "same_genome_requires_no_unmatched_records"
    if pair_class in {
        "combined_haplotype_legacy_candidate",
        "unverified_name_candidate",
    }:
        return "candidate_only_scope_unknown"
    return "deposit_unavailable"


def report_role_totals(
    report_by_id: dict[str, AssemblyReportRecord] | None,
) -> dict[str, int]:
    totals = {
        "assembled_molecule_records": 0,
        "assembled_molecule_bp": 0,
        "unlocalized_scaffold_records": 0,
        "unlocalized_scaffold_bp": 0,
        "unplaced_scaffold_records": 0,
        "unplaced_scaffold_bp": 0,
        "other_records": 0,
        "other_bp": 0,
    }
    if report_by_id is None:
        return totals
    for record in report_by_id.values():
        if record.role == "assembled-molecule":
            prefix = "assembled_molecule"
        elif record.role == "unlocalized-scaffold":
            prefix = "unlocalized_scaffold"
        elif record.role == "unplaced-scaffold":
            prefix = "unplaced_scaffold"
        else:
            prefix = "other"
        totals[f"{prefix}_records"] += 1
        totals[f"{prefix}_bp"] += record.length
    return totals


def classify_official_deposit_only(role: str) -> str:
    if role == "unplaced-scaffold":
        return "official_unplaced_scaffold"
    if role == "unlocalized-scaffold":
        return "official_deposit_only_unlocalized_scaffold"
    if role == "assembled-molecule":
        return "official_deposit_only_assembled_molecule"
    return "official_deposit_only_other_role"


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    try:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=columns,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row.get(column, "") for column in columns})
    except OSError as error:
        raise ScopeAuditError(f"Cannot write {path}: {error}") from error


def output_file_metadata(path: Path) -> dict[str, object]:
    size, digest = sha256_file(path)
    return {"path": path.name, "bytes": size, "sha256": digest}


def build_audit(
    rows: list[ManifestRow], manifest: Path, staging: Path
) -> dict[str, object]:
    fasta_cache = FastaIndexCache()
    report_cache: dict[Path, AssemblyReportIndex] = {}
    report_consumers: dict[Path, set[str]] = defaultdict(set)
    summary_rows: list[dict[str, object]] = []
    match_rows: list[dict[str, object]] = []
    analysis_only_rows: list[dict[str, object]] = []
    deposit_only_rows: list[dict[str, object]] = []
    conflict_rows: list[dict[str, object]] = []

    for row in rows:
        analysis = fasta_cache.get(
            row.analysis_fasta, row.assembly_unit_id, "analysis_fasta"
        )
        deposit: FastaIndex | None = None
        report: AssemblyReportIndex | None = None
        report_by_id: dict[str, AssemblyReportRecord] | None = None
        matches: list[RecordMatch] = []
        conflicts: list[RecordConflict] = []
        analysis_only: list[FastaRecord]
        deposit_only: list[FastaRecord] = []

        if row.deposit_fasta is None:
            analysis_only = list(analysis.records)
        else:
            deposit = fasta_cache.get(
                row.deposit_fasta, row.assembly_unit_id, "deposit_fasta"
            )
            if row.assembly_report is not None:
                report_path = row.assembly_report.resolve(strict=True)
                report_consumers[report_path].add(row.assembly_unit_id)
                if report_path not in report_cache:
                    report_cache[report_path] = parse_assembly_report(report_path)
                report = report_cache[report_path]
            if row.pair_class == "official_same_accession_full_deposit":
                if report is None:
                    raise AssertionError("Official report validation was bypassed")
                report_by_id = verify_official_report(
                    deposit, report, row.assembly_unit_id
                )
            matches, analysis_only, deposit_only, conflicts = reconcile_records(
                analysis,
                deposit,
                row.assembly_unit_id,
                allow_same_id_sequence_mismatch=(
                    row.pair_class == "unverified_name_candidate"
                ),
            )

        if row.pair_class == "official_same_accession_full_deposit" and analysis_only:
            raise ScopeAuditError(
                f"Official full deposit does not contain all analysis records for "
                f"{row.assembly_unit_id}: "
                f"{[record.record_id for record in analysis_only[:5]]}"
            )
        if row.pair_class == "same_genome_matched_source" and (
            analysis_only or deposit_only
        ):
            raise ScopeAuditError(
                "Same-genome matched source has unmatched records for "
                f"{row.assembly_unit_id}: "
                f"analysis_only={len(analysis_only)}, deposit_only={len(deposit_only)}"
            )

        for match in sorted(matches, key=lambda item: item.analysis.record_id):
            match_rows.append(
                {
                    "assembly_unit_id": row.assembly_unit_id,
                    "pair_class": row.pair_class,
                    "analysis_record_id": match.analysis.record_id,
                    "deposit_record_id": match.deposit.record_id,
                    "match_method": match.method,
                    "sequence_bp": match.analysis.length,
                    "sequence_sha256": match.analysis.sequence_sha256,
                }
            )

        conflict_ids = {conflict.analysis.record_id for conflict in conflicts}
        analysis_only_reason = (
            "deposit_unavailable"
            if row.pair_class == "unavailable"
            else "no_exact_or_unique_sequence_match"
        )
        for record in sorted(analysis_only, key=lambda item: item.record_id):
            analysis_only_rows.append(
                {
                    "assembly_unit_id": row.assembly_unit_id,
                    "pair_class": row.pair_class,
                    "analysis_record_id": record.record_id,
                    "sequence_bp": record.length,
                    "sequence_sha256": record.sequence_sha256,
                    "reason": (
                        "same_id_sequence_mismatch"
                        if record.record_id in conflict_ids
                        else analysis_only_reason
                    ),
                }
            )

        for conflict in sorted(conflicts, key=lambda item: item.analysis.record_id):
            conflict_rows.append(
                {
                    "assembly_unit_id": row.assembly_unit_id,
                    "pair_class": row.pair_class,
                    "record_id": conflict.analysis.record_id,
                    "analysis_sequence_bp": conflict.analysis.length,
                    "analysis_sequence_sha256": conflict.analysis.sequence_sha256,
                    "deposit_sequence_bp": conflict.deposit.length,
                    "deposit_sequence_sha256": conflict.deposit.sequence_sha256,
                    "interpretation": "candidate_same_id_sequence_mismatch",
                }
            )

        official_unplaced_records = 0
        official_unplaced_bp = 0
        for record in sorted(deposit_only, key=lambda item: item.record_id):
            raw_role = ""
            if report_by_id is not None:
                raw_role = report_by_id[record.record_id].role
                interpretation = classify_official_deposit_only(raw_role)
                if raw_role == "unplaced-scaffold":
                    official_unplaced_records += 1
                    official_unplaced_bp += record.length
            else:
                interpretation = "candidate_only_scope_unknown"
            deposit_only_rows.append(
                {
                    "assembly_unit_id": row.assembly_unit_id,
                    "pair_class": row.pair_class,
                    "deposit_record_id": record.record_id,
                    "sequence_bp": record.length,
                    "sequence_sha256": record.sequence_sha256,
                    "assembly_report_role": raw_role,
                    "scope_interpretation": interpretation,
                }
            )

        exact_matches = [
            match for match in matches if match.method == "exact_id_and_sequence"
        ]
        hash_matches = [
            match for match in matches if match.method == "unique_sequence_hash"
        ]
        matched_bp = sum(match.analysis.length for match in matches)
        analysis_only_bp = sum(record.length for record in analysis_only)
        deposit_only_bp = sum(record.length for record in deposit_only)
        role_totals = report_role_totals(report_by_id)

        if row.pair_class == "unavailable":
            reconciliation_status = "deposit_unavailable"
        elif row.pair_class == "official_same_accession_full_deposit":
            reconciliation_status = (
                "official_full_scope_exact"
                if not deposit_only
                else "official_analysis_subset_reconciled"
            )
        elif row.pair_class == "same_genome_matched_source":
            reconciliation_status = "same_genome_exactly_reconciled"
        elif conflicts and not matches:
            reconciliation_status = "candidate_incompatible_same_id_sequence_mismatch"
        elif conflicts and analysis_only:
            reconciliation_status = "candidate_partial_with_same_id_sequence_mismatch"
        elif conflicts:
            reconciliation_status = (
                "candidate_contains_all_analysis_records_after_cross_id_reconciliation_with_name_conflicts"
            )
        elif analysis_only:
            reconciliation_status = "candidate_partial_reconciliation"
        else:
            reconciliation_status = "candidate_contains_all_analysis_records"

        official_fields: dict[str, object]
        if report_by_id is None:
            official_fields = {
                "official_unplaced_record_count": "",
                "official_unplaced_bp": "",
                "official_unplaced_bp_percent_of_deposit": "",
                "report_assembled_molecule_record_count": "",
                "report_assembled_molecule_bp": "",
                "report_unlocalized_scaffold_record_count": "",
                "report_unlocalized_scaffold_bp": "",
                "report_unplaced_scaffold_record_count": "",
                "report_unplaced_scaffold_bp": "",
                "report_other_role_record_count": "",
                "report_other_role_bp": "",
            }
        else:
            assert deposit is not None
            official_fields = {
                "official_unplaced_record_count": official_unplaced_records,
                "official_unplaced_bp": official_unplaced_bp,
                "official_unplaced_bp_percent_of_deposit": percent(
                    official_unplaced_bp, deposit.total_bp
                ),
                "report_assembled_molecule_record_count": role_totals[
                    "assembled_molecule_records"
                ],
                "report_assembled_molecule_bp": role_totals[
                    "assembled_molecule_bp"
                ],
                "report_unlocalized_scaffold_record_count": role_totals[
                    "unlocalized_scaffold_records"
                ],
                "report_unlocalized_scaffold_bp": role_totals[
                    "unlocalized_scaffold_bp"
                ],
                "report_unplaced_scaffold_record_count": role_totals[
                    "unplaced_scaffold_records"
                ],
                "report_unplaced_scaffold_bp": role_totals[
                    "unplaced_scaffold_bp"
                ],
                "report_other_role_record_count": role_totals["other_records"],
                "report_other_role_bp": role_totals["other_bp"],
            }

        summary: dict[str, object] = {
            "assembly_unit_id": row.assembly_unit_id,
            "biological_species": row.biological_species,
            "pair_class": row.pair_class,
            "scope_gap_authority": scope_gap_authority(row.pair_class),
            "analysis_accession": row.analysis_accession,
            "deposit_accession": row.deposit_accession,
            "analysis_scope": row.analysis_scope,
            "deposit_scope": row.deposit_scope,
            "analysis_fasta": str(analysis.path),
            "deposit_fasta": "" if deposit is None else str(deposit.path),
            "assembly_report": "" if report is None else str(report.path),
            "analysis_record_count": len(analysis.records),
            "analysis_total_bp": analysis.total_bp,
            "deposit_record_count": "" if deposit is None else len(deposit.records),
            "deposit_total_bp": "" if deposit is None else deposit.total_bp,
            "exact_id_sequence_match_count": len(exact_matches),
            "exact_id_sequence_match_bp": sum(
                match.analysis.length for match in exact_matches
            ),
            "sequence_hash_only_match_count": len(hash_matches),
            "sequence_hash_only_match_bp": sum(
                match.analysis.length for match in hash_matches
            ),
            "same_id_sequence_mismatch_count": len(conflicts),
            "same_id_sequence_mismatch_analysis_bp": sum(
                conflict.analysis.length for conflict in conflicts
            ),
            "same_id_sequence_mismatch_deposit_bp": sum(
                conflict.deposit.length for conflict in conflicts
            ),
            "matched_analysis_record_count": len(matches),
            "matched_analysis_bp": matched_bp,
            "analysis_bp_reconciled_percent": (
                "" if deposit is None else percent(matched_bp, analysis.total_bp)
            ),
            "analysis_only_record_count": len(analysis_only),
            "analysis_only_bp": analysis_only_bp,
            "deposit_only_record_count": len(deposit_only),
            "deposit_only_bp": deposit_only_bp,
            "deposit_only_bp_percent_of_deposit": (
                "" if deposit is None else percent(deposit_only_bp, deposit.total_bp)
            ),
            "record_reconciliation_status": reconciliation_status,
            "deposit_only_interpretation": deposit_only_interpretation(
                row.pair_class
            ),
            "evidence_note": row.evidence_note,
        }
        summary.update(official_fields)
        summary_rows.append(summary)

    summary_rows.sort(key=lambda item: str(item["assembly_unit_id"]))
    match_rows.sort(
        key=lambda item: (
            str(item["assembly_unit_id"]), str(item["analysis_record_id"])
        )
    )
    analysis_only_rows.sort(
        key=lambda item: (
            str(item["assembly_unit_id"]), str(item["analysis_record_id"])
        )
    )
    deposit_only_rows.sort(
        key=lambda item: (
            str(item["assembly_unit_id"]), str(item["deposit_record_id"])
        )
    )
    conflict_rows.sort(
        key=lambda item: (str(item["assembly_unit_id"]), str(item["record_id"]))
    )

    summary_path = staging / "paired_scope_summary.tsv"
    matches_path = staging / "record_matches.tsv"
    analysis_only_path = staging / "analysis_only_records.tsv"
    deposit_only_path = staging / "deposit_only_records.tsv"
    conflicts_path = staging / "conflicting_id_records.tsv"
    inputs_path = staging / "input_checksums.tsv"
    write_tsv(summary_path, SUMMARY_COLUMNS, summary_rows)
    write_tsv(matches_path, MATCH_COLUMNS, match_rows)
    write_tsv(analysis_only_path, ANALYSIS_ONLY_COLUMNS, analysis_only_rows)
    write_tsv(deposit_only_path, DEPOSIT_ONLY_COLUMNS, deposit_only_rows)
    write_tsv(conflicts_path, CONFLICT_COLUMNS, conflict_rows)

    input_rows: list[dict[str, object]] = []
    for path, index in sorted(
        fasta_cache.indexes.items(), key=lambda item: str(item[0])
    ):
        input_rows.append(
            {
                "input_type": ";".join(sorted(fasta_cache.input_types[path])),
                "consumer_assembly_units": ";".join(
                    sorted(fasta_cache.consumers[path])
                ),
                "resolved_path": str(path),
                "physical_size_bytes": index.physical_size,
                "physical_sha256": index.physical_sha256,
                "record_count": len(index.records),
                "logical_sequence_bp": index.total_bp,
                "logical_index_sha256": index.logical_index_sha256,
            }
        )
    for path, report in sorted(report_cache.items(), key=lambda item: str(item[0])):
        input_rows.append(
            {
                "input_type": "assembly_report",
                "consumer_assembly_units": ";".join(
                    sorted(report_consumers[path])
                ),
                "resolved_path": str(path),
                "physical_size_bytes": report.physical_size,
                "physical_sha256": report.physical_sha256,
                "record_count": len(report.records),
                "logical_sequence_bp": sum(
                    record.length for record in report.records
                ),
                "logical_index_sha256": report.logical_index_sha256,
            }
        )
    write_tsv(inputs_path, INPUT_COLUMNS, input_rows)

    manifest_size, manifest_sha256 = sha256_file(manifest.resolve(strict=True))
    output_paths = [
        summary_path,
        matches_path,
        analysis_only_path,
        deposit_only_path,
        conflicts_path,
        inputs_path,
    ]
    metadata = {
        "status": "completed",
        "script": str(Path(__file__).resolve()),
        "script_version": SCRIPT_VERSION,
        "script_sha256": sha256_file(Path(__file__).resolve())[1],
        "command": [sys.executable, *sys.argv],
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "execution": {
            "processes": 1,
            "worker_threads_started": 0,
            "fasta_cache_hits": fasta_cache.hits,
            "fasta_cache_misses": fasta_cache.misses,
            "unique_fasta_files_indexed": len(fasta_cache.indexes),
            "unique_assembly_reports_parsed": len(report_cache),
        },
        "manifest": {
            "path": str(manifest.resolve(strict=True)),
            "bytes": manifest_size,
            "sha256": manifest_sha256,
            "assembly_unit_rows": len(rows),
            "pair_class_counts": dict(
                sorted(Counter(row.pair_class for row in rows).items())
            ),
        },
        "matching_contract": {
            "fasta_id": "first whitespace-delimited token after >; case-sensitive",
            "sequence_normalization": "remove whitespace and uppercase; retain all other characters",
            "phase_1": "exact ID plus exact normalized sequence SHA-256",
            "phase_2": "unique length plus normalized sequence SHA-256",
            "ambiguous_hash_match": "fatal",
            "same_id_different_sequence": (
                "recorded only for unverified_name_candidate; fatal for all stronger evidence classes"
            ),
            "official_unplaced_authority": "only official_same_accession_full_deposit with an exact assembly-report reconciliation",
            "candidate_only_scope": "candidate_only_scope_unknown",
        },
        "row_totals": {
            "summary": len(summary_rows),
            "record_matches": len(match_rows),
            "analysis_only_records": len(analysis_only_rows),
            "deposit_only_records": len(deposit_only_rows),
            "conflicting_id_records": len(conflict_rows),
            "input_files": len(input_rows),
        },
        "outputs": [output_file_metadata(path) for path in output_paths],
    }
    metadata_path = staging / "run_metadata.json"
    try:
        metadata_path.write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    except OSError as error:
        raise ScopeAuditError(f"Cannot write {metadata_path}: {error}") from error
    return metadata


def publish(rows: list[ManifestRow], manifest: Path, output_dir: Path) -> dict[str, object]:
    """Build in a sibling directory and atomically install the completed run."""
    output = absolute_lexical(output_dir)
    if output.exists():
        raise ScopeAuditError(f"Output directory already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        staging = Path(
            tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent)
        )
    except OSError as error:
        raise ScopeAuditError(
            f"Cannot create staging directory beside {output}: {error}"
        ) from error
    try:
        metadata = build_audit(rows, manifest, staging)
        os.replace(staging, output)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return metadata


def main() -> int:
    args = parse_args()
    try:
        manifest = absolute_lexical(args.manifest)
        rows = load_manifest(manifest, args.expected_assembly_unit_count)
        metadata = publish(rows, manifest, args.output_dir)
    except ScopeAuditError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(
        "PAIRED ASSEMBLY-SCOPE AUDIT COMPLETE: "
        f"{metadata['manifest']['assembly_unit_rows']} assembly units; "
        f"{metadata['execution']['unique_fasta_files_indexed']} unique FASTA files; "
        f"cache_hits={metadata['execution']['fasta_cache_hits']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
