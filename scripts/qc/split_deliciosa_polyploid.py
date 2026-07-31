#!/usr/bin/env python3
"""Audit and split an *Actinidia deliciosa* A--F public assembly bundle.

The program is deliberately fail-closed.  It accepts one already resolved
parent bundle, verifies the physical SHA-256 of genome, GFF3, CDS, and protein
inputs, and partitions all four assets with one explicit mapping table.  A
downstream-compatible six-row manifest is published only when chromosome
counts, sequence-ID scope, coordinates, CDS/protein equality, and exact GFF3
identifier joins all pass.

Qinmei chromosome records are assigned by the ``_1`` ... ``_6`` suffix.  Its
known ``scf`` annotation records are retained under ``audit/unassigned`` and
reported, but never assigned to A--F because the matching scaffold FASTA was
not released in the public chromosome bundle.  ADM chromosome accessions are
assigned from mutually agreeing ``OriSeqID=ChrNNhN`` and ``Chromosome ANN``
header evidence; GFF3 sequence IDs and annotation ``Position=`` accessions are
then joined exactly to that genome map.

Unexpected or conflicting identifiers are retained and make the run BLOCKED.
The audit directory is still installed atomically for diagnosis, but the
accepted ``resolved_assembly_units.tsv`` is omitted and the CLI exits with
status 2.  Structural input errors (including checksum mismatch) publish
nothing and exit with status 1.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO
from urllib.parse import unquote


SCRIPT_VERSION = "1.0.0"
READ_CHUNK = 8 * 1024 * 1024
GZIP_MAGIC = b"\x1f\x8b"

MAPPING_COLUMNS = (
    "bundle_id",
    "partition_scheme",
    "partition_token",
    "partition_label",
    "assembly_unit_id",
    "biological_species",
    "individual_id",
    "ploidy",
    "expected_chromosome_count",
    "allow_unplaced_annotations",
)
RESOLVED_COLUMNS = (
    "assembly_unit_id",
    "accession",
    "genome",
    "gff",
    "cds",
    "protein",
    "genome_local_sha256",
    "gff_local_sha256",
    "cds_local_sha256",
    "protein_local_sha256",
)
ROLES = ("genome", "gff", "cds", "protein")
FASTA_ROLES = ("genome", "cds", "protein")
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
QINMEI_CHROMOSOME_RE = re.compile(r"^chr0?([1-9]|[12][0-9])_([1-6])$", re.I)
QINMEI_ANNOTATION_RE = re.compile(r"^Ad_0?([1-9]|[12][0-9])_([1-6])g", re.I)
ADM_ORIGINAL_RE = re.compile(r"OriSeqID=Chr0?([1-9]|[12][0-9])h([1-6])(?:\b|$)", re.I)
ADM_CHROMOSOME_RE = re.compile(r"Chromosome[ =]+([A-F])0?([1-9]|[12][0-9])(?:\b|$)", re.I)
ADM_POSITION_RE = re.compile(r"(?:^|[\s;])Position=([^\s;]+)", re.I)
ADM_ANNOTATION_ID_RE = re.compile(r"^Achdmh([1-6])c0?([1-9]|[12][0-9])", re.I)
ADM_PROTEIN_LINK_RE = re.compile(r"(?:^|\s)Protein=([^\s]+)", re.I)
ADM_MRNA_LINK_RE = re.compile(r"(?:^|\s)mRNA=([^\s]+)", re.I)
UNPLACED_RE = re.compile(r"(?:scf|scaffold|contig|unplaced|unloc)", re.I)


class PartitionInputError(RuntimeError):
    """Raised when an input contract is invalid and no output may be published."""


@dataclass(frozen=True)
class Partition:
    bundle_id: str
    scheme: str
    token: str
    label: str
    unit_id: str
    species: str
    individual: str
    ploidy: str
    expected_chromosomes: int
    allow_unplaced: bool


@dataclass(frozen=True)
class Asset:
    role: str
    path: Path
    expected_sha256: str
    actual_sha256: str


@dataclass(frozen=True)
class Assignment:
    unit_id: str | None
    classification: str
    evidence: str
    chromosome_number: int | None


@dataclass(frozen=True)
class RunResult:
    status: str
    output_dir: Path
    error_count: int
    warning_count: int


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(READ_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_bool(value: str, *, location: str) -> bool:
    lowered = value.strip().lower()
    if lowered not in {"true", "false"}:
        raise PartitionInputError(f"{location}: expected true or false, found {value!r}")
    return lowered == "true"


def read_mapping(path: Path, bundle_id: str) -> tuple[Partition, ...]:
    resolved = path.expanduser().resolve()
    try:
        handle = resolved.open("r", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise PartitionInputError(f"Cannot open mapping table {resolved}: {error}") from error
    rows: list[Partition] = []
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PartitionInputError(f"Mapping table has no header: {resolved}")
        missing = [column for column in MAPPING_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise PartitionInputError(
                f"Mapping table is missing columns: {', '.join(missing)}"
            )
        for line_number, row in enumerate(reader, start=2):
            if (row.get("bundle_id") or "").strip() != bundle_id:
                continue
            location = f"{resolved.name}:{line_number}"
            token = (row.get("partition_token") or "").strip()
            label = (row.get("partition_label") or "").strip()
            unit_id = (row.get("assembly_unit_id") or "").strip()
            if token not in {"1", "2", "3", "4", "5", "6"}:
                raise PartitionInputError(f"{location}: invalid partition token {token!r}")
            if label not in {"A", "B", "C", "D", "E", "F"}:
                raise PartitionInputError(f"{location}: invalid partition label {label!r}")
            if not SAFE_ID.fullmatch(unit_id):
                raise PartitionInputError(f"{location}: unsafe assembly_unit_id {unit_id!r}")
            try:
                expected = int((row.get("expected_chromosome_count") or "").strip())
            except ValueError as error:
                raise PartitionInputError(
                    f"{location}: expected_chromosome_count must be an integer"
                ) from error
            if expected < 1:
                raise PartitionInputError(f"{location}: expected chromosome count must be positive")
            rows.append(
                Partition(
                    bundle_id=bundle_id,
                    scheme=(row.get("partition_scheme") or "").strip(),
                    token=token,
                    label=label,
                    unit_id=unit_id,
                    species=(row.get("biological_species") or "").strip(),
                    individual=(row.get("individual_id") or "").strip(),
                    ploidy=(row.get("ploidy") or "").strip(),
                    expected_chromosomes=expected,
                    allow_unplaced=parse_bool(
                        row.get("allow_unplaced_annotations") or "", location=location
                    ),
                )
            )
    if len(rows) != 6:
        raise PartitionInputError(
            f"{bundle_id}: mapping must contain exactly six A--F rows, found {len(rows)}"
        )
    for attribute in ("token", "label", "unit_id"):
        values = [getattr(row, attribute) for row in rows]
        if len(set(values)) != 6:
            raise PartitionInputError(f"{bundle_id}: duplicate mapping {attribute}")
    if {row.token for row in rows} != set("123456"):
        raise PartitionInputError(f"{bundle_id}: mapping tokens must be exactly 1--6")
    if {row.label for row in rows} != set("ABCDEF"):
        raise PartitionInputError(f"{bundle_id}: mapping labels must be exactly A--F")
    schemes = {row.scheme for row in rows}
    if schemes not in ({"qinmei_suffix"}, {"adm_header"}):
        raise PartitionInputError(
            f"{bundle_id}: partition_scheme must be uniformly qinmei_suffix or adm_header"
        )
    invariant_fields = ("species", "individual", "ploidy", "expected_chromosomes", "allow_unplaced")
    for attribute in invariant_fields:
        if len({getattr(row, attribute) for row in rows}) != 1:
            raise PartitionInputError(f"{bundle_id}: inconsistent mapping field {attribute}")
    return tuple(sorted(rows, key=lambda row: int(row.token)))


def _resolve_manifest_asset(raw: str, manifest_dir: Path, role: str) -> Path:
    if not raw.strip():
        raise PartitionInputError(f"Resolved bundle is missing required {role} path")
    candidate = Path(raw.strip()).expanduser()
    if not candidate.is_absolute():
        candidate = manifest_dir / candidate
    resolved = candidate.resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise PartitionInputError(f"Resolved {role} asset is missing or empty: {resolved}")
    return resolved


def read_assets(manifest_path: Path, bundle_id: str) -> tuple[str, dict[str, Asset]]:
    resolved = manifest_path.expanduser().resolve()
    try:
        handle = resolved.open("r", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise PartitionInputError(f"Cannot open resolved manifest {resolved}: {error}") from error
    matches: list[dict[str, str]] = []
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PartitionInputError(f"Resolved manifest has no header: {resolved}")
        missing = [column for column in RESOLVED_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise PartitionInputError(
                f"Resolved manifest is missing columns: {', '.join(missing)}"
            )
        for row in reader:
            if (row.get("assembly_unit_id") or "").strip() == bundle_id:
                matches.append({key: (value or "").strip() for key, value in row.items()})
    if len(matches) != 1:
        raise PartitionInputError(
            f"Resolved manifest must contain exactly one {bundle_id!r} row, found {len(matches)}"
        )
    row = matches[0]
    assets: dict[str, Asset] = {}
    for role in ROLES:
        expected = row[f"{role}_local_sha256"].lower()
        if not SHA256_RE.fullmatch(expected):
            raise PartitionInputError(
                f"{bundle_id}: {role}_local_sha256 must be a complete SHA-256"
            )
        asset_path = _resolve_manifest_asset(row[role], resolved.parent, role)
        actual = sha256_file(asset_path)
        if actual != expected:
            raise PartitionInputError(
                f"{bundle_id}: {role} checksum mismatch; expected {expected}, observed {actual}"
            )
        assets[role] = Asset(role, asset_path, expected, actual)
    return row.get("accession", ""), assets


@contextmanager
def open_text_auto(path: Path) -> Iterator[TextIO]:
    """Open gzip by magic bytes, not by filename suffix."""
    with path.open("rb") as probe:
        magic = probe.read(2)
    if magic == GZIP_MAGIC:
        with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline=None) as handle:
            yield handle
    else:
        with path.open("rt", encoding="utf-8", errors="strict", newline=None) as handle:
            yield handle


@contextmanager
def deterministic_gzip_writer(path: Path) -> Iterator[TextIO]:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as text:
                yield text


def parse_attributes(text: str) -> dict[str, list[str]]:
    values: dict[str, list[str]] = defaultdict(list)
    for item in text.strip().split(";"):
        if not item or "=" not in item:
            continue
        key, raw = item.split("=", 1)
        for value in raw.split(","):
            decoded = unquote(value.strip())
            if decoded:
                values[key.strip()].append(decoded)
    return dict(values)


def _unplaced(value: str) -> bool:
    return bool(UNPLACED_RE.search(value))


class SplitAudit:
    def __init__(
        self,
        *,
        partitions: tuple[Partition, ...],
        accession: str,
        assets: dict[str, Asset],
        stage: Path,
        manifest_name: str,
    ) -> None:
        self.partitions = partitions
        self.accession = accession
        self.assets = assets
        self.stage = stage
        self.manifest_name = manifest_name
        self.by_token = {row.token: row for row in partitions}
        self.by_label = {row.label: row for row in partitions}
        self.by_unit = {row.unit_id: row for row in partitions}
        self.scheme = partitions[0].scheme
        self.allow_unplaced = partitions[0].allow_unplaced
        self.genome_unit_by_seqid: dict[str, str] = {}
        self.genome_length_by_seqid: dict[str, int] = {}
        self.genome_chromosomes: dict[str, dict[int, str]] = defaultdict(dict)
        self.fasta_ids: dict[str, dict[str, set[str]]] = {
            role: defaultdict(set) for role in FASTA_ROLES
        }
        self.unassigned_fasta_ids: dict[str, set[str]] = defaultdict(set)
        self.gff_identifiers: dict[str, set[str]] = defaultdict(set)
        self.unassigned_gff_identifiers: set[str] = set()
        self.adm_cds_to_protein: dict[str, dict[str, str]] = defaultdict(dict)
        self.adm_protein_to_cds: dict[str, dict[str, str]] = defaultdict(dict)
        self.gff_seqids: dict[str, set[str]] = defaultdict(set)
        self.counts: Counter[tuple[str, str, str]] = Counter()
        self.issues: list[dict[str, str]] = []
        self.fasta_audit_rows: list[dict[str, str]] = []
        self.unassigned_rows: list[dict[str, str]] = []

    def issue(self, severity: str, code: str, unit: str, role: str, detail: str) -> None:
        self.issues.append(
            {
                "severity": severity,
                "code": code,
                "assembly_unit_id": unit,
                "asset_role": role,
                "detail": detail,
            }
        )

    def unit_for_token(self, token: str) -> Partition | None:
        return self.by_token.get(token)

    def assign_genome(self, record_id: str, header: str) -> Assignment:
        if self.scheme == "qinmei_suffix":
            match = QINMEI_CHROMOSOME_RE.fullmatch(record_id)
            if not match:
                return Assignment(None, "unmatched", "genome ID does not match chrNN_1--6", None)
            chrom, token = int(match.group(1)), match.group(2)
            return Assignment(self.by_token[token].unit_id, "assigned", "qinmei_seqid_suffix", chrom)

        evidence: list[tuple[str, int, str]] = []
        for match in ADM_ORIGINAL_RE.finditer(header):
            evidence.append((match.group(2), int(match.group(1)), "OriSeqID"))
        for match in ADM_CHROMOSOME_RE.finditer(header):
            partition = self.by_label.get(match.group(1).upper())
            if partition is not None:
                evidence.append((partition.token, int(match.group(2)), "Chromosome"))
        if not evidence:
            return Assignment(None, "unmatched", "ADM genome header lacks partition evidence", None)
        units = {self.by_token[token].unit_id for token, _chrom, _source in evidence}
        chromosomes = {chrom for _token, chrom, _source in evidence}
        if len(units) != 1 or len(chromosomes) != 1:
            return Assignment(None, "conflicting", "ADM genome header evidence conflicts", None)
        return Assignment(
            next(iter(units)),
            "assigned",
            "+".join(sorted({source for _token, _chrom, source in evidence})),
            next(iter(chromosomes)),
        )

    def assign_annotation_fasta(self, record_id: str, header: str) -> Assignment:
        if self.scheme == "qinmei_suffix":
            match = QINMEI_ANNOTATION_RE.match(record_id)
            if match:
                chrom, token = int(match.group(1)), match.group(2)
                return Assignment(
                    self.by_token[token].unit_id, "assigned", "qinmei_annotation_id", chrom
                )
            if _unplaced(record_id) or _unplaced(header):
                return Assignment(None, "unplaced", "recognized Qinmei scaffold annotation", None)
            return Assignment(None, "unmatched", "annotation ID does not match Qinmei rule", None)

        evidence: list[tuple[str, int | None, str]] = []
        for match in ADM_POSITION_RE.finditer(header):
            seqid = match.group(1).split(":", 1)[0]
            unit = self.genome_unit_by_seqid.get(seqid)
            if unit:
                token = self.by_unit[unit].token
                evidence.append((token, None, "Position"))
        for match in ADM_ORIGINAL_RE.finditer(header):
            evidence.append((match.group(2), int(match.group(1)), "OriSeqID"))
        identifier_match = ADM_ANNOTATION_ID_RE.match(record_id)
        if identifier_match:
            evidence.append(
                (identifier_match.group(1), int(identifier_match.group(2)), "record_id")
            )
        if not evidence:
            if _unplaced(record_id) or _unplaced(header):
                return Assignment(None, "unplaced", "recognized ADM unplaced annotation", None)
            return Assignment(None, "unmatched", "ADM annotation lacks exact partition evidence", None)
        units = {self.by_token[token].unit_id for token, _chrom, _source in evidence}
        chromosomes = {chrom for _token, chrom, _source in evidence if chrom is not None}
        if len(units) != 1 or len(chromosomes) > 1:
            return Assignment(None, "conflicting", "ADM annotation evidence conflicts", None)
        return Assignment(
            next(iter(units)),
            "assigned",
            "+".join(sorted({source for _token, _chrom, source in evidence})),
            next(iter(chromosomes)) if chromosomes else None,
        )

    def _output_path(self, unit: str, role: str) -> Path:
        suffix = "gff3.gz" if role == "gff" else f"{role}.fa.gz"
        return self.stage / "assembly_units" / unit / f"{unit}.{suffix}"

    def _unassigned_path(self, role: str) -> Path:
        suffix = "gff3.gz" if role == "gff" else f"{role}.fa.gz"
        return self.stage / "audit" / "unassigned" / f"{self.partitions[0].bundle_id}.{suffix}"

    def record_adm_cross_link(self, role: str, unit: str, record_id: str, header: str) -> None:
        """Freeze exact CDS--protein accession links declared in ADM headers."""
        if self.scheme != "adm_header" or unit == "UNASSIGNED":
            return
        if role == "cds":
            matches = ADM_PROTEIN_LINK_RE.findall(header)
            mapping = self.adm_cds_to_protein[unit]
            label = "Protein"
        elif role == "protein":
            matches = ADM_MRNA_LINK_RE.findall(header)
            mapping = self.adm_protein_to_cds[unit]
            label = "mRNA"
        else:
            return
        unique = sorted(set(matches))
        if len(unique) != 1:
            self.issue(
                "ERROR",
                "missing_or_ambiguous_adm_cross_link",
                unit,
                role,
                f"{record_id}: expected one exact {label}= accession, observed {unique}",
            )
            return
        previous = mapping.get(record_id)
        if previous is not None and previous != unique[0]:
            self.issue(
                "ERROR",
                "conflicting_adm_cross_link",
                unit,
                role,
                f"{record_id}: {previous} versus {unique[0]}",
            )
            return
        mapping[record_id] = unique[0]

    def split_fasta(self, role: str) -> None:
        asset = self.assets[role]
        seen: set[str] = set()
        with ExitStack() as stack:
            writers = {
                unit: stack.enter_context(deterministic_gzip_writer(self._output_path(unit, role)))
                for unit in self.by_unit
            }
            unassigned_writer = stack.enter_context(
                deterministic_gzip_writer(self._unassigned_path(role))
            )
            with open_text_auto(asset.path) as source:
                record_id: str | None = None
                header = ""
                assignment: Assignment | None = None
                writer: TextIO | None = None
                digest = hashlib.sha256()
                sequence_bp = 0
                wrap_buffer = ""

                def flush_sequence(final: bool = False) -> None:
                    nonlocal wrap_buffer
                    assert writer is not None
                    while len(wrap_buffer) >= 80:
                        writer.write(wrap_buffer[:80] + "\n")
                        wrap_buffer = wrap_buffer[80:]
                    if final and wrap_buffer:
                        writer.write(wrap_buffer + "\n")
                        wrap_buffer = ""

                def finish_record() -> None:
                    nonlocal record_id, header, assignment, writer, digest, sequence_bp, wrap_buffer
                    if record_id is None or assignment is None or writer is None:
                        return
                    flush_sequence(final=True)
                    unit = assignment.unit_id or "UNASSIGNED"
                    self.counts[(unit, role, assignment.classification)] += 1
                    if record_id in seen:
                        self.issue("ERROR", "duplicate_fasta_id", unit, role, record_id)
                    seen.add(record_id)
                    self.fasta_ids[role][unit].add(record_id)
                    self.record_adm_cross_link(role, unit, record_id, header)
                    if assignment.unit_id is None:
                        self.unassigned_fasta_ids[role].add(record_id)
                        self.unassigned_rows.append(
                            {
                                "asset_role": role,
                                "record_kind": "FASTA_record",
                                "record_id_or_seqid": record_id,
                                "classification": assignment.classification,
                                "retained_file": str(self._unassigned_path(role).relative_to(self.stage)),
                                "detail": assignment.evidence,
                            }
                        )
                        if assignment.classification != "unplaced" or not self.allow_unplaced:
                            self.issue(
                                "ERROR",
                                f"{assignment.classification}_{role}_record",
                                "UNASSIGNED",
                                role,
                                f"{record_id}: {assignment.evidence}",
                            )
                    self.fasta_audit_rows.append(
                        {
                            "asset_role": role,
                            "record_id": record_id,
                            "assembly_unit_id": unit,
                            "classification": assignment.classification,
                            "assignment_evidence": assignment.evidence,
                            "chromosome_number": str(assignment.chromosome_number or ""),
                            "sequence_bp": str(sequence_bp),
                            "sequence_sha256": digest.hexdigest(),
                        }
                    )
                    if role == "genome" and assignment.unit_id is not None:
                        if record_id in self.genome_unit_by_seqid:
                            self.issue("ERROR", "duplicate_genome_seqid", unit, role, record_id)
                        self.genome_unit_by_seqid[record_id] = assignment.unit_id
                        self.genome_length_by_seqid[record_id] = sequence_bp
                        if assignment.chromosome_number is None:
                            self.issue("ERROR", "missing_chromosome_number", unit, role, record_id)
                        elif assignment.chromosome_number in self.genome_chromosomes[unit]:
                            self.issue(
                                "ERROR",
                                "duplicate_chromosome_number",
                                unit,
                                role,
                                str(assignment.chromosome_number),
                            )
                        else:
                            self.genome_chromosomes[unit][assignment.chromosome_number] = record_id

                for raw_line in source:
                    if raw_line.startswith(">"):
                        finish_record()
                        header = raw_line[1:].strip()
                        record_id = header.split(None, 1)[0] if header else ""
                        if not record_id:
                            raise PartitionInputError(f"{role} FASTA contains an empty header")
                        assignment = (
                            self.assign_genome(record_id, header)
                            if role == "genome"
                            else self.assign_annotation_fasta(record_id, header)
                        )
                        writer = (
                            writers[assignment.unit_id]
                            if assignment.unit_id is not None
                            else unassigned_writer
                        )
                        writer.write(">" + header + "\n")
                        digest = hashlib.sha256()
                        sequence_bp = 0
                        wrap_buffer = ""
                        continue
                    if record_id is None:
                        if raw_line.strip():
                            raise PartitionInputError(f"{role} FASTA contains sequence before header")
                        continue
                    normalized = "".join(raw_line.split()).upper()
                    if normalized:
                        digest.update(normalized.encode("ascii"))
                        sequence_bp += len(normalized)
                        wrap_buffer += normalized
                        flush_sequence()
                finish_record()
        if not seen:
            raise PartitionInputError(f"{role} FASTA contains no records: {asset.path}")

    def _write_gff_global_line(self, writers: dict[str, TextIO], unassigned: TextIO, line: str) -> None:
        for writer in writers.values():
            writer.write(line)
        unassigned.write(line)

    def split_gff(self) -> None:
        asset = self.assets["gff"]
        unassigned_seqid_counts: Counter[str] = Counter()
        with ExitStack() as stack:
            writers = {
                unit: stack.enter_context(deterministic_gzip_writer(self._output_path(unit, "gff")))
                for unit in self.by_unit
            }
            unassigned = stack.enter_context(
                deterministic_gzip_writer(self._unassigned_path("gff"))
            )
            self._write_gff_global_line(writers, unassigned, "##gff-version 3\n")
            with open_text_auto(asset.path) as source:
                for line_number, raw_line in enumerate(source, start=1):
                    line = raw_line if raw_line.endswith("\n") else raw_line + "\n"
                    stripped = raw_line.strip()
                    if not stripped or stripped == "##gff-version 3":
                        continue
                    if stripped == "##FASTA":
                        raise PartitionInputError(
                            "Embedded FASTA in GFF3 is unsupported; provide a feature-only GFF3"
                        )
                    if stripped.startswith("##sequence-region"):
                        fields = stripped.split()
                        if len(fields) != 4:
                            unassigned.write(line)
                            self.issue(
                                "ERROR", "malformed_sequence_region", "UNASSIGNED", "gff", f"line {line_number}"
                            )
                            continue
                        seqid = fields[1]
                        unit = self.genome_unit_by_seqid.get(seqid)
                        if unit is None:
                            unassigned.write(line)
                            unassigned_seqid_counts[seqid] += 1
                            classification = "unplaced" if _unplaced(seqid) else "unmatched"
                            if classification != "unplaced" or not self.allow_unplaced:
                                self.issue(
                                    "ERROR", f"{classification}_gff_seqid", "UNASSIGNED", "gff", seqid
                                )
                        else:
                            writers[unit].write(line)
                        continue
                    if stripped.startswith("#"):
                        self._write_gff_global_line(writers, unassigned, line)
                        continue
                    fields = raw_line.rstrip("\r\n").split("\t")
                    if len(fields) != 9:
                        unassigned.write(line)
                        self.issue(
                            "ERROR", "malformed_gff_row", "UNASSIGNED", "gff", f"line {line_number}"
                        )
                        continue
                    seqid = fields[0]
                    unit = self.genome_unit_by_seqid.get(seqid)
                    attributes = parse_attributes(fields[8])
                    identifiers = {
                        value
                        for key in (
                            "ID",
                            "Parent",
                            "protein_id",
                            "transcript_id",
                            "gene_id",
                            "Accession",
                            "Parent_Accession",
                            "Protein_Accession",
                        )
                        for value in attributes.get(key, [])
                    }
                    if unit is None:
                        unassigned.write(line)
                        unassigned_seqid_counts[seqid] += 1
                        self.unassigned_gff_identifiers.update(identifiers)
                        classification = "unplaced" if _unplaced(seqid) else "unmatched"
                        if classification != "unplaced" or not self.allow_unplaced:
                            self.issue(
                                "ERROR", f"{classification}_gff_seqid", "UNASSIGNED", "gff", seqid
                            )
                        continue
                    writers[unit].write(line)
                    self.counts[(unit, "gff", "feature")] += 1
                    self.gff_seqids[unit].add(seqid)
                    self.gff_identifiers[unit].update(identifiers)
                    try:
                        start, end = int(fields[3]), int(fields[4])
                    except ValueError:
                        self.issue(
                            "ERROR", "nonnumeric_gff_coordinate", unit, "gff", f"line {line_number}"
                        )
                    else:
                        length = self.genome_length_by_seqid[seqid]
                        if start < 1 or end < start or end > length:
                            self.issue(
                                "ERROR",
                                "gff_coordinate_out_of_bounds",
                                unit,
                                "gff",
                                f"line {line_number}: {seqid}:{start}-{end}, length={length}",
                            )
        for seqid, count in sorted(unassigned_seqid_counts.items()):
            classification = "unplaced" if _unplaced(seqid) else "unmatched"
            self.unassigned_rows.append(
                {
                    "asset_role": "gff",
                    "record_kind": "GFF_seqid",
                    "record_id_or_seqid": seqid,
                    "classification": classification,
                    "retained_file": str(self._unassigned_path("gff").relative_to(self.stage)),
                    "detail": f"{count} retained rows/directives",
                }
            )

    def validate(self) -> list[dict[str, str]]:
        joins: list[dict[str, str]] = []
        expected_numbers_by_unit = {
            unit: set(range(1, partition.expected_chromosomes + 1))
            for unit, partition in self.by_unit.items()
        }
        for unit, partition in self.by_unit.items():
            observed_numbers = set(self.genome_chromosomes[unit])
            expected_numbers = expected_numbers_by_unit[unit]
            if observed_numbers != expected_numbers:
                self.issue(
                    "ERROR",
                    "chromosome_set_mismatch",
                    unit,
                    "genome",
                    f"expected={sorted(expected_numbers)} observed={sorted(observed_numbers)}",
                )
            if not self.gff_seqids[unit]:
                self.issue("ERROR", "empty_partition_gff", unit, "gff", "no assigned features")
            protein = self.fasta_ids["protein"][unit]
            cds = self.fasta_ids["cds"][unit]
            gff = self.gff_identifiers[unit]
            if self.scheme == "qinmei_suffix":
                comparisons = (
                    ("protein_without_cds", protein - cds),
                    ("cds_without_protein", cds - protein),
                    ("protein_without_exact_gff_identifier", protein - gff),
                    ("cds_without_exact_gff_identifier", cds - gff),
                )
            else:
                cds_links = self.adm_cds_to_protein[unit]
                protein_links = self.adm_protein_to_cds[unit]
                missing_proteins = {
                    f"{cds_id}->{protein_id}"
                    for cds_id, protein_id in cds_links.items()
                    if protein_id not in protein
                }
                missing_cds = {
                    f"{protein_id}->{cds_id}"
                    for protein_id, cds_id in protein_links.items()
                    if cds_id not in cds
                }
                reciprocal_conflicts = {
                    f"{cds_id}->{protein_id}->{protein_links.get(protein_id, '<missing>')}"
                    for cds_id, protein_id in cds_links.items()
                    if protein_links.get(protein_id) != cds_id
                }
                comparisons = (
                    ("adm_cds_without_protein_link", cds - set(cds_links)),
                    ("adm_protein_without_mrna_link", protein - set(protein_links)),
                    ("adm_cross_link_missing_protein_record", missing_proteins),
                    ("adm_cross_link_missing_cds_record", missing_cds),
                    ("adm_nonreciprocal_cds_protein_link", reciprocal_conflicts),
                    ("protein_without_exact_gff_identifier", protein - gff),
                    ("cds_without_exact_gff_identifier", cds - gff),
                )
            for relation, identifiers in comparisons:
                for identifier in sorted(identifiers):
                    joins.append(
                        {
                            "assembly_unit_id": unit,
                            "relation": relation,
                            "record_id": identifier,
                        }
                    )
                if identifiers:
                    self.issue(
                        "ERROR",
                        relation,
                        unit,
                        "identifier_join",
                        f"{len(identifiers)} exact identifiers",
                    )
            for role, values in (("genome", self.fasta_ids["genome"][unit]), ("cds", cds), ("protein", protein)):
                if not values:
                    self.issue("ERROR", f"empty_partition_{role}", unit, role, "no assigned records")

        unplaced_protein = self.unassigned_fasta_ids["protein"]
        unplaced_cds = self.unassigned_fasta_ids["cds"]
        for relation, identifiers in (
            ("unassigned_protein_without_cds", unplaced_protein - unplaced_cds),
            ("unassigned_cds_without_protein", unplaced_cds - unplaced_protein),
            (
                "unassigned_protein_without_exact_gff_identifier",
                unplaced_protein - self.unassigned_gff_identifiers,
            ),
            (
                "unassigned_cds_without_exact_gff_identifier",
                unplaced_cds - self.unassigned_gff_identifiers,
            ),
        ):
            for identifier in sorted(identifiers):
                joins.append(
                    {"assembly_unit_id": "UNASSIGNED", "relation": relation, "record_id": identifier}
                )
            if identifiers:
                self.issue(
                    "ERROR",
                    relation,
                    "UNASSIGNED",
                    "identifier_join",
                    f"{len(identifiers)} exact identifiers",
                )

        if self.unassigned_rows and self.allow_unplaced:
            recognized = sum(row["classification"] == "unplaced" for row in self.unassigned_rows)
            if recognized:
                self.issue(
                    "WARNING",
                    "known_unplaced_annotations_retained",
                    "UNASSIGNED",
                    "annotation",
                    f"{recognized} audit records retained outside A--F",
                )
        return joins

    @staticmethod
    def _write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def write_reports(self, joins: list[dict[str, str]]) -> tuple[int, int]:
        audit = self.stage / "audit"
        errors = sum(row["severity"] == "ERROR" for row in self.issues)
        warnings = sum(row["severity"] == "WARNING" for row in self.issues)
        self._write_tsv(
            audit / "fasta_records.tsv",
            (
                "asset_role",
                "record_id",
                "assembly_unit_id",
                "classification",
                "assignment_evidence",
                "chromosome_number",
                "sequence_bp",
                "sequence_sha256",
            ),
            self.fasta_audit_rows,
        )
        self._write_tsv(
            audit / "unassigned_records.tsv",
            (
                "asset_role",
                "record_kind",
                "record_id_or_seqid",
                "classification",
                "retained_file",
                "detail",
            ),
            self.unassigned_rows,
        )
        self._write_tsv(
            audit / "id_join_exceptions.tsv",
            ("assembly_unit_id", "relation", "record_id"),
            joins,
        )
        self._write_tsv(
            audit / "validation_issues.tsv",
            ("severity", "code", "assembly_unit_id", "asset_role", "detail"),
            self.issues,
        )
        summary_rows: list[dict[str, str]] = []
        for partition in self.partitions:
            unit = partition.unit_id
            unit_errors = sum(
                row["severity"] == "ERROR" and row["assembly_unit_id"] in {unit, "UNASSIGNED"}
                for row in self.issues
            )
            summary_rows.append(
                {
                    "assembly_unit_id": unit,
                    "partition_label": partition.label,
                    "partition_token": partition.token,
                    "genome_record_count": str(len(self.fasta_ids["genome"][unit])),
                    "gff_sequence_id_count": str(len(self.gff_seqids[unit])),
                    "gff_feature_count": str(self.counts[(unit, "gff", "feature")]),
                    "cds_record_count": str(len(self.fasta_ids["cds"][unit])),
                    "protein_record_count": str(len(self.fasta_ids["protein"][unit])),
                    "exact_gff_identifier_count": str(len(self.gff_identifiers[unit])),
                    "status": "PASS" if unit_errors == 0 else "BLOCKED",
                }
            )
        self._write_tsv(
            audit / "partition_summary.tsv",
            (
                "assembly_unit_id",
                "partition_label",
                "partition_token",
                "genome_record_count",
                "gff_sequence_id_count",
                "gff_feature_count",
                "cds_record_count",
                "protein_record_count",
                "exact_gff_identifier_count",
                "status",
            ),
            summary_rows,
        )
        return errors, warnings

    def write_manifest(self) -> None:
        rows: list[dict[str, str]] = []
        bundle = self.partitions[0].bundle_id
        for partition in self.partitions:
            unit = partition.unit_id
            rows.append(
                {
                    "sample": unit,
                    "current_or_alternative": "alternative",
                    "accession": self.accession,
                    "genome": str(self._output_path(unit, "genome").relative_to(self.stage)),
                    "gff": str(self._output_path(unit, "gff").relative_to(self.stage)),
                    "protein": str(self._output_path(unit, "protein").relative_to(self.stage)),
                    "source_url": "",
                    "assembly_unit_id": unit,
                    "biological_species": partition.species,
                    "individual_id": partition.individual,
                    "haplotype_or_subgenome": partition.label,
                    "ploidy": partition.ploidy,
                    "source_bundle_id": bundle,
                    "parent_bundle_id": bundle,
                    "partition_token": partition.token,
                    "partition_scheme": partition.scheme,
                    "cds": str(self._output_path(unit, "cds").relative_to(self.stage)),
                    "include_qc": "true",
                    "include_gene_loss": "false",
                    "include_species_tree": "false",
                    "selection_status": "candidate_pending_QC_and_JCVI",
                }
            )
        self._write_tsv(
            self.stage / "resolved_assembly_units.tsv",
            (
                "sample",
                "current_or_alternative",
                "accession",
                "genome",
                "gff",
                "protein",
                "source_url",
                "assembly_unit_id",
                "biological_species",
                "individual_id",
                "haplotype_or_subgenome",
                "ploidy",
                "source_bundle_id",
                "parent_bundle_id",
                "partition_token",
                "partition_scheme",
                "cds",
                "include_qc",
                "include_gene_loss",
                "include_species_tree",
                "selection_status",
            ),
            rows,
        )

    def write_checksums(self) -> None:
        rows: list[dict[str, str]] = []
        for role, asset in sorted(self.assets.items()):
            rows.append(
                {
                    "scope": "input",
                    "asset_role": role,
                    "assembly_unit_id": self.partitions[0].bundle_id,
                    "relative_path_or_filename": asset.path.name,
                    "bytes": str(asset.path.stat().st_size),
                    "sha256": asset.actual_sha256,
                    "expected_sha256": asset.expected_sha256,
                    "verification": "PASS",
                }
            )
        for role in ROLES:
            for unit in self.by_unit:
                path = self._output_path(unit, role)
                rows.append(
                    {
                        "scope": "derived",
                        "asset_role": role,
                        "assembly_unit_id": unit,
                        "relative_path_or_filename": str(path.relative_to(self.stage)),
                        "bytes": str(path.stat().st_size),
                        "sha256": sha256_file(path),
                        "expected_sha256": "",
                        "verification": "CALCULATED",
                    }
                )
            path = self._unassigned_path(role)
            rows.append(
                {
                    "scope": "retained_unassigned",
                    "asset_role": role,
                    "assembly_unit_id": "UNASSIGNED",
                    "relative_path_or_filename": str(path.relative_to(self.stage)),
                    "bytes": str(path.stat().st_size),
                    "sha256": sha256_file(path),
                    "expected_sha256": "",
                    "verification": "CALCULATED",
                }
            )
        self._write_tsv(
            self.stage / "audit" / "checksums.tsv",
            (
                "scope",
                "asset_role",
                "assembly_unit_id",
                "relative_path_or_filename",
                "bytes",
                "sha256",
                "expected_sha256",
                "verification",
            ),
            rows,
        )

    def execute(self) -> tuple[int, int]:
        self.split_fasta("genome")
        self.split_gff()
        self.split_fasta("cds")
        self.split_fasta("protein")
        joins = self.validate()
        errors, warnings = self.write_reports(joins)
        if errors == 0:
            self.write_manifest()
        self.write_checksums()
        metadata = {
            "schema_version": 1,
            "script_version": SCRIPT_VERSION,
            "bundle_id": self.partitions[0].bundle_id,
            "partition_scheme": self.scheme,
            "partition_count": len(self.partitions),
            "input_manifest_filename": self.manifest_name,
            "status": "PASS" if errors == 0 else "BLOCKED",
            "error_count": errors,
            "warning_count": warnings,
            "accepted_manifest": "resolved_assembly_units.tsv" if errors == 0 else None,
        }
        path = self.stage / "audit" / "run_metadata.json"
        path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return errors, warnings


def run_partition(
    *,
    resolved_manifest: Path,
    bundle_id: str,
    mapping: Path,
    output_dir: Path,
) -> RunResult:
    if not SAFE_ID.fullmatch(bundle_id):
        raise PartitionInputError(f"Unsafe --bundle-id {bundle_id!r}")
    output = output_dir.expanduser().resolve()
    if output.exists():
        raise PartitionInputError(f"Refusing to reuse existing output directory: {output}")
    partitions = read_mapping(mapping, bundle_id)
    accession, assets = read_assets(resolved_manifest, bundle_id)
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        audit = SplitAudit(
            partitions=partitions,
            accession=accession,
            assets=assets,
            stage=stage,
            manifest_name=resolved_manifest.name,
        )
        errors, warnings = audit.execute()
        os.replace(stage, output)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return RunResult("PASS" if errors == 0 else "BLOCKED", output, errors, warnings)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolved-manifest", required=True, type=Path)
    parser.add_argument("--bundle-id", required=True)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = run_partition(
            resolved_manifest=args.resolved_manifest,
            bundle_id=args.bundle_id,
            mapping=args.mapping,
            output_dir=args.output_dir,
        )
    except (PartitionInputError, OSError, UnicodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "status": result.status,
                "output_dir": str(result.output_dir),
                "error_count": result.error_count,
                "warning_count": result.warning_count,
            },
            sort_keys=True,
        )
    )
    return 0 if result.status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
