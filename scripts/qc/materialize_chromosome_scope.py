#!/usr/bin/env python3
"""Materialize an audited chromosome-only genome and matching GFF3.

The source genome and annotation are opened read-only.  Retained records are
selected and renamed by an explicit, one-to-one genome/GFF/canonical sequence
ID table.  Genome and GFF IDs may differ.  In the explicit three-column schema,
an empty ``gff_seqid`` declares that a retained genome record intentionally has
no GFF3 feature rows; all other retained genome/GFF IDs must be paired
explicitly.  The program never infers correspondence from record order.
Unmapped genome records and their annotation rows are omitted from the
materialized assets but remain visible in the audit.

The implementation is deliberately fail-closed.  It validates the complete
FASTA and GFF3 before publishing anything, including sequence-ID uniqueness,
map completeness, feature coordinate bounds, ``ID``/``Parent`` closure, and
embedded-FASTA absence.  Repeated feature IDs are permitted when they describe
compatible multipart rows (for example, a CDS split across exons).  Outputs
are staged in a sibling temporary directory and installed atomically only
after every check passes.  An existing output directory is never replaced.

This program uses one process and starts no worker threads.
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
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO
from urllib.parse import unquote


SCRIPT_VERSION = "1.2.0"
READ_CHUNK_BYTES = 8 * 1024 * 1024
GZIP_MAGIC = b"\x1f\x8b"
GZIP_COMPRESSLEVEL = 6
EXPLICIT_MAP_COLUMNS = ("genome_seqid", "gff_seqid", "canonical_seqid")
LEGACY_MAP_COLUMNS = ("source_seqid", "canonical_seqid")
SAFE_PREFIX = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


class ChromosomeScopeError(RuntimeError):
    """Raised when a validated chromosome scope cannot be published."""


@dataclass(frozen=True)
class SequenceRecord:
    genome_seqid: str
    canonical_seqid: str | None
    length_bp: int
    retained: bool


@dataclass(frozen=True)
class FeatureIdSignature:
    gff_seqid: str
    feature_type: str
    strand: str
    parents: tuple[str, ...]
    retained: bool


@dataclass(frozen=True)
class ParentReference:
    parent_id: str
    gff_seqid: str
    retained: bool
    line_number: int
    child_id: str


@dataclass(frozen=True)
class MaterializationResult:
    output_dir: Path
    retained_sequences: int
    excluded_sequences: int
    retained_bp: int
    excluded_bp: int
    retained_features: int
    excluded_features: int


@dataclass(frozen=True)
class SequenceIdMaps:
    genome_to_canonical: dict[str, str]
    gff_to_canonical: dict[str, str]
    canonical_to_genome: dict[str, str]
    canonical_to_gff: dict[str, str]
    genome_seqids_declared_without_gff: frozenset[str]
    schema: str

    def __len__(self) -> int:
        return len(self.genome_to_canonical)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create an atomic chromosome-only genome/GFF3 pair from an explicit "
            "one-to-one sequence-ID map, with complete sequence and feature audits."
        )
    )
    parser.add_argument("--genome", required=True, type=Path, help="Genome FASTA (plain or gzip).")
    parser.add_argument("--gff", required=True, type=Path, help="GFF3 annotation (plain or gzip).")
    parser.add_argument(
        "--seqid-map",
        required=True,
        type=Path,
        help=(
            "Headered TSV with genome_seqid, gff_seqid, and canonical_seqid columns; "
            "in this explicit schema only, an empty gff_seqid declares a retained "
            "genome record with no GFF3 feature rows. "
            "the legacy source_seqid/canonical_seqid form is accepted only when genome "
            "and GFF IDs are identical."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="New output directory; an existing path is refused.",
    )
    parser.add_argument(
        "--prefix",
        default="chromosome_scope",
        help="Safe output filename prefix (default: chromosome_scope).",
    )
    return parser.parse_args(argv)


def absolute_path(path: Path) -> Path:
    expanded = path.expanduser()
    return expanded if expanded.is_absolute() else Path.cwd() / expanded


def path_entry_exists(path: Path) -> bool:
    """Return true for every existing directory entry, including dangling symlinks."""

    return os.path.lexists(path)


def sha256_file(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(READ_CHUNK_BYTES), b""):
                size += len(block)
                digest.update(block)
    except OSError as error:
        raise ChromosomeScopeError(f"Cannot checksum {path}: {error}") from error
    return size, digest.hexdigest()


def detect_compression(path: Path) -> str:
    try:
        with path.open("rb") as handle:
            return "gzip" if handle.read(2) == GZIP_MAGIC else "plain"
    except OSError as error:
        raise ChromosomeScopeError(f"Cannot inspect {path}: {error}") from error


@contextmanager
def open_input_text(path: Path) -> Iterator[TextIO]:
    compression = detect_compression(path)
    try:
        if compression == "gzip":
            with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
                yield handle
        else:
            with path.open("rt", encoding="utf-8", errors="strict", newline="") as handle:
                yield handle
    except (OSError, EOFError, UnicodeError, gzip.BadGzipFile) as error:
        raise ChromosomeScopeError(f"Cannot read {path}: {error}") from error


@contextmanager
def open_reproducible_gzip_text(path: Path) -> Iterator[TextIO]:
    try:
        with path.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                compresslevel=GZIP_COMPRESSLEVEL,
                fileobj=raw_handle,
                mtime=0,
            ) as gzip_handle:
                with io.TextIOWrapper(
                    gzip_handle, encoding="utf-8", errors="strict", newline=""
                ) as text_handle:
                    yield text_handle
    except (OSError, UnicodeError) as error:
        raise ChromosomeScopeError(f"Cannot write staged output {path}: {error}") from error


def validate_input_paths(
    genome_path: Path, gff_path: Path, map_path: Path, output_dir: Path
) -> None:
    inputs = {"genome": genome_path, "gff": gff_path, "seqid map": map_path}
    resolved_inputs: dict[str, Path] = {}
    for label, path in inputs.items():
        if not path.is_file():
            raise ChromosomeScopeError(f"The {label} input is missing or not a file: {path}")
        if path.stat().st_size == 0:
            raise ChromosomeScopeError(f"The {label} input is empty: {path}")
        resolved_inputs[label] = path.resolve(strict=True)
    if len(set(resolved_inputs.values())) != len(resolved_inputs):
        raise ChromosomeScopeError("Genome, GFF3, and sequence map must be different files")
    if path_entry_exists(output_dir):
        raise ChromosomeScopeError(f"Output path already exists; refusing overwrite: {output_dir}")
    output_resolved = output_dir.resolve(strict=False)
    if output_resolved in resolved_inputs.values():
        raise ChromosomeScopeError("Output directory must not replace an input file")


def read_seqid_map(path: Path) -> SequenceIdMaps:
    genome_to_canonical: dict[str, str] = {}
    gff_to_canonical: dict[str, str] = {}
    canonical_to_genome: dict[str, str] = {}
    canonical_to_gff: dict[str, str] = {}
    genome_seqids_declared_without_gff: set[str] = set()
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames == list(EXPLICIT_MAP_COLUMNS):
                schema = "explicit_genome_gff_canonical_v1"
            elif reader.fieldnames == list(LEGACY_MAP_COLUMNS):
                schema = "legacy_identical_source_ids_v1"
            else:
                raise ChromosomeScopeError(
                    f"{path}: map header must be exactly "
                    f"{'<TAB>'.join(EXPLICIT_MAP_COLUMNS)} (preferred) or "
                    f"{'<TAB>'.join(LEGACY_MAP_COLUMNS)} (identical genome/GFF IDs only); "
                    f"found {reader.fieldnames or []}"
                )
            for line_number, row in enumerate(reader, start=2):
                if row.get(None):
                    raise ChromosomeScopeError(
                        f"{path.name}:{line_number}: unexpected extra tab-separated field(s)"
                    )
                if schema == "legacy_identical_source_ids_v1":
                    genome_seqid = (row.get("source_seqid") or "").strip()
                    gff_seqid = genome_seqid
                else:
                    genome_seqid = (row.get("genome_seqid") or "").strip()
                    raw_gff_seqid = row.get("gff_seqid") or ""
                    gff_seqid = raw_gff_seqid.strip()
                canonical = (row.get("canonical_seqid") or "").strip()
                location = f"{path.name}:{line_number}"
                if not genome_seqid or any(
                    character.isspace() for character in genome_seqid
                ):
                    raise ChromosomeScopeError(
                        f"{location}: genome_seqid must be one non-empty FASTA token"
                    )
                if schema == "legacy_identical_source_ids_v1" and not gff_seqid:
                    raise ChromosomeScopeError(
                        f"{location}: gff_seqid must be one non-empty GFF3 sequence ID"
                    )
                if schema == "explicit_genome_gff_canonical_v1" and (
                    not gff_seqid and raw_gff_seqid != ""
                ):
                    raise ChromosomeScopeError(
                        f"{location}: gff_seqid must be empty or one non-whitespace "
                        "GFF3 sequence ID; whitespace-only values are not accepted"
                    )
                if gff_seqid and any(character.isspace() for character in gff_seqid):
                    raise ChromosomeScopeError(
                        f"{location}: gff_seqid must be empty or one non-whitespace "
                        "GFF3 sequence ID"
                    )
                if not canonical or any(character.isspace() for character in canonical):
                    raise ChromosomeScopeError(
                        f"{location}: canonical_seqid must be one non-empty FASTA token"
                    )
                if genome_seqid in genome_to_canonical:
                    raise ChromosomeScopeError(
                        f"{location}: duplicate genome_seqid {genome_seqid!r}"
                    )
                if gff_seqid and gff_seqid in gff_to_canonical:
                    raise ChromosomeScopeError(
                        f"{location}: duplicate gff_seqid {gff_seqid!r}"
                    )
                if not gff_seqid and genome_seqid in gff_to_canonical:
                    raise ChromosomeScopeError(
                        f"{location}: genome_seqid {genome_seqid!r} cannot be declared "
                        "without GFF features because that ID is already a mapped "
                        "gff_seqid"
                    )
                if gff_seqid and gff_seqid in genome_seqids_declared_without_gff:
                    raise ChromosomeScopeError(
                        f"{location}: gff_seqid {gff_seqid!r} collides with a "
                        "genome_seqid already declared without GFF features"
                    )
                if canonical in canonical_to_genome:
                    raise ChromosomeScopeError(
                        f"{location}: canonical_seqid {canonical!r} is already mapped from "
                        f"genome_seqid {canonical_to_genome[canonical]!r} and gff_seqid "
                        f"{canonical_to_gff.get(canonical, '')!r}"
                    )
                genome_to_canonical[genome_seqid] = canonical
                canonical_to_genome[canonical] = genome_seqid
                if gff_seqid:
                    gff_to_canonical[gff_seqid] = canonical
                    canonical_to_gff[canonical] = gff_seqid
                else:
                    genome_seqids_declared_without_gff.add(genome_seqid)
    except ChromosomeScopeError:
        raise
    except (OSError, UnicodeError, csv.Error) as error:
        raise ChromosomeScopeError(f"Cannot parse sequence map {path}: {error}") from error
    if not genome_to_canonical:
        raise ChromosomeScopeError(f"Sequence map contains no data rows: {path}")
    return SequenceIdMaps(
        genome_to_canonical=genome_to_canonical,
        gff_to_canonical=gff_to_canonical,
        canonical_to_genome=canonical_to_genome,
        canonical_to_gff=canonical_to_gff,
        genome_seqids_declared_without_gff=frozenset(
            genome_seqids_declared_without_gff
        ),
        schema=schema,
    )


def parse_fasta_header(raw_line: str, *, path: Path, line_number: int) -> tuple[str, str]:
    header = raw_line[1:].rstrip("\r\n")
    if not header:
        raise ChromosomeScopeError(f"{path.name}:{line_number}: empty FASTA header")
    match = re.match(r"(\S+)(.*)$", header)
    if match is None:
        raise ChromosomeScopeError(f"{path.name}:{line_number}: invalid FASTA header")
    return match.group(1), match.group(2)


def materialize_genome(
    genome_path: Path,
    mapping: dict[str, str],
    output_path: Path,
) -> tuple[list[SequenceRecord], dict[str, int]]:
    sequence_records: list[SequenceRecord] = []
    lengths: dict[str, int] = {}
    observed: set[str] = set()
    current_source: str | None = None
    current_canonical: str | None = None
    current_length = 0
    current_retained = False
    record_count = 0

    def finish_record() -> None:
        nonlocal current_source, current_canonical, current_length, current_retained
        if current_source is None:
            return
        if current_length == 0:
            raise ChromosomeScopeError(
                f"{genome_path.name}: FASTA record {current_source!r} has no sequence"
            )
        sequence_records.append(
            SequenceRecord(
                genome_seqid=current_source,
                canonical_seqid=current_canonical,
                length_bp=current_length,
                retained=current_retained,
            )
        )
        lengths[current_source] = current_length

    with open_input_text(genome_path) as input_handle, open_reproducible_gzip_text(
        output_path
    ) as output_handle:
        for line_number, line in enumerate(input_handle, start=1):
            if line.startswith(">"):
                finish_record()
                source, description = parse_fasta_header(
                    line, path=genome_path, line_number=line_number
                )
                if source in observed:
                    raise ChromosomeScopeError(
                        f"{genome_path.name}:{line_number}: duplicate FASTA seqid {source!r}"
                    )
                observed.add(source)
                record_count += 1
                current_source = source
                current_canonical = mapping.get(source)
                current_retained = current_canonical is not None
                current_length = 0
                if current_retained:
                    output_handle.write(f">{current_canonical}{description}\n")
                continue
            sequence = "".join(line.split())
            if not sequence:
                continue
            if current_source is None:
                raise ChromosomeScopeError(
                    f"{genome_path.name}:{line_number}: sequence appears before the first header"
                )
            current_length += len(sequence)
            if current_retained:
                output_handle.write(sequence + "\n")
        finish_record()
    if record_count == 0:
        raise ChromosomeScopeError(f"Genome FASTA contains no records: {genome_path}")
    missing = sorted(set(mapping) - observed)
    if missing:
        preview = ", ".join(missing[:10])
        suffix = "" if len(missing) <= 10 else f" (and {len(missing) - 10} more)"
        raise ChromosomeScopeError(
            f"Sequence map contains {len(missing)} genome_seqid value(s) absent from the "
            f"genome: {preview}{suffix}"
        )
    return sequence_records, lengths


def parse_gff_attributes(
    raw_attributes: str, *, path: Path, line_number: int
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    if raw_attributes == ".":
        return (), ()
    ids: list[str] = []
    parents: list[str] = []
    for field in raw_attributes.split(";"):
        if not field:
            continue
        if "=" not in field:
            raise ChromosomeScopeError(
                f"{path.name}:{line_number}: malformed GFF3 attribute {field!r}"
            )
        key, raw_value = field.split("=", 1)
        if not key or raw_value == "":
            raise ChromosomeScopeError(
                f"{path.name}:{line_number}: empty GFF3 attribute key or value"
            )
        if key == "ID":
            value = unquote(raw_value)
            if not value:
                raise ChromosomeScopeError(
                    f"{path.name}:{line_number}: empty decoded ID attribute"
                )
            ids.append(value)
        elif key == "Parent":
            for raw_parent in raw_value.split(","):
                parent = unquote(raw_parent)
                if not parent:
                    raise ChromosomeScopeError(
                        f"{path.name}:{line_number}: empty Parent identifier"
                    )
                parents.append(parent)
    if len(ids) > 1:
        raise ChromosomeScopeError(
            f"{path.name}:{line_number}: a feature row may contain only one ID attribute"
        )
    if len(parents) != len(set(parents)):
        raise ChromosomeScopeError(
            f"{path.name}:{line_number}: duplicate Parent identifier in one feature row"
        )
    return tuple(ids), tuple(sorted(parents))


def materialize_gff(
    gff_path: Path,
    seqid_maps: SequenceIdMaps,
    lengths: dict[str, int],
    output_path: Path,
) -> tuple[
    Counter[tuple[bool, str, str]],
    int,
    int,
    int,
]:
    feature_counts: Counter[tuple[bool, str, str]] = Counter()
    feature_ids: dict[str, FeatureIdSignature] = {}
    parent_references: list[ParentReference] = []
    sequence_regions: set[str] = set()
    observed_feature_seqids: set[str] = set()
    retained_features = 0
    excluded_features = 0
    feature_rows = 0

    with open_input_text(gff_path) as input_handle, open_reproducible_gzip_text(
        output_path
    ) as output_handle:
        for line_number, line in enumerate(input_handle, start=1):
            stripped = line.rstrip("\r\n")
            if stripped == "##FASTA":
                raise ChromosomeScopeError(
                    f"{gff_path.name}:{line_number}: embedded GFF3 FASTA is not accepted; "
                    "provide a feature-only GFF3"
                )
            if stripped.startswith("##sequence-region"):
                fields = stripped.split()
                if len(fields) != 4:
                    raise ChromosomeScopeError(
                        f"{gff_path.name}:{line_number}: malformed ##sequence-region directive"
                    )
                _, gff_seqid, raw_start, raw_end = fields
                if gff_seqid in sequence_regions:
                    raise ChromosomeScopeError(
                        f"{gff_path.name}:{line_number}: duplicate ##sequence-region for "
                        f"{gff_seqid!r}"
                    )
                sequence_regions.add(gff_seqid)
                canonical = seqid_maps.gff_to_canonical.get(gff_seqid)
                genome_seqid = (
                    seqid_maps.canonical_to_genome[canonical]
                    if canonical is not None
                    else gff_seqid
                )
                if genome_seqid not in lengths:
                    raise ChromosomeScopeError(
                        f"{gff_path.name}:{line_number}: ##sequence-region seqid "
                        f"{gff_seqid!r} has no genome sequence with the same ID and no "
                        "explicit retained genome_seqid association"
                    )
                try:
                    start = int(raw_start)
                    end = int(raw_end)
                except ValueError as error:
                    raise ChromosomeScopeError(
                        f"{gff_path.name}:{line_number}: non-integer ##sequence-region bounds"
                    ) from error
                if start < 1 or end < start or end > lengths[genome_seqid]:
                    raise ChromosomeScopeError(
                        f"{gff_path.name}:{line_number}: ##sequence-region {start}-{end} "
                        f"is outside genome_seqid {genome_seqid!r} length "
                        f"{lengths[genome_seqid]}"
                    )
                if canonical is not None:
                    output_handle.write(
                        f"##sequence-region {canonical} {start} {end}\n"
                    )
                continue
            if not stripped or stripped.startswith("#"):
                output_handle.write(line if line.endswith(("\n", "\r")) else line + "\n")
                continue
            columns = stripped.split("\t")
            if len(columns) != 9:
                raise ChromosomeScopeError(
                    f"{gff_path.name}:{line_number}: expected 9 tab-separated GFF3 columns, "
                    f"found {len(columns)}"
                )
            gff_seqid = columns[0]
            if gff_seqid in seqid_maps.genome_seqids_declared_without_gff:
                raise ChromosomeScopeError(
                    f"{gff_path.name}:{line_number}: GFF3 feature uses seqid "
                    f"{gff_seqid!r}, but that retained genome_seqid was explicitly "
                    "declared without GFF features"
                )
            canonical = seqid_maps.gff_to_canonical.get(gff_seqid)
            genome_seqid = (
                seqid_maps.canonical_to_genome[canonical]
                if canonical is not None
                else gff_seqid
            )
            if genome_seqid not in lengths:
                raise ChromosomeScopeError(
                    f"{gff_path.name}:{line_number}: GFF3 seqid {gff_seqid!r} has no "
                    "genome sequence with the same ID and no explicit retained "
                    "genome_seqid association"
                )
            observed_feature_seqids.add(gff_seqid)
            try:
                start = int(columns[3])
                end = int(columns[4])
            except ValueError as error:
                raise ChromosomeScopeError(
                    f"{gff_path.name}:{line_number}: feature start/end must be integers"
                ) from error
            if start < 1 or end < start or end > lengths[genome_seqid]:
                raise ChromosomeScopeError(
                    f"{gff_path.name}:{line_number}: feature {start}-{end} is outside "
                    f"genome_seqid {genome_seqid!r} length {lengths[genome_seqid]}"
                )
            retained = canonical is not None
            feature_type = columns[2]
            if not feature_type or feature_type == ".":
                raise ChromosomeScopeError(
                    f"{gff_path.name}:{line_number}: feature type must be non-empty"
                )
            ids, parents = parse_gff_attributes(
                columns[8], path=gff_path, line_number=line_number
            )
            child_id = ids[0] if ids else ""
            if child_id:
                signature = FeatureIdSignature(
                    gff_seqid=gff_seqid,
                    feature_type=feature_type,
                    strand=columns[6],
                    parents=parents,
                    retained=retained,
                )
                previous = feature_ids.get(child_id)
                if previous is not None and previous != signature:
                    raise ChromosomeScopeError(
                        f"{gff_path.name}:{line_number}: repeated ID {child_id!r} has "
                        "incompatible seqid, type, strand, Parent, or retained scope; only "
                        "compatible multipart feature rows may repeat an ID"
                    )
                feature_ids[child_id] = signature
            for parent_id in parents:
                parent_references.append(
                    ParentReference(
                        parent_id=parent_id,
                        gff_seqid=gff_seqid,
                        retained=retained,
                        line_number=line_number,
                        child_id=child_id,
                    )
                )
            feature_counts[(retained, gff_seqid, feature_type)] += 1
            feature_rows += 1
            if retained:
                retained_features += 1
                columns[0] = canonical
                output_handle.write("\t".join(columns) + "\n")
            else:
                excluded_features += 1

    if feature_rows == 0:
        raise ChromosomeScopeError(f"GFF3 contains no feature rows: {gff_path}")
    missing_mapped_gff_seqids = sorted(
        set(seqid_maps.gff_to_canonical) - observed_feature_seqids
    )
    if missing_mapped_gff_seqids:
        preview = ", ".join(missing_mapped_gff_seqids[:10])
        suffix = (
            ""
            if len(missing_mapped_gff_seqids) <= 10
            else f" (and {len(missing_mapped_gff_seqids) - 10} more)"
        )
        raise ChromosomeScopeError(
            f"Sequence map contains {len(missing_mapped_gff_seqids)} gff_seqid value(s) "
            f"without any GFF3 feature row: {preview}{suffix}"
        )
    for reference in parent_references:
        parent = feature_ids.get(reference.parent_id)
        location = f"{gff_path.name}:{reference.line_number}"
        if parent is None:
            child = f" for child {reference.child_id!r}" if reference.child_id else ""
            raise ChromosomeScopeError(
                f"{location}: Parent {reference.parent_id!r}{child} has no matching ID"
            )
        if parent.retained != reference.retained:
            raise ChromosomeScopeError(
                f"{location}: Parent {reference.parent_id!r} crosses retained/excluded scope"
            )
        if parent.gff_seqid != reference.gff_seqid:
            raise ChromosomeScopeError(
                f"{location}: Parent {reference.parent_id!r} is on a different GFF3 seqid"
            )
    return feature_counts, retained_features, excluded_features, feature_rows


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    try:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(columns),
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            writer.writerows(rows)
    except (OSError, UnicodeError, csv.Error) as error:
        raise ChromosomeScopeError(f"Cannot write staged table {path}: {error}") from error


def build_sequence_audit(
    records: list[SequenceRecord],
    feature_counts: Counter[tuple[bool, str, str]],
    seqid_maps: SequenceIdMaps,
) -> list[dict[str, object]]:
    totals: Counter[tuple[bool, str]] = Counter()
    for (retained, gff_seqid, _feature_type), count in feature_counts.items():
        totals[(retained, gff_seqid)] += count
    rows: list[dict[str, object]] = []
    for record in records:
        if record.retained:
            assert record.canonical_seqid is not None
            gff_seqid = seqid_maps.canonical_to_gff.get(record.canonical_seqid, "")
        elif record.genome_seqid not in seqid_maps.gff_to_canonical and totals[
            (False, record.genome_seqid)
        ]:
            gff_seqid = record.genome_seqid
        else:
            gff_seqid = ""
        rows.append(
            {
                "scope": "retained" if record.retained else "excluded",
                "genome_seqid": record.genome_seqid,
                "gff_seqid": gff_seqid,
                "canonical_seqid": record.canonical_seqid or "",
                "length_bp": record.length_bp,
                "gff_feature_count": totals[(record.retained, gff_seqid)] if gff_seqid else 0,
            }
        )
    return rows


def build_feature_audit(
    feature_counts: Counter[tuple[bool, str, str]],
    seqid_maps: SequenceIdMaps,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (retained, gff_seqid, feature_type), count in sorted(
        feature_counts.items(),
        key=lambda item: (
            0 if item[0][0] else 1,
            seqid_maps.gff_to_canonical.get(item[0][1], item[0][1]),
            item[0][2],
        ),
    ):
        canonical_seqid = seqid_maps.gff_to_canonical.get(gff_seqid, "")
        genome_seqid = (
            seqid_maps.canonical_to_genome[canonical_seqid]
            if canonical_seqid
            else gff_seqid
        )
        rows.append(
            {
                "scope": "retained" if retained else "excluded",
                "genome_seqid": genome_seqid,
                "gff_seqid": gff_seqid,
                "canonical_seqid": canonical_seqid,
                "feature_type": feature_type,
                "feature_count": count,
            }
        )
    return rows


def write_json(path: Path, payload: dict[str, object]) -> None:
    try:
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except (OSError, UnicodeError) as error:
        raise ChromosomeScopeError(f"Cannot write staged JSON {path}: {error}") from error


def materialize_chromosome_scope(
    *,
    genome_path: Path,
    gff_path: Path,
    map_path: Path,
    output_dir: Path,
    prefix: str = "chromosome_scope",
) -> MaterializationResult:
    if not SAFE_PREFIX.fullmatch(prefix):
        raise ChromosomeScopeError(
            "--prefix must start with an alphanumeric character and contain only "
            "letters, numbers, periods, underscores, or hyphens (maximum 128 characters)"
        )
    genome_path = absolute_path(genome_path)
    gff_path = absolute_path(gff_path)
    map_path = absolute_path(map_path)
    output_dir = absolute_path(output_dir)
    validate_input_paths(genome_path, gff_path, map_path, output_dir)

    input_files: dict[str, dict[str, object]] = {}
    for role, path in (("genome", genome_path), ("gff", gff_path), ("seqid_map", map_path)):
        size, digest = sha256_file(path)
        input_files[role] = {
            "file_name": path.name,
            "size_bytes": size,
            "sha256": digest,
            "compression": detect_compression(path) if role != "seqid_map" else "plain",
        }
    seqid_maps = read_seqid_map(map_path)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staged_root = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", suffix=".tmp", dir=output_dir.parent)
    )
    try:
        audit_dir = staged_root / "audit"
        audit_dir.mkdir()
        genome_output = staged_root / f"{prefix}.genome.fa.gz"
        gff_output = staged_root / f"{prefix}.annotation.gff3.gz"
        sequence_audit_path = audit_dir / "sequence_scope.tsv"
        feature_audit_path = audit_dir / "feature_counts.tsv"
        validation_path = audit_dir / "validation.json"
        checksums_path = staged_root / "checksums.tsv"

        records, lengths = materialize_genome(
            genome_path, seqid_maps.genome_to_canonical, genome_output
        )
        feature_counts, retained_features, excluded_features, feature_rows = materialize_gff(
            gff_path, seqid_maps, lengths, gff_output
        )
        for role, path in (("genome", genome_path), ("gff", gff_path), ("seqid_map", map_path)):
            final_size, final_digest = sha256_file(path)
            if (
                final_size != input_files[role]["size_bytes"]
                or final_digest != input_files[role]["sha256"]
            ):
                raise ChromosomeScopeError(
                    f"The {role} input changed while it was being read; staged outputs were discarded"
                )
        sequence_rows = build_sequence_audit(records, feature_counts, seqid_maps)
        feature_rows_audit = build_feature_audit(feature_counts, seqid_maps)
        write_tsv(
            sequence_audit_path,
            (
                "scope",
                "genome_seqid",
                "gff_seqid",
                "canonical_seqid",
                "length_bp",
                "gff_feature_count",
            ),
            sequence_rows,
        )
        write_tsv(
            feature_audit_path,
            (
                "scope",
                "genome_seqid",
                "gff_seqid",
                "canonical_seqid",
                "feature_type",
                "feature_count",
            ),
            feature_rows_audit,
        )

        retained_records = [record for record in records if record.retained]
        excluded_records = [record for record in records if not record.retained]
        output_payloads: dict[str, dict[str, object]] = {}
        for role, path in (
            ("chromosome_genome", genome_output),
            ("chromosome_gff", gff_output),
            ("sequence_audit", sequence_audit_path),
            ("feature_audit", feature_audit_path),
        ):
            size, digest = sha256_file(path)
            output_payloads[role] = {
                "relative_path": path.relative_to(staged_root).as_posix(),
                "size_bytes": size,
                "sha256": digest,
            }

        validation_payload: dict[str, object] = {
            "schema_version": "1.2",
            "script_version": SCRIPT_VERSION,
            "status": "PASS",
            "policy": {
                "selection": (
                    "retain every mapped genome_seqid; retain GFF3 feature rows only for "
                    "non-empty gff_seqid values paired in the map"
                ),
                "renaming": (
                    "replace the first FASTA token and mapped GFF3 column 1 with canonical_seqid"
                ),
                "seqid_map_schema": seqid_maps.schema,
                "empty_gff_seqid": (
                    "allowed only in the explicit three-column schema and declares a "
                    "retained genome sequence with zero GFF3 feature rows"
                ),
                "seqid_correspondence_inferred_by_order": False,
                "raw_inputs_modified": False,
                "overwrite_existing_output": False,
                "processes": 1,
                "worker_threads_started": 0,
            },
            "checks": {
                "map_is_nonempty_one_to_one": "PASS",
                "input_files_unchanged_during_run": "PASS",
                "mapped_genome_seqids_exist_once_in_genome": "PASS",
                "mapped_gff_seqids_have_feature_rows": "PASS",
                "nonempty_mapped_gff_seqids_have_feature_rows": "PASS",
                "retained_genome_seqids_declared_without_gff_have_zero_feature_rows": "PASS",
                "genome_seqids_are_unique": "PASS",
                "gff_seqids_have_explicit_or_identical_genome_associations": "PASS",
                "feature_coordinates_within_sequence_bounds": "PASS",
                "id_parent_closure_and_scope": "PASS",
                "embedded_gff_fasta_absent": "PASS",
            },
            "counts": {
                "map_rows": len(seqid_maps),
                "nonempty_mapped_gff_seqids": len(seqid_maps.gff_to_canonical),
                "retained_genome_sequences_declared_without_gff": len(
                    seqid_maps.genome_seqids_declared_without_gff
                ),
                "input_genome_sequences": len(records),
                "retained_genome_sequences": len(retained_records),
                "excluded_genome_sequences": len(excluded_records),
                "retained_bp": sum(record.length_bp for record in retained_records),
                "excluded_bp": sum(record.length_bp for record in excluded_records),
                "input_gff_features": feature_rows,
                "retained_gff_features": retained_features,
                "excluded_gff_features": excluded_features,
            },
            "input_files": input_files,
            "outputs": output_payloads,
        }
        write_json(validation_path, validation_payload)

        checksum_rows: list[dict[str, object]] = []
        for role in ("genome", "gff", "seqid_map"):
            metadata = input_files[role]
            checksum_rows.append(
                {
                    "scope": "input",
                    "role": role,
                    "path": metadata["file_name"],
                    "size_bytes": metadata["size_bytes"],
                    "sha256": metadata["sha256"],
                }
            )
        for role, path in (
            ("chromosome_genome", genome_output),
            ("chromosome_gff", gff_output),
            ("sequence_audit", sequence_audit_path),
            ("feature_audit", feature_audit_path),
            ("validation", validation_path),
        ):
            size, digest = sha256_file(path)
            checksum_rows.append(
                {
                    "scope": "output",
                    "role": role,
                    "path": path.relative_to(staged_root).as_posix(),
                    "size_bytes": size,
                    "sha256": digest,
                }
            )
        write_tsv(
            checksums_path,
            ("scope", "role", "path", "size_bytes", "sha256"),
            checksum_rows,
        )

        if path_entry_exists(output_dir):
            raise ChromosomeScopeError(
                f"Output path appeared during processing; refusing overwrite: {output_dir}"
            )
        os.rename(staged_root, output_dir)
        return MaterializationResult(
            output_dir=output_dir,
            retained_sequences=len(retained_records),
            excluded_sequences=len(excluded_records),
            retained_bp=sum(record.length_bp for record in retained_records),
            excluded_bp=sum(record.length_bp for record in excluded_records),
            retained_features=retained_features,
            excluded_features=excluded_features,
        )
    except Exception:
        shutil.rmtree(staged_root, ignore_errors=True)
        raise


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = materialize_chromosome_scope(
            genome_path=args.genome,
            gff_path=args.gff,
            map_path=args.seqid_map,
            output_dir=args.output_dir,
            prefix=args.prefix,
        )
    except ChromosomeScopeError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(
        "PASS: "
        f"{result.retained_sequences} retained sequences ({result.retained_bp} bp), "
        f"{result.excluded_sequences} excluded sequences ({result.excluded_bp} bp), "
        f"{result.retained_features} retained features, "
        f"{result.excluded_features} excluded features -> {result.output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
