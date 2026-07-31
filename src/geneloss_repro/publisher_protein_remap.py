"""Fail-closed publisher-protein remapping for a selected primary transcript set.

Some publisher protein FASTA files use protein accessions as their first-token
identifiers even though the matched GFF3 uses transcript IDs.  This module
validates the complete GFF3-to-publisher mapping, subsets the publisher FASTA
to the already selected primary transcripts, and changes only the output
record identifiers.  Protein sequence characters are never normalized.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO
from urllib.parse import unquote


WORKFLOW_VERSION = "1.2.0"
GZIP_MAGIC = b"\x1f\x8b"
SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
SAFE_SCHEMA_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_.:-]*$")
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYXBJOUZ*")
PUBLISHER_HEADER_MODES = frozenset({"metadata", "first_token"})
TRANSCRIPT_ACCESSION_SOURCES = frozenset({"attribute", "transcript_id"})
PROTEIN_ACCESSION_SOURCES = frozenset({"attribute", "cds_parent"})


class PublisherProteinRemapError(RuntimeError):
    """Raised when an exact publisher-primary remap cannot be published."""


@dataclass(frozen=True)
class PublisherProteinRemapResult:
    """Counts from one atomically published publisher-primary remap."""

    output_dir: Path
    source_publisher_record_count: int
    selected_primary_record_count: int
    excluded_nonprimary_record_count: int
    output_protein_path: Path


@dataclass(frozen=True)
class _FastaRecord:
    identifier: str
    metadata: dict[str, str]
    sequence: str

    @property
    def sequence_sha256(self) -> str:
        return hashlib.sha256(self.sequence.encode("ascii")).hexdigest()


@dataclass(frozen=True)
class _GffProteinMap:
    transcript_accessions: dict[str, str]
    protein_by_transcript: dict[str, str]
    transcript_by_protein: dict[str, str]
    transcript_count: int
    gene_count: int
    source_transcript_row_count: int
    noncoding_model_count: int
    graph_mode: str


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == GZIP_MAGIC
    except OSError as error:
        raise PublisherProteinRemapError(
            f"Cannot inspect input {path.name}: {error}"
        ) from error


def _open_text(path: Path) -> TextIO:
    try:
        if _is_gzip(path):
            return gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="")
        return path.open("rt", encoding="utf-8", errors="strict", newline="")
    except OSError as error:
        raise PublisherProteinRemapError(
            f"Cannot open input {path.name}: {error}"
        ) from error


def _checksum(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                size += len(block)
                digest.update(block)
    except OSError as error:
        raise PublisherProteinRemapError(
            f"Cannot checksum input {path.name}: {error}"
        ) from error
    return size, digest.hexdigest()


def _input_signature(path: Path) -> tuple[int, int, int, int]:
    try:
        status = path.stat()
    except OSError as error:
        raise PublisherProteinRemapError(
            f"Cannot stat input {path.name}: {error}"
        ) from error
    return status.st_dev, status.st_ino, status.st_size, status.st_mtime_ns


def _validate_identifier(identifier: str, context: str) -> None:
    if not identifier or any(
        character.isspace() or ord(character) < 32 for character in identifier
    ):
        raise PublisherProteinRemapError(
            f"{context} must be a non-empty identifier without whitespace"
        )


def _validate_protein_sequence(sequence: str, identifier: str, role: str) -> None:
    if not sequence:
        raise PublisherProteinRemapError(
            f"{role} FASTA record {identifier!r} has an empty sequence"
        )
    invalid = sorted(set(sequence) - PROTEIN_ALPHABET)
    if invalid:
        raise PublisherProteinRemapError(
            f"{role} FASTA record {identifier!r} contains unsupported symbols: "
            + ",".join(invalid[:10])
        )
    stop_positions = [index for index, residue in enumerate(sequence) if residue == "*"]
    if stop_positions and stop_positions != [len(sequence) - 1]:
        raise PublisherProteinRemapError(
            f"{role} FASTA record {identifier!r} contains an internal or repeated stop codon"
        )


def _parse_header_metadata(
    header: str,
    *,
    required_keys: tuple[str, ...],
    role: str,
    line_number: int,
) -> tuple[str, dict[str, str]]:
    fields = header.split("\t")
    tokens = fields[0].split()
    if not tokens:
        raise PublisherProteinRemapError(
            f"{role} FASTA line {line_number} has an empty header"
        )
    identifier = tokens[0]
    _validate_identifier(identifier, f"{role} FASTA line {line_number} identifier")
    metadata: dict[str, str] = {}
    for field in fields[1:]:
        if not field or "=" not in field:
            raise PublisherProteinRemapError(
                f"{role} FASTA line {line_number} has malformed tab-delimited metadata"
            )
        key, value = field.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key or key in metadata:
            raise PublisherProteinRemapError(
                f"{role} FASTA line {line_number} has an empty or duplicate metadata key"
            )
        metadata[key] = value
    for key in required_keys:
        value = metadata.get(key, "")
        if not value:
            raise PublisherProteinRemapError(
                f"{role} FASTA line {line_number} is missing required header field {key!r}"
            )
        _validate_identifier(value, f"{role} FASTA header field {key!r}")
    return identifier, metadata


def _read_protein_fasta(
    path: Path,
    role: str,
    *,
    required_header_keys: tuple[str, ...] = (),
) -> dict[str, _FastaRecord]:
    records: dict[str, _FastaRecord] = {}
    identifier: str | None = None
    metadata: dict[str, str] = {}
    parts: list[str] = []

    def finish() -> None:
        nonlocal identifier, metadata, parts
        if identifier is None:
            return
        sequence = "".join(parts)
        _validate_protein_sequence(sequence, identifier, role)
        if identifier in records:
            raise PublisherProteinRemapError(
                f"{role} protein FASTA repeats identifier {identifier!r}"
            )
        records[identifier] = _FastaRecord(identifier, dict(metadata), sequence)

    try:
        with _open_text(path) as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith(">"):
                    finish()
                    identifier, metadata = _parse_header_metadata(
                        line[1:],
                        required_keys=required_header_keys,
                        role=role,
                        line_number=line_number,
                    )
                    if identifier in records:
                        raise PublisherProteinRemapError(
                            f"{role} protein FASTA repeats identifier {identifier!r}"
                        )
                    parts = []
                elif identifier is None:
                    raise PublisherProteinRemapError(
                        f"{role} protein FASTA line {line_number} has sequence before a header"
                    )
                else:
                    if any(character.isspace() for character in line):
                        raise PublisherProteinRemapError(
                            f"{role} FASTA record {identifier!r} contains embedded whitespace"
                        )
                    parts.append(line)
        finish()
    except (UnicodeError, EOFError, gzip.BadGzipFile) as error:
        raise PublisherProteinRemapError(
            f"Cannot read {role} protein FASTA {path.name}: {error}"
        ) from error
    if not records:
        raise PublisherProteinRemapError(
            f"{role} protein FASTA {path.name} contains no records"
        )
    return records


def _parse_gff_attributes(raw: str, source: str, line_number: int) -> dict[str, str]:
    attributes: dict[str, str] = {}
    if raw == ".":
        return attributes
    for field in raw.split(";"):
        if not field:
            continue
        if "=" not in field:
            raise PublisherProteinRemapError(
                f"{source}:{line_number}: expected GFF3 key=value attributes"
            )
        encoded_key, encoded_value = field.split("=", 1)
        key = unquote(encoded_key.strip())
        value = unquote(encoded_value.strip())
        if not key or key in attributes:
            raise PublisherProteinRemapError(
                f"{source}:{line_number}: empty or duplicate GFF3 attribute key"
            )
        attributes[key] = value
    return attributes


def _parse_gff_protein_map(
    path: Path,
    *,
    transcript_features: frozenset[str],
    gene_features: frozenset[str],
    transcript_id_attribute: str,
    transcript_accession_attribute: str,
    transcript_accession_source: str,
    cds_parent_attribute: str,
    protein_accession_attribute: str,
    protein_accession_source: str,
    gene_as_transcript: bool,
) -> _GffProteinMap:
    transcript_rows: dict[str, tuple[dict[str, str], int]] = {}
    gene_rows: dict[str, tuple[str, dict[str, str], int]] = {}
    cds_rows: list[tuple[tuple[str, ...], dict[str, str], int]] = []
    gene_row_count = 0
    saw_feature = False
    try:
        with _open_text(path) as handle:
            for line_number, raw in enumerate(handle, start=1):
                if raw.startswith("##FASTA"):
                    raise PublisherProteinRemapError(
                        f"{path.name}:{line_number}: embedded FASTA is not accepted"
                    )
                if not raw.strip() or raw.startswith("#"):
                    continue
                saw_feature = True
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 9:
                    raise PublisherProteinRemapError(
                        f"{path.name}:{line_number}: expected 9 GFF3 columns, found {len(fields)}"
                    )
                feature_type = fields[2]
                if feature_type in gene_features and not gene_as_transcript:
                    # Gene rows were intentionally outside the original normal-
                    # transcript remap contract. Count them without introducing
                    # new attribute/ID validation into that established mode.
                    gene_row_count += 1
                    continue
                if (
                    feature_type not in transcript_features
                    and feature_type not in gene_features
                    and feature_type != "CDS"
                ):
                    continue
                attributes = _parse_gff_attributes(fields[8], path.name, line_number)
                if feature_type in transcript_features:
                    transcript_id = attributes.get(transcript_id_attribute, "")
                    _validate_identifier(
                        transcript_id,
                        f"{path.name}:{line_number} transcript {transcript_id_attribute}",
                    )
                    if transcript_id in transcript_rows:
                        raise PublisherProteinRemapError(
                            f"{path.name}:{line_number}: duplicate transcript ID {transcript_id!r}"
                        )
                    transcript_rows[transcript_id] = (attributes, line_number)
                elif feature_type in gene_features:
                    gene_row_count += 1
                    gene_id = attributes.get(transcript_id_attribute, "")
                    _validate_identifier(
                        gene_id,
                        f"{path.name}:{line_number} gene {transcript_id_attribute}",
                    )
                    if gene_id in gene_rows:
                        raise PublisherProteinRemapError(
                            f"{path.name}:{line_number}: duplicate gene-level ID {gene_id!r}"
                        )
                    gene_rows[gene_id] = (feature_type, attributes, line_number)
                else:
                    parent_value = attributes.get(cds_parent_attribute, "")
                    parents = [item.strip() for item in parent_value.split(",") if item.strip()]
                    if not parents:
                        raise PublisherProteinRemapError(
                            f"{path.name}:{line_number}: CDS lacks {cds_parent_attribute!r}"
                        )
                    for parent in parents:
                        _validate_identifier(parent, f"{path.name}:{line_number} CDS parent")
                    cds_rows.append((tuple(parents), attributes, line_number))
    except (UnicodeError, EOFError, gzip.BadGzipFile) as error:
        raise PublisherProteinRemapError(
            f"Cannot read GFF3 {path.name}: {error}"
        ) from error
    if not saw_feature:
        raise PublisherProteinRemapError(
            f"{path.name}: no GFF3 feature rows were found"
        )

    def one_attribute_value(
        attributes: dict[str, str], attribute: str, context: str
    ) -> str:
        raw_value = attributes.get(attribute, "")
        values = [value.strip() for value in raw_value.split(",") if value.strip()]
        if len(values) != 1:
            raise PublisherProteinRemapError(
                f"{context} must declare exactly one {attribute!r} accession"
            )
        _validate_identifier(values[0], f"{context} {attribute}")
        return values[0]

    def model_accession(
        model_id: str, attributes: dict[str, str], line_number: int, role: str
    ) -> str:
        if transcript_accession_source == "transcript_id":
            return model_id
        return one_attribute_value(
            attributes,
            transcript_accession_attribute,
            f"{path.name}:{line_number} {role}",
        )

    def cds_protein_accession(
        parents: tuple[str, ...], attributes: dict[str, str], line_number: int
    ) -> str:
        if protein_accession_source == "cds_parent":
            if len(parents) != 1:
                raise PublisherProteinRemapError(
                    f"{path.name}:{line_number}: CDS-parent protein accession mode "
                    "requires exactly one CDS parent"
                )
            return parents[0]
        return one_attribute_value(
            attributes,
            protein_accession_attribute,
            f"{path.name}:{line_number} CDS",
        )

    transcript_accessions: dict[str, str] = {}
    accession_to_transcript: dict[str, str] = {}
    proteins_by_transcript: dict[str, set[str]] = defaultdict(set)

    if gene_as_transcript:
        if transcript_rows:
            raise PublisherProteinRemapError(
                f"{path.name}: --gene-as-transcript requires zero accepted transcript rows; "
                f"found {len(transcript_rows)}"
            )
        if not gene_rows:
            raise PublisherProteinRemapError(
                f"{path.name}: --gene-as-transcript requires declared gene-level rows"
            )
        for gene_id, (_, attributes, line_number) in gene_rows.items():
            accession = model_accession(gene_id, attributes, line_number, "gene")
            if accession in accession_to_transcript:
                raise PublisherProteinRemapError(
                    f"{path.name}:{line_number}: transcript accession {accession!r} "
                    "maps to more than one synthesized gene transcript"
                )
            transcript_accessions[gene_id] = accession
            accession_to_transcript[accession] = gene_id
        for parents, attributes, line_number in cds_rows:
            if len(parents) != 1:
                raise PublisherProteinRemapError(
                    f"{path.name}:{line_number}: gene-as-transcript CDS must have exactly "
                    "one declared gene parent"
                )
            parent = parents[0]
            if parent not in gene_rows:
                raise PublisherProteinRemapError(
                    f"{path.name}:{line_number}: gene-as-transcript CDS parent {parent!r} "
                    "is not a declared gene-level ID"
                )
            if gene_rows[parent][0].lower() == "pseudogene":
                raise PublisherProteinRemapError(
                    f"{path.name}:{line_number}: CDS is attached to declared pseudogene "
                    f"{parent!r}"
                )
            proteins_by_transcript[parent].add(
                cds_protein_accession(parents, attributes, line_number)
            )
        graph_mode = "gene_as_transcript"
        model_ids = set(gene_rows)
    else:
        if not transcript_rows:
            raise PublisherProteinRemapError(
                f"{path.name}: no accepted transcript rows were found; the explicit "
                "--gene-as-transcript option is required for a gene-to-CDS graph"
            )
        for transcript_id, (attributes, line_number) in transcript_rows.items():
            accession = model_accession(
                transcript_id, attributes, line_number, "transcript"
            )
            if accession in accession_to_transcript:
                raise PublisherProteinRemapError(
                    f"{path.name}:{line_number}: transcript accession {accession!r} "
                    "maps to more than one transcript"
                )
            transcript_accessions[transcript_id] = accession
            accession_to_transcript[accession] = transcript_id
        for parents, attributes, line_number in cds_rows:
            protein_accession = cds_protein_accession(
                parents, attributes, line_number
            )
            for parent in parents:
                proteins_by_transcript[parent].add(protein_accession)
        graph_mode = "declared_transcripts"
        model_ids = set(transcript_rows)

    if not proteins_by_transcript:
        raise PublisherProteinRemapError(
            f"{path.name}: no usable transcript-to-protein annotation was found"
        )

    orphan_parents = sorted(set(proteins_by_transcript) - model_ids)
    if orphan_parents:
        raise PublisherProteinRemapError(
            f"{path.name}: CDS parents are absent from the accepted {graph_mode} models: "
            + ",".join(orphan_parents[:5])
        )
    ambiguous = sorted(
        transcript_id
        for transcript_id, accessions in proteins_by_transcript.items()
        if len(accessions) != 1
    )
    if ambiguous:
        raise PublisherProteinRemapError(
            f"{path.name}: transcripts map to multiple protein accessions: "
            + ",".join(ambiguous[:5])
        )
    protein_by_transcript = {
        transcript_id: next(iter(accessions))
        for transcript_id, accessions in proteins_by_transcript.items()
    }
    transcripts_by_protein: dict[str, list[str]] = defaultdict(list)
    for transcript_id, protein_id in protein_by_transcript.items():
        transcripts_by_protein[protein_id].append(transcript_id)
    repeated_proteins = sorted(
        protein_id
        for protein_id, transcript_ids in transcripts_by_protein.items()
        if len(transcript_ids) != 1
    )
    if repeated_proteins:
        raise PublisherProteinRemapError(
            f"{path.name}: protein accessions map to multiple transcripts: "
            + ",".join(repeated_proteins[:5])
        )
    return _GffProteinMap(
        transcript_accessions=transcript_accessions,
        protein_by_transcript=protein_by_transcript,
        transcript_by_protein={
            protein_id: transcript_ids[0]
            for protein_id, transcript_ids in transcripts_by_protein.items()
        },
        transcript_count=len(model_ids),
        gene_count=gene_row_count,
        source_transcript_row_count=len(transcript_rows),
        noncoding_model_count=len(model_ids - set(protein_by_transcript)),
        graph_mode=graph_mode,
    )


def _write_tsv(
    path: Path, rows: Iterable[dict[str, object]], columns: tuple[str, ...]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_remapped_fasta(
    path: Path,
    selected_ids: Iterable[str],
    gff_map: _GffProteinMap,
    publisher_records: dict[str, _FastaRecord],
) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for transcript_id in selected_ids:
            protein_id = gff_map.protein_by_transcript[transcript_id]
            record = publisher_records[protein_id]
            mrna_accession = gff_map.transcript_accessions[transcript_id]
            handle.write(
                f">{transcript_id}\tpublisher_protein_id={protein_id}"
                f"\tpublisher_mRNA_accession={mrna_accession}\n"
            )
            handle.write(record.sequence + "\n")


def _validate_schema_names(names: Iterable[str], label: str) -> None:
    values = tuple(names)
    if not values or any(not SAFE_SCHEMA_NAME.fullmatch(value) for value in values):
        raise PublisherProteinRemapError(
            f"{label} values must be non-empty GFF3/header-safe names"
        )
    if len(set(values)) != len(values):
        raise PublisherProteinRemapError(f"{label} values must be unique")


def _validate_inputs(
    selected: Path,
    gff: Path,
    publisher: Path,
    output: Path,
    sample_id: str,
) -> None:
    if not SAFE_SAMPLE_ID.fullmatch(sample_id):
        raise PublisherProteinRemapError(
            "sample_id must start with an alphanumeric character and contain only letters, "
            "numbers, periods, underscores, or hyphens"
        )
    for role, path in (
        ("selected primary proteins", selected),
        ("GFF3", gff),
        ("publisher proteins", publisher),
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise PublisherProteinRemapError(
                f"The {role} input is missing or empty: {path.name}"
            )
    resolved = {selected.resolve(), gff.resolve(), publisher.resolve()}
    if len(resolved) != 3:
        raise PublisherProteinRemapError("The three input files must be distinct")
    if output.exists() or output.is_symlink():
        raise PublisherProteinRemapError(
            f"Output directory already exists; refusing overwrite: {output.name}"
        )


def remap_publisher_primary_proteins(
    selected_primary_proteins: str | Path,
    gff_path: str | Path,
    publisher_proteins: str | Path,
    output_dir: str | Path,
    sample_id: str,
    *,
    transcript_features: Iterable[str] = ("mRNA", "transcript"),
    gene_features: Iterable[str] = ("gene", "pseudogene"),
    transcript_id_attribute: str = "ID",
    transcript_accession_attribute: str = "Accession",
    transcript_accession_source: str = "attribute",
    cds_parent_attribute: str = "Parent",
    protein_accession_attribute: str = "Protein_Accession",
    protein_accession_source: str = "attribute",
    gene_as_transcript: bool = False,
    publisher_transcript_key: str = "OriID",
    publisher_mrna_accession_key: str = "mRNA",
    publisher_header_mode: str = "metadata",
) -> PublisherProteinRemapResult:
    """Create an exact, sequence-preserving publisher-primary protein subset."""
    selected_path = Path(selected_primary_proteins).expanduser()
    gff = Path(gff_path).expanduser()
    publisher_path = Path(publisher_proteins).expanduser()
    output = Path(output_dir).expanduser()
    transcript_types = frozenset(
        value.strip() for value in transcript_features if value.strip()
    )
    gene_types = frozenset(value.strip() for value in gene_features if value.strip())
    _validate_schema_names(transcript_types, "transcript feature")
    _validate_schema_names(gene_types, "gene feature")
    if transcript_types & gene_types:
        raise PublisherProteinRemapError(
            "transcript and gene feature types must be disjoint"
        )
    if transcript_accession_source not in TRANSCRIPT_ACCESSION_SOURCES:
        raise PublisherProteinRemapError(
            "transcript_accession_source must be one of: "
            + ", ".join(sorted(TRANSCRIPT_ACCESSION_SOURCES))
        )
    if protein_accession_source not in PROTEIN_ACCESSION_SOURCES:
        raise PublisherProteinRemapError(
            "protein_accession_source must be one of: "
            + ", ".join(sorted(PROTEIN_ACCESSION_SOURCES))
        )
    active_gff_schema_names = [transcript_id_attribute, cds_parent_attribute]
    if transcript_accession_source == "attribute":
        active_gff_schema_names.append(transcript_accession_attribute)
    if protein_accession_source == "attribute":
        active_gff_schema_names.append(protein_accession_attribute)
    if publisher_header_mode not in PUBLISHER_HEADER_MODES:
        raise PublisherProteinRemapError(
            "publisher_header_mode must be one of: "
            + ", ".join(sorted(PUBLISHER_HEADER_MODES))
        )
    publisher_schema_names = (publisher_transcript_key, publisher_mrna_accession_key)
    _validate_schema_names(active_gff_schema_names, "active GFF3 attribute")
    if publisher_header_mode == "metadata":
        _validate_schema_names(publisher_schema_names, "publisher header key")
    _validate_inputs(selected_path, gff, publisher_path, output, sample_id)

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        inputs = (
            ("selected_primary_proteins", selected_path),
            ("matched_GFF3", gff),
            ("publisher_proteins", publisher_path),
        )
        signatures = {path: _input_signature(path) for _, path in inputs}
        input_rows: list[dict[str, object]] = []
        for role, path in inputs:
            size, digest = _checksum(path)
            input_rows.append(
                {"role": role, "file_name": path.name, "bytes": size, "sha256": digest}
            )

        selected_records = _read_protein_fasta(selected_path, "selected primary")
        selected_ids = tuple(selected_records)
        gff_map = _parse_gff_protein_map(
            gff,
            transcript_features=transcript_types,
            gene_features=gene_types,
            transcript_id_attribute=transcript_id_attribute,
            transcript_accession_attribute=transcript_accession_attribute,
            transcript_accession_source=transcript_accession_source,
            cds_parent_attribute=cds_parent_attribute,
            protein_accession_attribute=protein_accession_attribute,
            protein_accession_source=protein_accession_source,
            gene_as_transcript=gene_as_transcript,
        )
        publisher_records = _read_protein_fasta(
            publisher_path,
            "publisher",
            required_header_keys=(
                publisher_schema_names
                if publisher_header_mode == "metadata"
                else ()
            ),
        )

        gff_protein_ids = set(gff_map.transcript_by_protein)
        publisher_ids = set(publisher_records)
        if gff_protein_ids != publisher_ids:
            gff_only = sorted(gff_protein_ids - publisher_ids)
            publisher_only = sorted(publisher_ids - gff_protein_ids)
            raise PublisherProteinRemapError(
                "GFF3/publisher protein accession sets are not exactly equal: "
                f"gff_only={len(gff_only)} ({','.join(gff_only[:5])}); "
                f"publisher_only={len(publisher_only)} ({','.join(publisher_only[:5])})"
            )

        if publisher_header_mode == "metadata":
            observed_header_transcripts: dict[str, str] = {}
            observed_header_mrnas: dict[str, str] = {}
            for protein_id, record in publisher_records.items():
                transcript_id = record.metadata[publisher_transcript_key]
                mrna_accession = record.metadata[publisher_mrna_accession_key]
                if transcript_id in observed_header_transcripts:
                    raise PublisherProteinRemapError(
                        f"Publisher header transcript ID {transcript_id!r} occurs for multiple proteins"
                    )
                if mrna_accession in observed_header_mrnas:
                    raise PublisherProteinRemapError(
                        f"Publisher header mRNA accession {mrna_accession!r} occurs "
                        "for multiple proteins"
                    )
                observed_header_transcripts[transcript_id] = protein_id
                observed_header_mrnas[mrna_accession] = protein_id
                expected_transcript = gff_map.transcript_by_protein[protein_id]
                expected_mrna = gff_map.transcript_accessions[expected_transcript]
                if transcript_id != expected_transcript or mrna_accession != expected_mrna:
                    raise PublisherProteinRemapError(
                        f"Publisher header mapping disagrees with GFF3 for protein {protein_id!r}: "
                        f"expected {expected_transcript!r}/{expected_mrna!r}, "
                        f"observed {transcript_id!r}/{mrna_accession!r}"
                    )

        missing_selected = sorted(set(selected_ids) - set(gff_map.protein_by_transcript))
        if missing_selected:
            raise PublisherProteinRemapError(
                "Selected primary transcript IDs lack exact GFF3 protein mappings: "
                + ",".join(missing_selected[:5])
            )

        prefix = sample_id
        protein_name = f"{prefix}.publisher_primary.remapped.protein.faa"
        mapping_name = f"{prefix}.publisher_primary.mapping.tsv"
        inventory_name = f"{prefix}.publisher_protein.source_inventory.tsv"
        summary_name = f"{prefix}.publisher_primary.summary.tsv"
        protein_output = staging / protein_name
        _write_remapped_fasta(
            protein_output, selected_ids, gff_map, publisher_records
        )

        output_records = _read_protein_fasta(protein_output, "remapped publisher-primary")
        if tuple(output_records) != selected_ids or set(output_records) != set(selected_ids):
            raise PublisherProteinRemapError(
                "Remapped output does not have exact selected-primary ID closure"
            )
        mapping_rows: list[dict[str, object]] = []
        for transcript_id in selected_ids:
            protein_id = gff_map.protein_by_transcript[transcript_id]
            publisher_record = publisher_records[protein_id]
            output_record = output_records[transcript_id]
            if publisher_record.sequence != output_record.sequence:
                raise PublisherProteinRemapError(
                    f"Publisher sequence changed during remapping for {transcript_id!r}"
                )
            mapping_rows.append(
                {
                    "selected_transcript_id": transcript_id,
                    "gff_transcript_accession": gff_map.transcript_accessions[transcript_id],
                    "publisher_protein_id": protein_id,
                    "publisher_header_transcript_id": publisher_record.metadata.get(
                        publisher_transcript_key, ""
                    ),
                    "publisher_header_mRNA_accession": publisher_record.metadata.get(
                        publisher_mrna_accession_key, ""
                    ),
                    "sequence_length": len(publisher_record.sequence),
                    "source_sequence_sha256": publisher_record.sequence_sha256,
                    "output_sequence_sha256": output_record.sequence_sha256,
                    "status": "PASS_EXACT_MAPPING_AND_SEQUENCE_PRESERVATION",
                }
            )
        _write_tsv(
            staging / mapping_name,
            mapping_rows,
            (
                "selected_transcript_id",
                "gff_transcript_accession",
                "publisher_protein_id",
                "publisher_header_transcript_id",
                "publisher_header_mRNA_accession",
                "sequence_length",
                "source_sequence_sha256",
                "output_sequence_sha256",
                "status",
            ),
        )

        selected_set = set(selected_ids)
        inventory_rows: list[dict[str, object]] = []
        for protein_id, publisher_record in publisher_records.items():
            transcript_id = gff_map.transcript_by_protein[protein_id]
            is_selected = transcript_id in selected_set
            inventory_rows.append(
                {
                    "publisher_protein_id": protein_id,
                    "gff_transcript_id": transcript_id,
                    "gff_transcript_accession": gff_map.transcript_accessions[transcript_id],
                    "publisher_header_transcript_id": publisher_record.metadata.get(
                        publisher_transcript_key, ""
                    ),
                    "publisher_header_mRNA_accession": publisher_record.metadata.get(
                        publisher_mrna_accession_key, ""
                    ),
                    "selected_primary": str(is_selected).lower(),
                    "disposition": (
                        "SELECTED_AND_REMAPPED" if is_selected else "EXCLUDED_NONPRIMARY"
                    ),
                    "sequence_length": len(publisher_record.sequence),
                    "sequence_sha256": publisher_record.sequence_sha256,
                    "status": "PASS",
                }
            )
        _write_tsv(
            staging / inventory_name,
            inventory_rows,
            (
                "publisher_protein_id",
                "gff_transcript_id",
                "gff_transcript_accession",
                "publisher_header_transcript_id",
                "publisher_header_mRNA_accession",
                "selected_primary",
                "disposition",
                "sequence_length",
                "sequence_sha256",
                "status",
            ),
        )

        excluded_count = len(publisher_records) - len(selected_ids)
        summary_rows = [
            {
                "sample_id": sample_id,
                "status": "PASS",
                "publication_gate": "PASS",
                "workflow_version": WORKFLOW_VERSION,
                "annotation_graph_mode": gff_map.graph_mode,
                "gene_as_transcript_requested": str(gene_as_transcript).lower(),
                "publisher_header_mode": publisher_header_mode,
                "publisher_header_mapping_check": (
                    "true" if publisher_header_mode == "metadata" else "not_applicable"
                ),
                "gff_transcript_record_count": gff_map.transcript_count,
                "gff_gene_record_count": gff_map.gene_count,
                "gff_source_transcript_row_count": gff_map.source_transcript_row_count,
                "gff_coding_transcript_count": len(gff_map.protein_by_transcript),
                "gff_noncoding_model_count": gff_map.noncoding_model_count,
                "source_publisher_record_count": len(publisher_records),
                "selected_primary_record_count": len(selected_ids),
                "output_record_count": len(output_records),
                "excluded_nonprimary_record_count": excluded_count,
                "exact_source_accession_closure": "true",
                "exact_selected_primary_closure": "true",
                "one_to_one_mapping": "true",
                "sequence_preservation": "true",
            }
        ]
        _write_tsv(
            staging / summary_name,
            summary_rows,
            (
                "sample_id",
                "status",
                "publication_gate",
                "workflow_version",
                "annotation_graph_mode",
                "gene_as_transcript_requested",
                "publisher_header_mode",
                "publisher_header_mapping_check",
                "gff_transcript_record_count",
                "gff_gene_record_count",
                "gff_source_transcript_row_count",
                "gff_coding_transcript_count",
                "gff_noncoding_model_count",
                "source_publisher_record_count",
                "selected_primary_record_count",
                "output_record_count",
                "excluded_nonprimary_record_count",
                "exact_source_accession_closure",
                "exact_selected_primary_closure",
                "one_to_one_mapping",
                "sequence_preservation",
            ),
        )
        _write_tsv(
            staging / "input_checksums.tsv",
            input_rows,
            ("role", "file_name", "bytes", "sha256"),
        )

        manifest = {
            "schema_version": 1,
            "workflow": "publisher_primary_protein_remap",
            "workflow_version": WORKFLOW_VERSION,
            "status": "PASS",
            "publication_gate": "PASS",
            "sample_id": sample_id,
            "execution": {"processes": 1, "worker_threads": 0},
            "inputs": input_rows,
            "schema": {
                "transcript_feature_types": sorted(transcript_types),
                "gene_feature_types": sorted(gene_types),
                "transcript_ID_attribute": transcript_id_attribute,
                "transcript_accession_attribute": transcript_accession_attribute,
                "transcript_accession_source": transcript_accession_source,
                "CDS_parent_attribute": cds_parent_attribute,
                "protein_accession_attribute": protein_accession_attribute,
                "protein_accession_source": protein_accession_source,
                "publisher_transcript_header_key": publisher_transcript_key,
                "publisher_mRNA_accession_header_key": publisher_mrna_accession_key,
                "publisher_header_mode": publisher_header_mode,
            },
            "policy": {
                "gene_as_transcript_requested": gene_as_transcript,
                "annotation_graph_mode": gff_map.graph_mode,
                "gene_as_transcript_policy": (
                    "explicit_only;zero_accepted_transcript_rows;one_declared_gene_parent_per_CDS;"
                    "preserve_gene_ID_as_self_transcript;reject_pseudogene_CDS"
                    if gene_as_transcript
                    else "disabled"
                ),
                "source_accession_policy": (
                    "exact_GFF3_coding_protein_accession_set_equals_publisher_FASTA_first_token_set"
                ),
                "mapping_policy": (
                    "one_to_one_GFF3_transcript_ID_to_transcript_accession_to_protein_accession"
                ),
                "header_policy": (
                    "publisher_header_transcript_and_mRNA_fields_must_match_GFF3"
                    if publisher_header_mode == "metadata"
                    else "publisher_first_token_is_protein_accession;transcript_mapping_comes_from_GFF3"
                ),
                "selection_policy": "exact_selected_primary_FASTA_first_token_ID_set",
                "sequence_policy": "preserve_every_sequence_character_without_normalization",
                "source_nonprimary_policy": "retain_in_audit_and_exclude_from_remapped_subset",
                "missing_ambiguous_or_extra_policy": "reject_complete_run",
            },
            "counts": {
                "GFF3_transcripts": gff_map.transcript_count,
                "GFF3_declared_gene_rows": gff_map.gene_count,
                "GFF3_source_transcript_rows": gff_map.source_transcript_row_count,
                "GFF3_coding_transcripts": len(gff_map.protein_by_transcript),
                "GFF3_noncoding_models": gff_map.noncoding_model_count,
                "source_publisher_proteins": len(publisher_records),
                "selected_primary_transcripts": len(selected_ids),
                "output_publisher_primary_proteins": len(output_records),
                "excluded_nonprimary_source_proteins": excluded_count,
            },
            "output_protein_file": protein_name,
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

        for path, signature in signatures.items():
            if _input_signature(path) != signature:
                raise PublisherProteinRemapError(
                    f"Input {path.name} changed during remapping; no output was published"
                )
        checksum_rows: list[dict[str, object]] = []
        for path in sorted(staging.iterdir(), key=lambda item: item.name):
            if path.is_file() and path.name != "checksums.tsv":
                size, digest = _checksum(path)
                checksum_rows.append({"file": path.name, "bytes": size, "sha256": digest})
        _write_tsv(staging / "checksums.tsv", checksum_rows, ("file", "bytes", "sha256"))
        if output.exists() or output.is_symlink():
            raise PublisherProteinRemapError(
                f"Output directory appeared during the run; refusing overwrite: {output.name}"
            )
        os.replace(staging, output)
        return PublisherProteinRemapResult(
            output_dir=output,
            source_publisher_record_count=len(publisher_records),
            selected_primary_record_count=len(selected_ids),
            excluded_nonprimary_record_count=excluded_count,
            output_protein_path=output / protein_name,
        )
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
