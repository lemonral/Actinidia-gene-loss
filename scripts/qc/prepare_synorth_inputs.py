#!/usr/bin/env python3
"""Build deterministic primary-protein and coordinate inputs for SynOrths.

The script accepts one version-matched GFF3 and protein FASTA (plain text or
gzip compressed), resolves the annotation relationships described below, and
selects the longest non-empty mapped protein for every gene.  A lexical sort
of source protein identifiers is the deterministic tie-break.

Supported relationship rules
----------------------------
1. ``CDS Parent=<gene-container>;protein_id=<protein>`` (direct CDS-to-gene).
2. A gene-container ID exactly matching an input FASTA identifier (the
   MAKER/no-transcript convention).
3. ``gene-container <- mRNA/transcript <- CDS protein_id``
   (transcript-mediated).

For CDS features, official GWH ``Protein_Accession`` is accepted as an exact
alias of NCBI/DDBJ ``protein_id``.  If both attributes occur on one row, their
single values must agree.  GWH transcript ``geneIDParent`` is likewise an alias
for ``Parent``.  Gene/transcript ``Accession`` values must be unique when
present, and ``Parent_Accession`` must match the resolved parent.

Both ``gene`` and ``pseudogene`` are explicit gene-container feature types.
They receive identical coordinate, duplicate-ID, and parent validation.  A
container enters the SynOrths outputs only when it has a mapped non-empty
protein record.

Identifiers are exact and case-sensitive.  FASTA identifiers are the first
whitespace-delimited token after ``>``; the program does not strip accession
versions or guess aliases.  Every non-excluded input FASTA record must resolve
to exactly one GFF3 gene.  A protein may not resolve through several
transcripts or genes.  GFF3-referenced proteins absent from the FASTA are
reported, while GFF3 genes without protein relationships are retained as QC
rows but cannot be written to the SynOrths inputs.

Output FASTA headers are replaced by the selected gene ID.  The header IDs and
the first column of the headerless coordinate table are checked for exact set
equality before any output is committed.

One or more explicitly named sequence IDs may be removed with repeated
``--exclude-seqid`` options.  Exclusion happens before feature-coordinate and
duplicate-ID validation, which permits known organellar multipart annotations
to be removed without weakening validation of the retained nuclear features.

Gene-container/transcript containment is strict by default.  The explicit
``--repair-gene-coordinates-from-children`` mode first requires every child
transcript to match its parent's sequence ID and strand, then expands each
gene/pseudogene interval to the deterministic union of the declared interval
and all child-transcript intervals.  Transcript/CDS containment remains strict.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, TextIO
from urllib.parse import unquote


TRANSCRIPT_FEATURE_TYPES = {"mrna", "transcript"}
GENE_CONTAINER_TYPES = {"gene", "pseudogene"}
MAPPING_COLUMNS = (
    "gene_id",
    "gene_feature_type",
    "source_transcript_id",
    "source_protein_id",
    "protein_length_aa",
    "selected",
    "status",
    "reason",
    "mapping_rule",
    "gene_coordinate_status",
    "gene_original_coordinates",
    "gene_output_coordinates",
    "excluded_seqids",
)


class PreparationError(RuntimeError):
    """Raised when an input or relationship is unsafe for SynOrths."""


@dataclass(frozen=True)
class Location:
    """A validated GFF3 genomic interval."""

    seqid: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class Gene:
    """One unique GFF3 gene feature."""

    gene_id: str
    feature_type: str
    accession: str
    original_location: Location
    location: Location
    line_number: int


@dataclass(frozen=True)
class Transcript:
    """One mRNA/transcript feature and its declared gene parent."""

    transcript_id: str
    parent_gene_id: str
    accession: str
    parent_accession: str
    location: Location
    line_number: int


@dataclass(frozen=True)
class CDS:
    """One protein-bearing CDS row."""

    parent_id: str
    parent_accession: str
    protein_id: str
    location: Location
    line_number: int


@dataclass(frozen=True)
class Protein:
    """One exact FASTA identifier and its whitespace-free sequence."""

    protein_id: str
    sequence: str
    header: str


@dataclass(frozen=True)
class Relationship:
    """A uniquely resolved protein-to-gene relationship."""

    gene_id: str
    transcript_id: str
    protein_id: str
    mapping_rule: str


@dataclass(frozen=True)
class MappingRow:
    """One auditable selection or exclusion record."""

    gene_id: str
    gene_feature_type: str
    source_transcript_id: str
    source_protein_id: str
    protein_length_aa: str
    selected: str
    status: str
    reason: str
    mapping_rule: str
    gene_coordinate_status: str
    gene_original_coordinates: str
    gene_output_coordinates: str
    excluded_seqids: str


@dataclass
class ParseCounts:
    """Counts captured while parsing the GFF3."""

    gff_feature_rows: int = 0
    gff_gene_rows: int = 0
    gff_pseudogene_rows: int = 0
    gff_gene_container_rows: int = 0
    gff_gene_container_rows_with_Accession: int = 0
    gff_transcript_rows: int = 0
    gff_transcript_rows_with_Accession: int = 0
    gff_transcript_rows_with_Parent: int = 0
    gff_transcript_rows_with_geneIDParent: int = 0
    gff_transcript_rows_with_Parent_Accession: int = 0
    gff_cds_rows: int = 0
    gff_cds_rows_with_protein_id: int = 0
    gff_cds_rows_with_lowercase_protein_id: int = 0
    gff_cds_rows_with_Protein_Accession: int = 0
    gff_cds_rows_with_Parent_Accession: int = 0
    gff_cds_rows_without_protein_id: int = 0
    gff_excluded_feature_rows: int = 0
    gff_excluded_gene_rows: int = 0
    gff_excluded_pseudogene_rows: int = 0
    gff_excluded_gene_container_rows: int = 0
    gff_excluded_transcript_rows: int = 0
    gff_excluded_cds_rows: int = 0
    gff_excluded_cds_rows_with_protein_id: int = 0
    gff_excluded_cds_rows_with_lowercase_protein_id: int = 0
    gff_excluded_cds_rows_with_Protein_Accession: int = 0
    gff_excluded_rows_with_unreadable_attributes: int = 0
    gene_containers_repaired_from_child_transcripts: int = 0
    transcript_rows_outside_declared_gene: int = 0


@dataclass
class ExcludedGffEvidence:
    """Identifiers observed on sequence IDs removed before strict validation."""

    observed_seqids: set[str]
    protein_seqids: dict[str, set[str]]
    direct_gene_seqids: dict[str, set[str]]
    direct_gene_types: dict[str, set[str]]


def utc_now() -> str:
    """Return a timezone-explicit timestamp."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    """Hash a file in bounded memory."""
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        raise PreparationError(f"Cannot hash {path}: {error}") from error
    return digest.hexdigest()


def is_gzip(path: Path) -> bool:
    """Detect gzip input by magic bytes, independent of filename suffix."""
    try:
        with path.open("rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError as error:
        raise PreparationError(f"Cannot inspect {path}: {error}") from error


def open_text_auto(path: Path) -> TextIO:
    """Open a plain-text or gzip-compressed input as UTF-8."""
    try:
        if is_gzip(path):
            return gzip.open(path, "rt", encoding="utf-8")
        return path.open("rt", encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise PreparationError(f"Cannot open text input {path}: {error}") from error


def require_input(path: Path, label: str) -> Path:
    """Resolve and validate a non-empty regular input file."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise PreparationError(f"{label} is not a regular file: {resolved}")
    try:
        if resolved.stat().st_size == 0:
            raise PreparationError(f"{label} is empty: {resolved}")
    except OSError as error:
        raise PreparationError(f"Cannot stat {label} {resolved}: {error}") from error
    return resolved


def validate_token(value: str, label: str, *, source: str) -> str:
    """Reject identifiers that cannot be represented safely in SynOrths text."""
    if not value:
        raise PreparationError(f"{source}: empty {label}")
    if any(character.isspace() for character in value):
        raise PreparationError(
            f"{source}: {label} contains whitespace and is unsafe for SynOrths: {value!r}"
        )
    return value


def parse_attributes(text: str, *, source: str) -> dict[str, list[str]]:
    """Parse GFF3 attributes into decoded, potentially multi-valued fields."""
    attributes: dict[str, list[str]] = defaultdict(list)
    if text == ".":
        return {}
    for item in text.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise PreparationError(f"{source}: malformed GFF3 attribute without '=': {item!r}")
        raw_key, raw_value = item.split("=", 1)
        key = unquote(raw_key)
        if not key:
            raise PreparationError(f"{source}: empty GFF3 attribute key")
        values = [unquote(value) for value in raw_value.split(",")]
        if any(value == "" for value in values):
            raise PreparationError(f"{source}: empty value for GFF3 attribute {key!r}")
        attributes[key].extend(values)
    return dict(attributes)


def parse_excluded_attributes(text: str) -> tuple[dict[str, list[str]], bool]:
    """Best-effort attribute parsing for rows already removed by sequence ID.

    Excluded rows must not fail the retained-annotation validation.  The
    boolean records whether any malformed attribute fragment was skipped.
    """
    attributes: dict[str, list[str]] = defaultdict(list)
    malformed = False
    if text == ".":
        return {}, False
    for item in text.split(";"):
        if not item:
            continue
        if "=" not in item:
            malformed = True
            continue
        raw_key, raw_value = item.split("=", 1)
        key = unquote(raw_key)
        values = [unquote(value) for value in raw_value.split(",")]
        if not key or any(value == "" for value in values):
            malformed = True
            continue
        attributes[key].extend(values)
    return dict(attributes), malformed


def one_attribute(
    attributes: dict[str, list[str]],
    key: str,
    *,
    source: str,
    required: bool,
) -> str | None:
    """Return one unique attribute value or reject an ambiguous field."""
    values = attributes.get(key, [])
    unique = sorted(set(values))
    if not unique:
        if required:
            raise PreparationError(f"{source}: required GFF3 attribute {key!r} is missing")
        return None
    if len(unique) != 1 or len(values) != 1:
        raise PreparationError(
            f"{source}: GFF3 attribute {key!r} must contain exactly one value; found {values!r}"
        )
    return unique[0]


def one_protein_identifier(
    attributes: dict[str, list[str]],
    *,
    source: str,
) -> tuple[str | None, set[str]]:
    """Resolve NCBI/DDBJ or GWH CDS protein-identifier attributes.

    The accepted attributes are ``protein_id`` and ``Protein_Accession``.
    When both occur they must contain the same single exact identifier.
    """
    observed_attributes = {
        key for key in ("protein_id", "Protein_Accession") if key in attributes
    }
    if not observed_attributes:
        return None, set()

    identifiers: dict[str, str] = {}
    for key in sorted(observed_attributes):
        value = one_attribute(attributes, key, source=source, required=True)
        assert value is not None
        identifiers[key] = value
    unique_identifiers = set(identifiers.values())
    if len(unique_identifiers) != 1:
        raise PreparationError(
            f"{source}: protein_id and Protein_Accession disagree: {identifiers!r}"
        )
    return next(iter(unique_identifiers)), observed_attributes


def one_aliased_identifier(
    attributes: dict[str, list[str]],
    keys: tuple[str, ...],
    *,
    label: str,
    source: str,
    required: bool,
) -> tuple[str | None, set[str]]:
    """Resolve one identifier from equivalent GFF3 attribute names."""
    observed_attributes = {key for key in keys if key in attributes}
    if not observed_attributes:
        if required:
            raise PreparationError(
                f"{source}: required {label} attribute is missing; accepted names are "
                f"{', '.join(keys)}"
            )
        return None, set()

    identifiers: dict[str, str] = {}
    for key in sorted(observed_attributes):
        value = one_attribute(attributes, key, source=source, required=True)
        assert value is not None
        identifiers[key] = value
    unique_identifiers = set(identifiers.values())
    if len(unique_identifiers) != 1:
        raise PreparationError(
            f"{source}: equivalent {label} attributes disagree: {identifiers!r}"
        )
    return next(iter(unique_identifiers)), observed_attributes


def parse_location(fields: list[str], *, source: str) -> Location:
    """Validate a gene/transcript/CDS interval and strand."""
    seqid = validate_token(fields[0], "sequence ID", source=source)
    try:
        start = int(fields[3])
        end = int(fields[4])
    except ValueError as error:
        raise PreparationError(
            f"{source}: start/end must be integers, found {fields[3]!r}/{fields[4]!r}"
        ) from error
    if start < 1 or end < start:
        raise PreparationError(
            f"{source}: invalid coordinates; require 1 <= start <= end, found {start}-{end}"
        )
    strand = fields[6]
    if strand not in {"+", "-"}:
        raise PreparationError(
            f"{source}: strand must be '+' or '-' for SynOrths, found {strand!r}"
        )
    return Location(seqid=seqid, start=start, end=end, strand=strand)


def require_child_location(child: Location, parent: Location, *, source: str) -> None:
    """Require a child feature to be contained on its parent's locus."""
    if child.seqid != parent.seqid or child.strand != parent.strand:
        raise PreparationError(
            f"{source}: child and parent sequence/strand disagree: "
            f"{child.seqid}:{child.strand} versus {parent.seqid}:{parent.strand}"
        )
    if child.start < parent.start or child.end > parent.end:
        raise PreparationError(
            f"{source}: child interval {child.start}-{child.end} is outside parent "
            f"interval {parent.start}-{parent.end}"
        )


def require_same_seqid_and_strand(
    child: Location,
    parent: Location,
    *,
    source: str,
) -> None:
    """Require a child and parent to share sequence ID and strand."""
    if child.seqid != parent.seqid or child.strand != parent.strand:
        raise PreparationError(
            f"{source}: child and parent sequence/strand disagree: "
            f"{child.seqid}:{child.strand} versus {parent.seqid}:{parent.strand}"
        )


def format_location(location: Location) -> str:
    """Format one exact interval for compact mapping/JSON audit fields."""
    return (
        f"{location.seqid}:{location.start}-{location.end}:{location.strand}"
    )


def gene_coordinate_mapping_fields(gene: Gene) -> dict[str, str]:
    """Return original/output coordinate fields for one mapping-QC row."""
    repaired = gene.location != gene.original_location
    return {
        "gene_coordinate_status": (
            "expanded_from_child_transcripts" if repaired else "unchanged"
        ),
        "gene_original_coordinates": format_location(gene.original_location),
        "gene_output_coordinates": format_location(gene.location),
    }


def require_parent_accession(
    provided: str,
    expected: str,
    *,
    source: str,
) -> None:
    """Validate an optional GWH Parent_Accession against its resolved parent."""
    if not provided:
        return
    if not expected:
        raise PreparationError(
            f"{source}: Parent_Accession {provided!r} is present but the resolved "
            "parent has no Accession"
        )
    if provided != expected:
        raise PreparationError(
            f"{source}: Parent_Accession {provided!r} does not match resolved parent "
            f"Accession {expected!r}"
        )


def parse_gff(
    path: Path,
    excluded_seqids: set[str],
    *,
    repair_gene_coordinates_from_children: bool,
) -> tuple[
    dict[str, Gene],
    dict[str, Transcript],
    list[CDS],
    ParseCounts,
    ExcludedGffEvidence,
]:
    """Parse supported GFF3 schemas and optionally repair gene intervals."""
    genes: dict[str, Gene] = {}
    transcripts: dict[str, Transcript] = {}
    cds_rows: list[CDS] = []
    counts = ParseCounts()
    excluded = ExcludedGffEvidence(
        observed_seqids=set(),
        protein_seqids=defaultdict(set),
        direct_gene_seqids=defaultdict(set),
        direct_gene_types=defaultdict(set),
    )

    with open_text_auto(path) as handle:
        try:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                source = f"{path}: line {line_number}"
                raw_seqid = fields[0] if fields else ""
                if raw_seqid in excluded_seqids:
                    excluded.observed_seqids.add(raw_seqid)
                    counts.gff_excluded_feature_rows += 1
                    feature_type = fields[2].lower() if len(fields) >= 3 else ""
                    if feature_type in GENE_CONTAINER_TYPES:
                        counts.gff_excluded_gene_container_rows += 1
                        if feature_type == "gene":
                            counts.gff_excluded_gene_rows += 1
                        else:
                            counts.gff_excluded_pseudogene_rows += 1
                    elif feature_type in TRANSCRIPT_FEATURE_TYPES:
                        counts.gff_excluded_transcript_rows += 1
                    elif feature_type == "cds":
                        counts.gff_excluded_cds_rows += 1

                    if len(fields) != 9:
                        counts.gff_excluded_rows_with_unreadable_attributes += 1
                        continue
                    attributes, malformed = parse_excluded_attributes(fields[8])
                    if malformed:
                        counts.gff_excluded_rows_with_unreadable_attributes += 1
                    if feature_type in GENE_CONTAINER_TYPES:
                        for gene_id in attributes.get("ID", []):
                            excluded.direct_gene_seqids[gene_id].add(raw_seqid)
                            excluded.direct_gene_types[gene_id].add(feature_type)
                    elif feature_type == "cds":
                        lowercase_ids = attributes.get("protein_id", [])
                        gwh_ids = attributes.get("Protein_Accession", [])
                        protein_ids = sorted(set(lowercase_ids).union(gwh_ids))
                        if protein_ids:
                            counts.gff_excluded_cds_rows_with_protein_id += 1
                        if lowercase_ids:
                            counts.gff_excluded_cds_rows_with_lowercase_protein_id += 1
                        if gwh_ids:
                            counts.gff_excluded_cds_rows_with_Protein_Accession += 1
                        for protein_id in protein_ids:
                            excluded.protein_seqids[protein_id].add(raw_seqid)
                    continue
                if len(fields) != 9:
                    raise PreparationError(
                        f"{source}: expected 9 tab-separated GFF3 columns, found {len(fields)}"
                    )
                counts.gff_feature_rows += 1
                feature_type = fields[2].lower()
                if feature_type not in {
                    *GENE_CONTAINER_TYPES,
                    "cds",
                    *TRANSCRIPT_FEATURE_TYPES,
                }:
                    continue
                location = parse_location(fields, source=source)
                attributes = parse_attributes(fields[8], source=source)

                if feature_type in GENE_CONTAINER_TYPES:
                    counts.gff_gene_container_rows += 1
                    if feature_type == "gene":
                        counts.gff_gene_rows += 1
                    else:
                        counts.gff_pseudogene_rows += 1
                    gene_id = one_attribute(attributes, "ID", source=source, required=True)
                    accession = one_attribute(
                        attributes,
                        "Accession",
                        source=source,
                        required=False,
                    )
                    assert gene_id is not None
                    validate_token(gene_id, "gene ID", source=source)
                    if accession is not None:
                        counts.gff_gene_container_rows_with_Accession += 1
                        validate_token(accession, "gene Accession", source=source)
                    if gene_id in genes:
                        previous = genes[gene_id]
                        raise PreparationError(
                            f"{source}: duplicate gene ID {gene_id!r}; first seen at line "
                            f"{previous.line_number}"
                        )
                    genes[gene_id] = Gene(
                        gene_id,
                        feature_type,
                        accession or "",
                        location,
                        location,
                        line_number,
                    )
                    continue

                if feature_type in TRANSCRIPT_FEATURE_TYPES:
                    counts.gff_transcript_rows += 1
                    transcript_id = one_attribute(attributes, "ID", source=source, required=True)
                    parent_gene, parent_attributes = one_aliased_identifier(
                        attributes,
                        ("Parent", "geneIDParent"),
                        label="transcript parent",
                        source=source,
                        required=True,
                    )
                    accession = one_attribute(
                        attributes,
                        "Accession",
                        source=source,
                        required=False,
                    )
                    parent_accession = one_attribute(
                        attributes,
                        "Parent_Accession",
                        source=source,
                        required=False,
                    )
                    assert transcript_id is not None and parent_gene is not None
                    validate_token(transcript_id, "transcript ID", source=source)
                    validate_token(parent_gene, "transcript Parent", source=source)
                    if "Parent" in parent_attributes:
                        counts.gff_transcript_rows_with_Parent += 1
                    if "geneIDParent" in parent_attributes:
                        counts.gff_transcript_rows_with_geneIDParent += 1
                    if accession is not None:
                        counts.gff_transcript_rows_with_Accession += 1
                        validate_token(accession, "transcript Accession", source=source)
                    if parent_accession is not None:
                        counts.gff_transcript_rows_with_Parent_Accession += 1
                        validate_token(
                            parent_accession,
                            "transcript Parent_Accession",
                            source=source,
                        )
                    if transcript_id in transcripts:
                        previous = transcripts[transcript_id]
                        raise PreparationError(
                            f"{source}: duplicate transcript ID {transcript_id!r}; first seen at "
                            f"line {previous.line_number}"
                        )
                    transcripts[transcript_id] = Transcript(
                        transcript_id,
                        parent_gene,
                        accession or "",
                        parent_accession or "",
                        location,
                        line_number,
                    )
                    continue

                counts.gff_cds_rows += 1
                protein_id, protein_attributes = one_protein_identifier(
                    attributes,
                    source=source,
                )
                if protein_id is None:
                    counts.gff_cds_rows_without_protein_id += 1
                    continue
                counts.gff_cds_rows_with_protein_id += 1
                if "protein_id" in protein_attributes:
                    counts.gff_cds_rows_with_lowercase_protein_id += 1
                if "Protein_Accession" in protein_attributes:
                    counts.gff_cds_rows_with_Protein_Accession += 1
                parent = one_attribute(attributes, "Parent", source=source, required=True)
                parent_accession = one_attribute(
                    attributes,
                    "Parent_Accession",
                    source=source,
                    required=False,
                )
                assert parent is not None
                validate_token(parent, "CDS Parent", source=source)
                if parent_accession is not None:
                    counts.gff_cds_rows_with_Parent_Accession += 1
                    validate_token(
                        parent_accession,
                        "CDS Parent_Accession",
                        source=source,
                    )
                validate_token(protein_id, "protein_id", source=source)
                cds_rows.append(
                    CDS(
                        parent,
                        parent_accession or "",
                        protein_id,
                        location,
                        line_number,
                    )
                )
        except (OSError, UnicodeError) as error:
            raise PreparationError(f"Cannot read GFF3 {path}: {error}") from error

    if not genes:
        raise PreparationError(f"GFF3 contains no gene or pseudogene features: {path}")

    gene_accessions: dict[str, str] = {}
    for gene in genes.values():
        if not gene.accession:
            continue
        previous_gene_id = gene_accessions.get(gene.accession)
        if previous_gene_id is not None:
            raise PreparationError(
                f"{path}: duplicate gene-container Accession {gene.accession!r} on "
                f"{previous_gene_id!r} and {gene.gene_id!r}"
            )
        gene_accessions[gene.accession] = gene.gene_id

    transcript_accessions: dict[str, str] = {}
    children_by_gene: dict[str, list[Transcript]] = defaultdict(list)
    for transcript in sorted(
        transcripts.values(),
        key=lambda item: item.transcript_id,
    ):
        source = f"{path}: line {transcript.line_number} transcript {transcript.transcript_id!r}"
        gene = genes.get(transcript.parent_gene_id)
        if gene is None:
            raise PreparationError(
                f"{source}: Parent gene {transcript.parent_gene_id!r} does not exist"
            )
        if repair_gene_coordinates_from_children:
            require_same_seqid_and_strand(
                transcript.location,
                gene.original_location,
                source=source,
            )
            children_by_gene[gene.gene_id].append(transcript)
            if (
                transcript.location.start < gene.original_location.start
                or transcript.location.end > gene.original_location.end
            ):
                counts.transcript_rows_outside_declared_gene += 1
        else:
            require_child_location(
                transcript.location,
                gene.original_location,
                source=source,
            )
        require_parent_accession(
            transcript.parent_accession,
            gene.accession,
            source=source,
        )
        if transcript.accession:
            previous_transcript_id = transcript_accessions.get(transcript.accession)
            if previous_transcript_id is not None:
                raise PreparationError(
                    f"{path}: duplicate transcript Accession {transcript.accession!r} on "
                    f"{previous_transcript_id!r} and {transcript.transcript_id!r}"
                )
            transcript_accessions[transcript.accession] = transcript.transcript_id

    if repair_gene_coordinates_from_children:
        for gene_id in sorted(children_by_gene):
            gene = genes[gene_id]
            children = children_by_gene[gene_id]
            repaired_start = min(
                [gene.original_location.start]
                + [transcript.location.start for transcript in children]
            )
            repaired_end = max(
                [gene.original_location.end]
                + [transcript.location.end for transcript in children]
            )
            repaired_location = Location(
                seqid=gene.original_location.seqid,
                start=repaired_start,
                end=repaired_end,
                strand=gene.original_location.strand,
            )
            if repaired_location != gene.original_location:
                genes[gene_id] = replace(gene, location=repaired_location)
                counts.gene_containers_repaired_from_child_transcripts += 1

        for transcript in sorted(
            transcripts.values(),
            key=lambda item: item.transcript_id,
        ):
            source = (
                f"{path}: line {transcript.line_number} transcript "
                f"{transcript.transcript_id!r}"
            )
            require_child_location(
                transcript.location,
                genes[transcript.parent_gene_id].location,
                source=source,
            )

    return genes, transcripts, cds_rows, counts, excluded


def read_proteins(path: Path) -> dict[str, Protein]:
    """Read protein FASTA records and reject duplicate/unsafe identifiers."""
    proteins: dict[str, Protein] = {}
    current_id: str | None = None
    current_header = ""
    sequence_parts: list[str] = []

    def commit() -> None:
        nonlocal current_id, current_header, sequence_parts
        if current_id is None:
            return
        if current_id in proteins:
            raise PreparationError(f"{path}: duplicate FASTA identifier {current_id!r}")
        proteins[current_id] = Protein(
            protein_id=current_id,
            sequence="".join(sequence_parts),
            header=current_header,
        )

    with open_text_auto(path) as handle:
        try:
            for line_number, raw_line in enumerate(handle, start=1):
                stripped = raw_line.strip()
                if not stripped:
                    continue
                if stripped.startswith(">"):
                    commit()
                    header = stripped[1:].strip()
                    if not header:
                        raise PreparationError(f"{path}: empty FASTA header at line {line_number}")
                    protein_id = header.split(maxsplit=1)[0]
                    validate_token(
                        protein_id,
                        "protein FASTA identifier",
                        source=f"{path}: line {line_number}",
                    )
                    current_id = protein_id
                    current_header = header
                    sequence_parts = []
                    continue
                if current_id is None:
                    raise PreparationError(
                        f"{path}: sequence data precedes the first FASTA header at line "
                        f"{line_number}"
                    )
                sequence_parts.append("".join(stripped.split()))
        except (OSError, UnicodeError) as error:
            raise PreparationError(f"Cannot read protein FASTA {path}: {error}") from error
    commit()
    if not proteins:
        raise PreparationError(f"Protein FASTA contains no records: {path}")
    return proteins


def resolve_relationships(
    gff_path: Path,
    genes: dict[str, Gene],
    transcripts: dict[str, Transcript],
    cds_rows: list[CDS],
    proteins: dict[str, Protein],
    excluded_evidence: ExcludedGffEvidence,
) -> tuple[
    dict[str, Relationship],
    dict[str, int],
    dict[str, set[str]],
    set[str],
]:
    """Resolve CDS/direct mappings and reject every ambiguous input protein."""
    explicit: dict[str, set[tuple[str, str, str]]] = defaultdict(set)
    relationship_counts: dict[str, int] = defaultdict(int)

    for cds in cds_rows:
        source = f"{gff_path}: line {cds.line_number} CDS protein_id={cds.protein_id!r}"
        if cds.parent_id in genes:
            gene = genes[cds.parent_id]
            require_child_location(cds.location, gene.location, source=source)
            require_parent_accession(
                cds.parent_accession,
                gene.accession,
                source=source,
            )
            explicit[cds.protein_id].add((gene.gene_id, "", "cds_parent_gene"))
            relationship_counts["cds_parent_gene_rows"] += 1
            continue
        transcript = transcripts.get(cds.parent_id)
        if transcript is None:
            raise PreparationError(
                f"{source}: Parent {cds.parent_id!r} is neither a gene nor an mRNA/transcript"
            )
        gene = genes[transcript.parent_gene_id]
        require_child_location(cds.location, transcript.location, source=source)
        require_child_location(cds.location, gene.location, source=source)
        require_parent_accession(
            cds.parent_accession,
            transcript.accession,
            source=source,
        )
        explicit[cds.protein_id].add(
            (gene.gene_id, transcript.transcript_id, "cds_parent_transcript")
        )
        relationship_counts["cds_parent_transcript_rows"] += 1

    direct_ids = set(proteins).intersection(genes)
    relationships: dict[str, Relationship] = {}
    for protein_id in sorted(set(explicit).union(direct_ids)):
        candidates = explicit.get(protein_id, set())
        if candidates:
            gene_ids = {candidate[0] for candidate in candidates}
            transcript_ids = {candidate[1] for candidate in candidates if candidate[1]}
            if len(gene_ids) != 1:
                raise PreparationError(
                    f"Protein {protein_id!r} maps through CDS rows to multiple genes: "
                    f"{sorted(gene_ids)!r}"
                )
            if len(transcript_ids) > 1:
                raise PreparationError(
                    f"Protein {protein_id!r} maps through several transcripts: "
                    f"{sorted(transcript_ids)!r}; each source protein must identify one isoform"
                )
            gene_id = next(iter(gene_ids))
            if protein_id in direct_ids and protein_id != gene_id:
                raise PreparationError(
                    f"Protein {protein_id!r} exactly matches gene ID {protein_id!r} but its "
                    f"CDS relationship points to different gene {gene_id!r}"
                )
            if protein_id in direct_ids and protein_id == gene_id:
                relationship_counts["direct_gene_id_confirmed_by_cds"] += 1
            transcript_id = next(iter(transcript_ids), "")
            rule = "cds_parent_transcript" if transcript_id else "cds_parent_gene"
            relationships[protein_id] = Relationship(
                gene_id=gene_id,
                transcript_id=transcript_id,
                protein_id=protein_id,
                mapping_rule=rule,
            )
        else:
            relationships[protein_id] = Relationship(
                gene_id=protein_id,
                transcript_id="",
                protein_id=protein_id,
                mapping_rule="exact_gene_id",
            )
            relationship_counts["exact_gene_id_proteins"] += 1

    excluded_protein_seqids = {
        protein_id: set(seqids)
        for protein_id, seqids in excluded_evidence.protein_seqids.items()
    }
    for gene_id, seqids in excluded_evidence.direct_gene_seqids.items():
        if gene_id in proteins:
            excluded_protein_seqids.setdefault(gene_id, set()).update(seqids)

    excluded_only_input_ids = (
        set(proteins).intersection(excluded_protein_seqids).difference(relationships)
    )
    unmapped = sorted(
        set(proteins).difference(relationships).difference(excluded_only_input_ids)
    )
    if unmapped:
        preview = ", ".join(repr(identifier) for identifier in unmapped[:10])
        suffix = "" if len(unmapped) <= 10 else f" ... and {len(unmapped) - 10} more"
        raise PreparationError(
            f"{len(unmapped)} FASTA protein identifier(s) do not map by any supported GFF3 "
            f"rule: {preview}{suffix}"
        )

    return (
        relationships,
        dict(relationship_counts),
        excluded_protein_seqids,
        excluded_only_input_ids,
    )


def select_primary_proteins(
    genes: dict[str, Gene],
    proteins: dict[str, Protein],
    relationships: dict[str, Relationship],
    excluded_protein_seqids: dict[str, set[str]],
    excluded_only_input_ids: set[str],
    excluded_direct_gene_types: dict[str, set[str]],
) -> tuple[dict[str, tuple[Protein, Relationship]], list[MappingRow], dict[str, int]]:
    """Select one longest non-empty protein per gene and construct the QC audit."""
    by_gene: dict[str, list[tuple[Protein, Relationship]]] = defaultdict(list)
    for protein_id, protein in proteins.items():
        if protein_id in excluded_only_input_ids:
            continue
        relationship = relationships[protein_id]
        by_gene[relationship.gene_id].append((protein, relationship))

    selected: dict[str, tuple[Protein, Relationship]] = {}
    rows: list[MappingRow] = []
    genes_with_only_empty = 0
    discarded_isoforms = 0

    for gene_id in sorted(by_gene):
        candidates = by_gene[gene_id]
        nonempty = [candidate for candidate in candidates if candidate[0].sequence]
        winner: tuple[Protein, Relationship] | None = None
        if nonempty:
            winner = sorted(
                nonempty,
                key=lambda item: (-len(item[0].sequence), item[0].protein_id),
            )[0]
            selected[gene_id] = winner
        else:
            genes_with_only_empty += 1

        winner_length = len(winner[0].sequence) if winner else 0
        winner_id = winner[0].protein_id if winner else ""
        for protein, relationship in sorted(candidates, key=lambda item: item[0].protein_id):
            length = len(protein.sequence)
            if winner is not None and protein.protein_id == winner_id:
                selected_flag = "yes"
                status = "selected"
                reason = "longest_nonempty"
            elif length == 0:
                selected_flag = "no"
                status = "not_selected"
                reason = "empty_sequence"
                discarded_isoforms += 1
            elif length < winner_length:
                selected_flag = "no"
                status = "not_selected"
                reason = "shorter_than_selected"
                discarded_isoforms += 1
            else:
                selected_flag = "no"
                status = "not_selected"
                reason = "lexicographic_tie_break"
                discarded_isoforms += 1
            rows.append(
                MappingRow(
                    gene_id=gene_id,
                    gene_feature_type=genes[gene_id].feature_type,
                    source_transcript_id=relationship.transcript_id,
                    source_protein_id=protein.protein_id,
                    protein_length_aa=str(length),
                    selected=selected_flag,
                    status=status,
                    reason=reason,
                    mapping_rule=relationship.mapping_rule,
                    **gene_coordinate_mapping_fields(genes[gene_id]),
                    excluded_seqids=";".join(
                        sorted(excluded_protein_seqids.get(protein.protein_id, set()))
                    ),
                )
            )

    referenced_protein_ids = set(relationships)
    missing_referenced = sorted(referenced_protein_ids.difference(proteins))
    for protein_id in missing_referenced:
        relationship = relationships[protein_id]
        rows.append(
            MappingRow(
                gene_id=relationship.gene_id,
                gene_feature_type=genes[relationship.gene_id].feature_type,
                source_transcript_id=relationship.transcript_id,
                source_protein_id=protein_id,
                protein_length_aa="",
                selected="no",
                status="not_available",
                reason="referenced_protein_missing_from_fasta",
                mapping_rule=relationship.mapping_rule,
                **gene_coordinate_mapping_fields(genes[relationship.gene_id]),
                excluded_seqids=";".join(
                    sorted(excluded_protein_seqids.get(protein_id, set()))
                ),
            )
        )

    genes_with_relationship = {relationship.gene_id for relationship in relationships.values()}
    genes_without_relationship = sorted(set(genes).difference(genes_with_relationship))
    for gene_id in genes_without_relationship:
        rows.append(
            MappingRow(
                gene_id=gene_id,
                gene_feature_type=genes[gene_id].feature_type,
                source_transcript_id="",
                source_protein_id="",
                protein_length_aa="",
                selected="no",
                status="not_available",
                reason="gene_has_no_protein_relationship",
                mapping_rule="",
                **gene_coordinate_mapping_fields(genes[gene_id]),
                excluded_seqids="",
            )
        )

    excluded_only_annotation_ids = sorted(
        set(excluded_protein_seqids).difference(relationships)
    )
    excluded_missing_from_fasta = 0
    for protein_id in excluded_only_annotation_ids:
        protein = proteins.get(protein_id)
        if protein is None:
            excluded_missing_from_fasta += 1
        direct_gene_match = protein_id in excluded_direct_gene_types
        rows.append(
            MappingRow(
                gene_id=protein_id if direct_gene_match else "",
                gene_feature_type=(
                    ";".join(sorted(excluded_direct_gene_types[protein_id]))
                    if direct_gene_match
                    else ""
                ),
                source_transcript_id="",
                source_protein_id=protein_id,
                protein_length_aa=str(len(protein.sequence)) if protein else "",
                selected="no",
                status="excluded",
                reason=(
                    "excluded_seqid"
                    if protein is not None
                    else "excluded_seqid_referenced_protein_missing_from_fasta"
                ),
                mapping_rule=(
                    "excluded_exact_gene_id" if direct_gene_match else "excluded_seqid"
                ),
                gene_coordinate_status="",
                gene_original_coordinates="",
                gene_output_coordinates="",
                excluded_seqids=";".join(
                    sorted(excluded_protein_seqids[protein_id])
                ),
            )
        )

    if len(selected) != len(set(selected)):
        raise PreparationError("Duplicate output gene identifiers were produced")

    rows.sort(
        key=lambda row: (
            row.gene_id,
            row.source_protein_id,
            row.source_transcript_id,
            row.reason,
        )
    )
    selection_counts = {
        "selected_primary_proteins": len(selected),
        "genes_with_mapped_proteins": len(by_gene),
        "genes_with_only_empty_proteins": genes_with_only_empty,
        "genes_without_protein_relationship": len(genes_without_relationship),
        "gff_referenced_proteins_missing_from_fasta": len(missing_referenced),
        "discarded_input_isoforms": discarded_isoforms,
        "excluded_annotation_protein_ids": len(excluded_only_annotation_ids),
        "excluded_input_protein_records": len(excluded_only_input_ids),
        "excluded_referenced_proteins_missing_from_fasta": excluded_missing_from_fasta,
        "included_input_proteins_also_seen_on_excluded_seqids": len(
            set(proteins).intersection(relationships).intersection(excluded_protein_seqids)
        ),
    }
    return selected, rows, selection_counts


def resolve_outputs(args: argparse.Namespace, inputs: Iterable[Path]) -> dict[str, Path]:
    """Resolve output paths, enforce uniqueness, and create parent directories."""
    outputs = {
        "protein_fasta": args.output_protein.expanduser().resolve(),
        "coordinate_table": args.output_coords.expanduser().resolve(),
        "mapping_qc": args.output_mapping.expanduser().resolve(),
        "summary_json": args.summary_json.expanduser().resolve(),
    }
    input_set = set(inputs)
    if len(set(outputs.values())) != len(outputs):
        raise PreparationError("Every output option must name a different file")
    collisions = sorted(path for path in outputs.values() if path in input_set)
    if collisions:
        raise PreparationError(
            "Output paths must not overwrite input files: "
            + ", ".join(str(path) for path in collisions)
        )
    for label, path in outputs.items():
        if path.exists() and not args.force:
            raise PreparationError(
                f"Refusing to replace existing {label} without --force: {path}"
            )
        if path.exists() and path.is_dir():
            raise PreparationError(f"Output path is a directory: {path}")
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            raise PreparationError(
                f"Cannot create output directory {path.parent}: {error}"
            ) from error
    return outputs


def temporary_path(path: Path) -> Path:
    """Allocate a same-directory temporary filename for atomic replacement."""
    return path.with_name(f".{path.name}.tmp.{os.getpid()}")


def write_primary_fasta(path: Path, selected: dict[str, tuple[Protein, Relationship]]) -> None:
    """Write selected sequences with exact gene-only headers."""
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for gene_id in sorted(selected):
            protein, _relationship = selected[gene_id]
            handle.write(f">{gene_id}\n")
            for offset in range(0, len(protein.sequence), 60):
                handle.write(protein.sequence[offset : offset + 60] + "\n")


def write_coords(
    path: Path,
    selected: dict[str, tuple[Protein, Relationship]],
    genes: dict[str, Gene],
) -> None:
    """Write the five-column, headerless SynOrths coordinate table."""
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for gene_id in sorted(selected):
            location = genes[gene_id].location
            handle.write(
                f"{gene_id}\t{location.seqid}\t{location.start}\t{location.end}\t"
                f"{location.strand}\n"
            )


def write_mapping(path: Path, rows: list[MappingRow]) -> None:
    """Write the complete mapping and selection audit."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MAPPING_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def exact_output_id_check(
    selected: dict[str, tuple[Protein, Relationship]], genes: dict[str, Gene]
) -> None:
    """Assert the intended FASTA and coordinate ID sets are identical."""
    fasta_ids = set(selected)
    coordinate_ids = {gene_id for gene_id in selected if gene_id in genes}
    if len(fasta_ids) != len(selected):
        raise PreparationError("Duplicate output FASTA identifiers were produced")
    if fasta_ids != coordinate_ids:
        missing_coords = sorted(fasta_ids.difference(coordinate_ids))
        missing_fasta = sorted(coordinate_ids.difference(fasta_ids))
        raise PreparationError(
            "Output FASTA and coordinate identifiers differ; "
            f"missing coordinates={missing_coords!r}, missing FASTA={missing_fasta!r}"
        )


def commit_outputs(
    outputs: dict[str, Path],
    selected: dict[str, tuple[Protein, Relationship]],
    genes: dict[str, Gene],
    rows: list[MappingRow],
    summary: dict[str, object],
) -> None:
    """Write temporary files, add checksums, then atomically publish all outputs."""
    temporary = {label: temporary_path(path) for label, path in outputs.items()}
    try:
        write_primary_fasta(temporary["protein_fasta"], selected)
        write_coords(temporary["coordinate_table"], selected, genes)
        write_mapping(temporary["mapping_qc"], rows)
        summary["outputs"] = {
            label: {
                "path": str(outputs[label]),
                "size_bytes": temporary[label].stat().st_size,
                "sha256": sha256_file(temporary[label]),
            }
            for label in ("protein_fasta", "coordinate_table", "mapping_qc")
        }
        with temporary["summary_json"].open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
        for label in ("protein_fasta", "coordinate_table", "mapping_qc", "summary_json"):
            os.replace(temporary[label], outputs[label])
    except (OSError, UnicodeError) as error:
        raise PreparationError(f"Cannot write outputs: {error}") from error
    finally:
        for path in temporary.values():
            try:
                path.unlink(missing_ok=True)
            except OSError:
                pass


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Select one deterministic primary protein per GFF3 gene and create "
            "ID-matched SynOrths protein/coordinate inputs."
        )
    )
    parser.add_argument(
        "--gff",
        required=True,
        type=Path,
        help="Version-matched GFF3 (plain/gzip)",
    )
    parser.add_argument(
        "--protein-fasta", required=True, type=Path, help="Version-matched proteins (plain/gzip)"
    )
    parser.add_argument(
        "--output-protein", required=True, type=Path, help="Selected primary FASTA"
    )
    parser.add_argument(
        "--output-coords", required=True, type=Path, help="Headerless SynOrths coords"
    )
    parser.add_argument("--output-mapping", required=True, type=Path, help="Mapping/QC TSV")
    parser.add_argument(
        "--summary-json", required=True, type=Path, help="Counts and provenance JSON"
    )
    parser.add_argument(
        "--exclude-seqid",
        action="append",
        default=[],
        metavar="SEQID",
        help=(
            "Exact GFF3 sequence ID to remove before feature validation; repeat for "
            "multiple organellar/other sequences"
        ),
    )
    parser.add_argument(
        "--repair-gene-coordinates-from-children",
        action="store_true",
        help=(
            "Opt in to expanding each gene/pseudogene interval to the union of its "
            "declared interval and all same-sequence/same-strand child transcripts"
        ),
    )
    parser.add_argument("--force", action="store_true", help="Replace all four existing outputs")
    return parser


def run(args: argparse.Namespace) -> dict[str, object]:
    """Execute preparation and return the summary dictionary."""
    gff_path = require_input(args.gff, "GFF3")
    protein_path = require_input(args.protein_fasta, "Protein FASTA")
    if gff_path == protein_path:
        raise PreparationError("GFF3 and protein FASTA must be different files")
    excluded_seqids = set()
    for seqid in args.exclude_seqid:
        validate_token(seqid, "excluded sequence ID", source="--exclude-seqid")
        excluded_seqids.add(seqid)
    outputs = resolve_outputs(args, (gff_path, protein_path))

    genes, transcripts, cds_rows, parse_counts, excluded_evidence = parse_gff(
        gff_path,
        excluded_seqids,
        repair_gene_coordinates_from_children=(
            args.repair_gene_coordinates_from_children
        ),
    )
    offending_children_by_gene: dict[str, list[str]] = defaultdict(list)
    if args.repair_gene_coordinates_from_children:
        for transcript in transcripts.values():
            gene = genes[transcript.parent_gene_id]
            if (
                transcript.location.start < gene.original_location.start
                or transcript.location.end > gene.original_location.end
            ):
                offending_children_by_gene[gene.gene_id].append(
                    transcript.transcript_id
                )
    coordinate_repairs = [
        {
            "gene_id": gene_id,
            "gene_feature_type": genes[gene_id].feature_type,
            "original_coordinates": format_location(
                genes[gene_id].original_location
            ),
            "output_coordinates": format_location(genes[gene_id].location),
            "child_transcripts_outside_declared_gene": sorted(
                offending_children_by_gene[gene_id]
            ),
        }
        for gene_id in sorted(genes)
        if genes[gene_id].location != genes[gene_id].original_location
    ]
    if len(coordinate_repairs) != parse_counts.gene_containers_repaired_from_child_transcripts:
        raise AssertionError("Coordinate-repair gene count is internally inconsistent")
    if (
        sum(len(items) for items in offending_children_by_gene.values())
        != parse_counts.transcript_rows_outside_declared_gene
    ):
        raise AssertionError("Coordinate-repair transcript count is internally inconsistent")
    proteins = read_proteins(protein_path)
    (
        relationships,
        relationship_counts,
        excluded_protein_seqids,
        excluded_only_input_ids,
    ) = resolve_relationships(
        gff_path,
        genes,
        transcripts,
        cds_rows,
        proteins,
        excluded_evidence,
    )
    selected, mapping_rows, selection_counts = select_primary_proteins(
        genes,
        proteins,
        relationships,
        excluded_protein_seqids,
        excluded_only_input_ids,
        excluded_evidence.direct_gene_types,
    )
    if not selected:
        raise PreparationError(
            "No non-empty mapped proteins remain after primary-isoform selection"
        )
    exact_output_id_check(selected, genes)

    summary_counts: dict[str, int] = {
        **asdict(parse_counts),
        **relationship_counts,
        "unique_gff_gene_containers": len(genes),
        "unique_gff_genes": sum(
            gene.feature_type == "gene" for gene in genes.values()
        ),
        "unique_gff_pseudogenes": sum(
            gene.feature_type == "pseudogene" for gene in genes.values()
        ),
        "unique_gff_transcripts": len(transcripts),
        "input_protein_records": len(proteins),
        "empty_input_protein_records": sum(not protein.sequence for protein in proteins.values()),
        "mapped_input_protein_records": len(proteins) - len(excluded_only_input_ids),
        "unmapped_input_protein_records": 0,
        "requested_excluded_seqid_count": len(excluded_seqids),
        "observed_excluded_seqid_count": len(excluded_evidence.observed_seqids),
        "gff_excluded_unique_cds_protein_ids": len(
            excluded_evidence.protein_seqids
        ),
        "gff_excluded_unique_gene_container_ids": len(
            excluded_evidence.direct_gene_types
        ),
        "gff_excluded_unique_gene_ids": sum(
            "gene" in feature_types
            for feature_types in excluded_evidence.direct_gene_types.values()
        ),
        "gff_excluded_unique_pseudogene_ids": sum(
            "pseudogene" in feature_types
            for feature_types in excluded_evidence.direct_gene_types.values()
        ),
        "protein_ids_associated_with_excluded_seqids": len(
            excluded_protein_seqids
        ),
        "input_protein_records_seen_on_excluded_seqids": len(
            set(proteins).intersection(excluded_protein_seqids)
        ),
        **selection_counts,
        "mapping_qc_rows": len(mapping_rows),
    }
    summary: dict[str, object] = {
        "status": "completed",
        "created_utc": utc_now(),
        "selection_rule": (
            "longest non-empty mapped protein per gene; lexicographically smallest exact "
            "source protein identifier breaks equal-length ties"
        ),
        "identifier_rule": (
            "exact case-sensitive identifiers; FASTA ID is the first token after '>'; "
            "output FASTA and coordinate IDs are GFF3 gene-container IDs"
        ),
        "gene_container_rule": {
            "supported_feature_types": sorted(GENE_CONTAINER_TYPES),
            "selection_requirement": "mapped non-empty protein record",
            "validation": (
                "gene and pseudogene containers share strict duplicate-ID, coordinate, "
                "strand, transcript-parent, and CDS-parent checks"
            ),
        },
        "coordinate_repair": {
            "enabled": args.repair_gene_coordinates_from_children,
            "default_policy": "strict containment; no repair",
            "opt_in_rule": (
                "on the same sequence ID and strand, expand each gene-container "
                "interval to the union of its declared interval and all child "
                "transcript intervals; do not alter transcript or CDS intervals"
            ),
            "repaired_gene_container_count": len(coordinate_repairs),
            "offending_child_transcript_count": sum(
                len(transcript_ids)
                for transcript_ids in offending_children_by_gene.values()
            ),
            "repairs": coordinate_repairs,
        },
        "exclusions": {
            "requested_seqids": sorted(excluded_seqids),
            "observed_seqids": sorted(excluded_evidence.observed_seqids),
            "requested_but_not_observed_seqids": sorted(
                excluded_seqids.difference(excluded_evidence.observed_seqids)
            ),
            "rule": (
                "exact sequence-ID match before duplicate, coordinate, strand, and "
                "parent-location validation"
            ),
        },
        "inputs": {
            "gff": {
                "path": str(gff_path),
                "size_bytes": gff_path.stat().st_size,
                "sha256": sha256_file(gff_path),
            },
            "protein_fasta": {
                "path": str(protein_path),
                "size_bytes": protein_path.stat().st_size,
                "sha256": sha256_file(protein_path),
            },
        },
        "counts": summary_counts,
    }
    commit_outputs(outputs, selected, genes, mapping_rows, summary)
    return summary


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        summary = run(args)
    except PreparationError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    counts = summary["counts"]
    assert isinstance(counts, dict)
    print(
        "Prepared "
        f"{counts['selected_primary_proteins']} primary proteins from "
        f"{counts['input_protein_records']} input records and "
        f"{counts['unique_gff_gene_containers']} GFF3 gene containers."
    )
    print(
        "QC: "
        f"{counts['discarded_input_isoforms']} non-primary/empty input records; "
        f"{counts['excluded_input_protein_records']} explicitly excluded input records; "
        f"{counts['gff_referenced_proteins_missing_from_fasta']} referenced proteins absent; "
        f"{counts['genes_without_protein_relationship']} genes without protein relationships."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
