"""Fail-closed primary-isoform extraction from chromosome-scope FASTA/GFF3.

This module is the production annotation standardizer.  It validates every
transcript before isoform selection, selects one valid coding transcript per
gene deterministically, and publishes a complete output directory atomically.
When gffread is available, it runs on the normalized selected-primary GFF3
(normally gene, mRNA, and CDS rows; top-level self-mRNA and CDS rows in the
explicit gene-as-transcript compatibility mode). Its independently extracted
CDS and protein sequences must agree with the Python implementation for every
selected ID.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, TextIO
from urllib.parse import quote, unquote

from .annotation import build_spliced_cds, translate_standard
from .gff import GffFeature, Transcript
from .io_utils import natural_key


WORKFLOW_VERSION = "1.2.0"
GZIP_MAGIC = b"\x1f\x8b"
SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
DNA_ALPHABET = set("ACGTRYSWKMBDHVN")


class PrimaryAnnotationError(RuntimeError):
    """Raised when a standardized annotation bundle cannot be published."""


@dataclass(frozen=True)
class CanonicalRule:
    attribute: str
    value: str | None

    @property
    def label(self) -> str:
        return self.attribute if self.value is None else f"{self.attribute}={self.value}"

    def matches(self, attributes: dict[str, str]) -> bool:
        if self.attribute not in attributes:
            return False
        observed = attributes[self.attribute].strip()
        if self.value is None:
            return bool(observed)
        return self.value == observed or self.value in {
            token.strip() for token in observed.split(",") if token.strip()
        }


@dataclass(frozen=True)
class TranscriptRecord:
    gene_id: str
    transcript_id: str
    seqid: str
    start: int
    end: int
    strand: str
    attributes: dict[str, str]
    line_number: int
    cds_features: tuple[GffFeature, ...]
    gene_as_transcript: bool = False

    @property
    def genomic_span(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class ValidTranscript:
    model: TranscriptRecord
    cds: str
    protein: str
    flags: tuple[str, ...]
    canonical_rule_index: int | None
    canonical_rule_label: str


@dataclass(frozen=True)
class StandardizationResult:
    output_dir: Path
    source_gene_count: int
    selected_gene_count: int
    invalid_coding_gene_count: int
    selected_transcript_count: int
    gffread_status: str


@dataclass(frozen=True)
class ParsedAnnotation:
    transcripts: tuple[TranscriptRecord, ...]
    gene_ids: frozenset[str]
    gene_locations: dict[str, tuple[str, int, int, str]]
    graph_mode: str
    source_transcript_row_count: int


def parse_canonical_rule(value: str) -> CanonicalRule:
    """Parse ``ATTRIBUTE`` or ``ATTRIBUTE=VALUE`` from the command line."""
    raw = value.strip()
    if not raw:
        raise PrimaryAnnotationError("A canonical-tag rule must not be empty")
    attribute, separator, expected = raw.partition("=")
    attribute = attribute.strip()
    if not attribute or any(character.isspace() for character in attribute):
        raise PrimaryAnnotationError(
            f"Invalid canonical-tag attribute in {value!r}; use ATTRIBUTE or ATTRIBUTE=VALUE"
        )
    if separator and not expected.strip():
        raise PrimaryAnnotationError(f"Canonical-tag rule {value!r} has an empty value")
    return CanonicalRule(attribute, expected.strip() if separator else None)


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == GZIP_MAGIC
    except OSError as error:
        raise PrimaryAnnotationError(f"Cannot inspect input {path.name}: {error}") from error


def _open_text(path: Path) -> TextIO:
    try:
        if _is_gzip(path):
            return gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="")
        return path.open("rt", encoding="utf-8", errors="strict", newline="")
    except OSError as error:
        raise PrimaryAnnotationError(f"Cannot open input {path.name}: {error}") from error


def _checksum(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                size += len(block)
                digest.update(block)
    except OSError as error:
        raise PrimaryAnnotationError(f"Cannot checksum input {path.name}: {error}") from error
    return size, digest.hexdigest()


def _input_signature(path: Path) -> tuple[int, int, int, int]:
    try:
        status = path.stat()
    except OSError as error:
        raise PrimaryAnnotationError(f"Cannot stat input {path.name}: {error}") from error
    return status.st_dev, status.st_ino, status.st_size, status.st_mtime_ns


def _validate_identifier(identifier: str, label: str, line_number: int, source: str) -> None:
    if not identifier or any(character.isspace() or ord(character) < 32 for character in identifier):
        raise PrimaryAnnotationError(
            f"{source}:{line_number}: {label} must be a non-empty identifier without whitespace"
        )


def _parse_attributes(raw: str, line_number: int, source: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    if raw == ".":
        return attributes
    for field in raw.split(";"):
        if not field:
            continue
        if "=" not in field:
            raise PrimaryAnnotationError(
                f"{source}:{line_number}: expected GFF3 key=value attributes; found {field!r}"
            )
        encoded_key, encoded_value = field.split("=", 1)
        key = unquote(encoded_key.strip())
        value = unquote(encoded_value.strip())
        if not key:
            raise PrimaryAnnotationError(f"{source}:{line_number}: empty GFF3 attribute key")
        if key in attributes:
            raise PrimaryAnnotationError(
                f"{source}:{line_number}: duplicate GFF3 attribute {key!r}"
            )
        attributes[key] = value
    return attributes


def _parents(attributes: dict[str, str]) -> list[str]:
    return [item.strip() for item in attributes.get("Parent", "").split(",") if item.strip()]


def load_genome(path: Path) -> dict[str, str]:
    """Load a chromosome-scope FASTA with strict ID and sequence validation."""
    records: dict[str, str] = {}
    identifier: str | None = None
    parts: list[str] = []

    def finish() -> None:
        nonlocal identifier, parts
        if identifier is None:
            return
        sequence = "".join(parts).upper()
        if not sequence:
            raise PrimaryAnnotationError(f"{path.name}: FASTA record {identifier!r} is empty")
        invalid = sorted(set(sequence) - DNA_ALPHABET)
        if invalid:
            raise PrimaryAnnotationError(
                f"{path.name}: FASTA record {identifier!r} contains invalid DNA symbols: "
                + ",".join(invalid[:10])
            )
        records[identifier] = sequence

    try:
        with _open_text(path) as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    finish()
                    identifier = line[1:].split()[0] if line[1:].split() else ""
                    _validate_identifier(identifier, "FASTA identifier", line_number, path.name)
                    if identifier in records:
                        raise PrimaryAnnotationError(
                            f"{path.name}:{line_number}: duplicate FASTA identifier {identifier!r}"
                        )
                    parts = []
                elif identifier is None:
                    raise PrimaryAnnotationError(
                        f"{path.name}:{line_number}: sequence occurs before the first FASTA header"
                    )
                else:
                    parts.append("".join(line.split()))
        finish()
    except (UnicodeError, EOFError, gzip.BadGzipFile) as error:
        raise PrimaryAnnotationError(f"Cannot read FASTA {path.name}: {error}") from error
    if not records:
        raise PrimaryAnnotationError(f"{path.name}: no FASTA records found")
    return records


def parse_annotation(
    path: Path,
    genome_lengths: dict[str, int],
    transcript_features: frozenset[str],
    gene_features: frozenset[str],
    gene_as_transcript: bool = False,
) -> ParsedAnnotation:
    """Parse a GFF3 coding graph and reject structural ambiguity.

    The normal graph is gene -> transcript -> CDS.  An explicitly enabled
    compatibility mode accepts gene -> CDS only when the *entire* source GFF3
    has no accepted transcript rows.  In that mode each declared gene becomes
    one self-transcript whose ID is the unchanged publisher gene ID.
    """
    source = path.name
    transcript_rows: dict[str, tuple[str, int, int, str, dict[str, str], int]] = {}
    cds_by_parent: dict[str, list[GffFeature]] = defaultdict(list)
    gene_ids: set[str] = set()
    gene_locations: dict[str, tuple[str, int, int, str]] = {}
    gene_rows: dict[
        str, tuple[str, str, int, int, str, dict[str, str], int]
    ] = {}
    cds_parent_groups: list[tuple[int, tuple[str, ...]]] = []
    saw_feature = False
    saw_embedded_fasta = False
    try:
        with _open_text(path) as handle:
            for line_number, raw in enumerate(handle, start=1):
                if raw.startswith("##FASTA"):
                    saw_embedded_fasta = True
                    break
                if not raw.strip() or raw.startswith("#"):
                    continue
                saw_feature = True
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 9:
                    raise PrimaryAnnotationError(
                        f"{source}:{line_number}: expected 9 GFF3 columns, found {len(fields)}"
                    )
                seqid, _, feature_type, start_text, end_text, _, strand, phase, raw_attrs = fields
                if seqid not in genome_lengths:
                    raise PrimaryAnnotationError(
                        f"{source}:{line_number}: GFF3 seqid {seqid!r} is absent from chromosome FASTA"
                    )
                try:
                    start, end = int(start_text), int(end_text)
                except ValueError as error:
                    raise PrimaryAnnotationError(
                        f"{source}:{line_number}: non-integer GFF3 interval {start_text!r}-{end_text!r}"
                    ) from error
                if start < 1 or end < start or end > genome_lengths[seqid]:
                    raise PrimaryAnnotationError(
                        f"{source}:{line_number}: interval {seqid}:{start}-{end} is outside the FASTA"
                    )
                if strand not in {"+", "-", ".", "?"}:
                    raise PrimaryAnnotationError(
                        f"{source}:{line_number}: invalid strand {strand!r}"
                    )
                attrs = _parse_attributes(raw_attrs, line_number, source)
                identifier = attrs.get("ID", "")
                if feature_type in gene_features:
                    _validate_identifier(identifier, f"{feature_type} ID", line_number, source)
                    if identifier in gene_locations:
                        raise PrimaryAnnotationError(
                            f"{source}:{line_number}: duplicate gene-level ID {identifier!r}"
                        )
                    gene_ids.add(identifier)
                    gene_locations[identifier] = (seqid, start, end, strand)
                    gene_rows[identifier] = (
                        feature_type,
                        seqid,
                        start,
                        end,
                        strand,
                        attrs,
                        line_number,
                    )
                elif feature_type in transcript_features:
                    _validate_identifier(identifier, "transcript ID", line_number, source)
                    parents = _parents(attrs)
                    if len(parents) != 1:
                        raise PrimaryAnnotationError(
                            f"{source}:{line_number}: transcript {identifier!r} must have exactly one Parent"
                        )
                    _validate_identifier(parents[0], "transcript Parent", line_number, source)
                    if identifier in transcript_rows:
                        raise PrimaryAnnotationError(
                            f"{source}:{line_number}: duplicate transcript ID {identifier!r}"
                        )
                    transcript_rows[identifier] = (seqid, start, end, strand, attrs, line_number)
                    gene_ids.add(parents[0])
                elif feature_type == "CDS":
                    parents = _parents(attrs)
                    if not parents:
                        raise PrimaryAnnotationError(
                            f"{source}:{line_number}: CDS must have at least one transcript Parent"
                        )
                    cds_parent_groups.append((line_number, tuple(parents)))
                    for parent in parents:
                        _validate_identifier(parent, "CDS Parent", line_number, source)
                        cds_by_parent[parent].append(
                            GffFeature(
                                sequence_id=seqid,
                                feature_type="CDS",
                                start=start,
                                end=end,
                                strand=strand,
                                phase=phase,
                                attributes=attrs,
                                line_number=line_number,
                            )
                        )
    except (UnicodeError, EOFError, gzip.BadGzipFile) as error:
        raise PrimaryAnnotationError(f"Cannot read GFF3 {source}: {error}") from error
    if saw_embedded_fasta:
        raise PrimaryAnnotationError(
            f"{source}: embedded FASTA is not accepted; provide a separate chromosome-scope FASTA"
        )
    if not saw_feature:
        raise PrimaryAnnotationError(f"{source}: no GFF3 feature rows found")
    if not transcript_rows and not gene_as_transcript:
        raise PrimaryAnnotationError(
            f"{source}: no transcript rows found for feature types {sorted(transcript_features)}"
        )
    if not transcript_rows:
        if not gene_rows:
            raise PrimaryAnnotationError(
                f"{source}: --gene-as-transcript requires declared gene-level rows"
            )
        for line_number, parents in cds_parent_groups:
            if len(parents) != 1:
                raise PrimaryAnnotationError(
                    f"{source}:{line_number}: gene-as-transcript CDS must have exactly one "
                    "declared gene Parent"
                )
            parent = parents[0]
            if parent not in gene_rows:
                raise PrimaryAnnotationError(
                    f"{source}:{line_number}: gene-as-transcript CDS Parent {parent!r} "
                    "does not identify a declared gene-level feature"
                )
            feature_type = gene_rows[parent][0]
            if feature_type == "pseudogene":
                raise PrimaryAnnotationError(
                    f"{source}:{line_number}: CDS Parent {parent!r} is declared as a "
                    "pseudogene; refusing to synthesize a coding transcript"
                )

        synthesized: list[TranscriptRecord] = []
        for gene_id, (
            _, seqid, start, end, strand, attrs, line_number
        ) in gene_rows.items():
            synthesized.append(
                TranscriptRecord(
                    gene_id=gene_id,
                    transcript_id=gene_id,
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    attributes=attrs,
                    line_number=line_number,
                    cds_features=tuple(cds_by_parent.get(gene_id, [])),
                    gene_as_transcript=True,
                )
            )
        return ParsedAnnotation(
            transcripts=tuple(synthesized),
            gene_ids=frozenset(gene_ids),
            gene_locations=gene_locations,
            graph_mode="gene_as_transcript",
            source_transcript_row_count=0,
        )

    identifier_collisions = sorted(set(gene_locations) & set(transcript_rows), key=natural_key)
    if identifier_collisions:
        raise PrimaryAnnotationError(
            f"{source}: gene and transcript IDs share the same GFF3 ID namespace: "
            + ", ".join(identifier_collisions[:5])
        )
    orphan_cds = sorted(set(cds_by_parent) - set(transcript_rows), key=natural_key)
    if orphan_cds:
        preview = ", ".join(orphan_cds[:5])
        raise PrimaryAnnotationError(
            f"{source}: CDS Parent IDs do not identify declared transcripts: {preview}"
        )

    transcripts: list[TranscriptRecord] = []
    for transcript_id, (seqid, start, end, strand, attrs, line_number) in transcript_rows.items():
        gene_id = _parents(attrs)[0]
        gene_location = gene_locations.get(gene_id)
        if gene_location is not None:
            gene_seqid, gene_start, gene_end, gene_strand = gene_location
            if gene_seqid != seqid or start < gene_start or end > gene_end:
                raise PrimaryAnnotationError(
                    f"{source}:{line_number}: transcript {transcript_id!r} lies outside parent gene {gene_id!r}"
                )
            if gene_strand in {"+", "-"} and gene_strand != strand:
                raise PrimaryAnnotationError(
                    f"{source}:{line_number}: transcript {transcript_id!r} strand differs from parent gene"
                )
        transcripts.append(
            TranscriptRecord(
                gene_id=gene_id,
                transcript_id=transcript_id,
                seqid=seqid,
                start=start,
                end=end,
                strand=strand,
                attributes=attrs,
                line_number=line_number,
                cds_features=tuple(cds_by_parent.get(transcript_id, [])),
                gene_as_transcript=False,
            )
        )
    return ParsedAnnotation(
        transcripts=tuple(transcripts),
        gene_ids=frozenset(gene_ids),
        gene_locations=gene_locations,
        graph_mode="declared_transcripts",
        source_transcript_row_count=len(transcript_rows),
    )


def _canonical_match(
    attributes: dict[str, str], rules: tuple[CanonicalRule, ...]
) -> tuple[int | None, str]:
    for index, rule in enumerate(rules):
        if rule.matches(attributes):
            return index, rule.label
    return None, ""


def validate_transcript(
    record: TranscriptRecord,
    genome: dict[str, str],
    canonical_rules: tuple[CanonicalRule, ...],
    missing_phase_policy: str,
) -> tuple[ValidTranscript | None, str, tuple[str, ...]]:
    """Return a valid sequence-bearing candidate or an explicit rejection reason."""
    rule_index, rule_label = _canonical_match(record.attributes, canonical_rules)
    model_flags = ("gene_as_transcript",) if record.gene_as_transcript else ()
    if not record.cds_features:
        return None, "no_CDS", model_flags
    if record.strand not in {"+", "-"}:
        return None, "transcript_strand_is_not_plus_or_minus", model_flags
    normalized_features: list[GffFeature] = []
    for feature in record.cds_features:
        if feature.sequence_id != record.seqid:
            return None, "CDS_seqid_differs_from_transcript", model_flags
        if feature.strand != record.strand:
            return None, "CDS_strand_differs_from_transcript", model_flags
        if feature.start < record.start or feature.end > record.end:
            return None, "CDS_interval_outside_transcript", model_flags
        phase = feature.phase
        if phase in {"", "."}:
            if missing_phase_policy == "fail":
                return None, "missing_CDS_phase", model_flags
            phase = "0"
        if phase not in {"0", "1", "2"}:
            return None, f"invalid_CDS_phase:{phase}", model_flags
        normalized_features.append(replace(feature, phase=phase))
    genomic_order = sorted(normalized_features, key=lambda item: (item.start, item.end))
    for previous, current in zip(genomic_order, genomic_order[1:]):
        if current.start <= previous.end:
            return None, "overlapping_or_duplicate_CDS_intervals", model_flags
    transcript = Transcript(
        gene_id=record.gene_id,
        transcript_id=record.transcript_id,
        chrom=record.seqid,
        start=record.start,
        end=record.end,
        strand=record.strand,
        cds_features=tuple(normalized_features),
    )
    try:
        cds = build_spliced_cds(transcript, genome)
    except (ValueError, RuntimeError) as error:
        return None, f"sequence_extraction_failed:{error}", model_flags
    if not cds:
        return None, "empty_spliced_CDS", model_flags
    invalid_symbols = sorted(set(cds) - DNA_ALPHABET)
    if invalid_symbols:
        return None, "invalid_CDS_symbols:" + ",".join(invalid_symbols), model_flags
    if len(cds) % 3:
        return None, "spliced_CDS_length_not_divisible_by_three", model_flags
    protein_with_stop = translate_standard(cds)
    if not protein_with_stop:
        return None, "empty_translation", model_flags
    if "*" in protein_with_stop[:-1]:
        return None, "internal_stop_codon", model_flags
    protein = protein_with_stop[:-1] if protein_with_stop.endswith("*") else protein_with_stop
    if not protein:
        return None, "translation_contains_only_stop", model_flags
    flags = list(model_flags)
    ordered = sorted(
        normalized_features, key=lambda item: item.start, reverse=record.strand == "-"
    )
    if ordered[0].phase != "0":
        flags.append("five_prime_partial_phase")
    if not protein.startswith("M"):
        flags.append("non_methionine_start")
    if not protein_with_stop.endswith("*"):
        flags.append("no_terminal_stop")
    if "X" in protein:
        flags.append("ambiguous_codon")
    if missing_phase_policy == "zero" and any(
        feature.phase in {"", "."} for feature in record.cds_features
    ):
        flags.append("missing_phase_treated_as_zero")
    return (
        ValidTranscript(record, cds, protein, tuple(flags), rule_index, rule_label),
        "valid",
        tuple(flags),
    )


def select_primary_isoforms(
    annotation: ParsedAnnotation,
    genome: dict[str, str],
    canonical_rules: tuple[CanonicalRule, ...],
    missing_phase_policy: str,
) -> tuple[list[ValidTranscript], list[dict[str, object]], list[dict[str, object]], int]:
    """Validate all isoforms first, then choose one valid candidate per gene."""
    by_gene: dict[str, list[TranscriptRecord]] = defaultdict(list)
    evaluations: dict[str, tuple[ValidTranscript | None, str, tuple[str, ...]]] = {}
    for transcript in annotation.transcripts:
        by_gene[transcript.gene_id].append(transcript)
        evaluations[transcript.transcript_id] = validate_transcript(
            transcript, genome, canonical_rules, missing_phase_policy
        )

    selected: list[ValidTranscript] = []
    selected_by_gene: dict[str, ValidTranscript] = {}
    selection_reason_by_gene: dict[str, str] = {}
    invalid_coding_gene_count = 0
    gene_rows: list[dict[str, object]] = []
    for gene_id in sorted(annotation.gene_ids, key=natural_key):
        models = by_gene.get(gene_id, [])
        valid = [
            evaluations[model.transcript_id][0]
            for model in models
            if evaluations[model.transcript_id][0] is not None
        ]
        valid = [item for item in valid if item is not None]
        models_with_cds = sum(bool(model.cds_features) for model in models)
        chosen: ValidTranscript | None = None
        selection_reason = ""
        if valid:
            canonical_indices = [
                item.canonical_rule_index
                for item in valid
                if item.canonical_rule_index is not None
            ]
            if canonical_indices:
                best_index = min(index for index in canonical_indices if index is not None)
                choices = [item for item in valid if item.canonical_rule_index == best_index]
                selection_reason = (
                    f"canonical_tag:{canonical_rules[best_index].label};"
                    "then_longest_valid_spliced_CDS;then_genomic_span;then_transcript_ID"
                )
            else:
                choices = valid
                selection_reason = (
                    "longest_valid_spliced_CDS;then_genomic_span;then_transcript_ID"
                )
            chosen = sorted(
                choices,
                key=lambda item: (
                    -len(item.cds),
                    -item.model.genomic_span,
                    item.model.transcript_id,
                ),
            )[0]
            selected.append(chosen)
            selected_by_gene[gene_id] = chosen
            selection_reason_by_gene[gene_id] = selection_reason
            gene_status = "selected"
        elif models_with_cds:
            invalid_coding_gene_count += 1
            gene_status = "no_valid_coding_transcript"
        elif models:
            gene_status = "no_CDS_transcript"
        else:
            gene_status = "no_transcript"
        gene_rows.append(
            {
                "gene_id": gene_id,
                "annotated_transcript_count": len(models),
                "transcripts_with_CDS_count": models_with_cds,
                "valid_coding_transcript_count": len(valid),
                "invalid_coding_transcript_count": sum(
                    bool(model.cds_features) and evaluations[model.transcript_id][0] is None
                    for model in models
                ),
                "selected_transcript_id": chosen.model.transcript_id if chosen else "",
                "selection_rule": selection_reason,
                "status": gene_status,
            }
        )

    transcript_rows: list[dict[str, object]] = []
    for model in sorted(
        annotation.transcripts,
        key=lambda item: (natural_key(item.seqid), item.start, item.transcript_id),
    ):
        candidate, validation_status, flags = evaluations[model.transcript_id]
        chosen = selected_by_gene.get(model.gene_id)
        if candidate is None:
            disposition = "invalid" if model.cds_features else "noncoding_or_no_CDS"
        elif chosen and chosen.model.transcript_id == model.transcript_id:
            disposition = "selected"
        else:
            disposition = "valid_not_selected"
        transcript_rows.append(
            {
                "gene_id": model.gene_id,
                "transcript_id": model.transcript_id,
                "seqid": model.seqid,
                "start": model.start,
                "end": model.end,
                "strand": model.strand,
                "raw_CDS_length": sum(
                    feature.end - feature.start + 1 for feature in model.cds_features
                ),
                "validated_spliced_CDS_length": len(candidate.cds) if candidate else "",
                "protein_length": len(candidate.protein) if candidate else "",
                "canonical_rule": candidate.canonical_rule_label if candidate else _canonical_match(
                    model.attributes, canonical_rules
                )[1],
                "validation_status": validation_status,
                "QC_flags": ";".join(flags),
                "disposition": disposition,
                "gene_selection_rule": selection_reason_by_gene.get(model.gene_id, ""),
            }
        )
    selected.sort(
        key=lambda item: (
            natural_key(item.model.seqid),
            item.model.start,
            item.model.transcript_id,
        )
    )
    return selected, transcript_rows, gene_rows, invalid_coding_gene_count


def _write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index:index + 60] + "\n")


def _write_tsv(path: Path, rows: Iterable[dict[str, object]], columns: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _encode_attribute(value: str) -> str:
    return quote(value, safe="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.:_+-|@")


def _write_primary_gff(
    path: Path,
    selected: list[ValidTranscript],
    gene_locations: dict[str, tuple[str, int, int, str]],
    graph_mode: str,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("##gff-version 3\n")
        for selected_model in selected:
            model = selected_model.model
            gene = _encode_attribute(model.gene_id)
            transcript = _encode_attribute(model.transcript_id)
            source = "primary_isoform_standardizer"
            gene_seqid, gene_start, gene_end, gene_strand = gene_locations.get(
                model.gene_id,
                (model.seqid, model.start, model.end, model.strand),
            )
            if graph_mode == "gene_as_transcript":
                # The publisher gene ID is also the self-transcript ID.  A gene
                # row with the same ID would violate the global GFF3 ID
                # namespace, so the normalized comparison object deliberately
                # uses a top-level mRNA plus CDS rows.  This retains the exact
                # publisher ID used in FASTA and audit outputs.
                handle.write(
                    f"{model.seqid}\t{source}\tmRNA\t{model.start}\t{model.end}\t.\t"
                    f"{model.strand}\t.\tID={transcript};gene_as_transcript=true\n"
                )
            else:
                handle.write(
                    f"{gene_seqid}\t{source}\tgene\t{gene_start}\t{gene_end}\t.\t"
                    f"{gene_strand}\t.\tID={gene}\n"
                )
                handle.write(
                    f"{model.seqid}\t{source}\tmRNA\t{model.start}\t{model.end}\t.\t"
                    f"{model.strand}\t.\tID={transcript};Parent={gene}\n"
                )
            for feature in sorted(model.cds_features, key=lambda item: (item.start, item.end)):
                phase = "0" if feature.phase in {"", "."} else feature.phase
                handle.write(
                    f"{model.seqid}\t{source}\tCDS\t{feature.start}\t{feature.end}\t.\t"
                    f"{model.strand}\t{phase}\tParent={transcript}\n"
                )


def _resolve_executable(specification: str) -> Path | None:
    if specification == "none":
        return None
    if specification == "auto":
        discovered = shutil.which("gffread")
        return Path(discovered) if discovered else None
    candidate = Path(specification).expanduser()
    if candidate.parent != Path(".") or candidate.is_absolute():
        if not candidate.is_file() or not os.access(candidate, os.X_OK):
            raise PrimaryAnnotationError("The configured gffread executable is missing or not executable")
        return candidate
    discovered = shutil.which(specification)
    if not discovered:
        raise PrimaryAnnotationError("The configured gffread command was not found on PATH")
    return Path(discovered)


def _materialize_plain_input(source: Path, destination: Path) -> Path:
    if not _is_gzip(source):
        return source
    try:
        with gzip.open(source, "rb") as input_handle, destination.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle, length=8 * 1024 * 1024)
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise PrimaryAnnotationError(f"Cannot stage compressed input {source.name}: {error}") from error
    return destination


def _read_output_fasta(path: Path, label: str) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    parts: list[str] = []

    def finish() -> None:
        if identifier is None:
            return
        if identifier in records:
            raise PrimaryAnnotationError(f"gffread produced duplicate {label} ID {identifier!r}")
        records[identifier] = "".join(parts).upper()

    try:
        with path.open("r", encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    finish()
                    tokens = line[1:].split()
                    identifier = tokens[0] if tokens else ""
                    if not identifier:
                        raise PrimaryAnnotationError(f"gffread produced an empty {label} FASTA ID")
                    parts = []
                elif identifier is None:
                    raise PrimaryAnnotationError(
                        f"gffread {label} FASTA contains sequence before a header"
                    )
                else:
                    parts.append("".join(line.split()))
        finish()
    except OSError as error:
        raise PrimaryAnnotationError(f"Cannot read gffread {label} output: {error}") from error
    return records


def compare_with_gffread(
    executable: Path,
    genome_path: Path,
    selected_primary_gff_path: Path,
    selected: list[ValidTranscript],
    staging_dir: Path,
    comparison_annotation_scope: str,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Run gffread on selected-primary GFF3 and require exact sequence agreement."""
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    version_process = subprocess.run(
        [str(executable), "--version"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=environment,
        check=False,
    )
    version_text = (version_process.stdout + "\n" + version_process.stderr).strip()
    if version_process.returncode != 0 or not version_text:
        raise PrimaryAnnotationError("gffread version query failed; no output was published")
    version = " ".join(version_text.splitlines()[0].split())[:200]
    # Tool versions belong in the public manifest, executable locations do not.
    version = version.replace(str(executable), executable.name)
    version = re.sub(r"(?<![A-Za-z0-9])/(?:[^\s]+)", "<path>", version)
    selected_gff_signature = _input_signature(selected_primary_gff_path)
    selected_gff_bytes, selected_gff_sha256 = _checksum(selected_primary_gff_path)
    plain_genome = _materialize_plain_input(genome_path, staging_dir / "gffread_input.fa")
    plain_gff = _materialize_plain_input(
        selected_primary_gff_path, staging_dir / "gffread_selected_primary.gff3"
    )
    cds_output = staging_dir / "gffread.cds.fa"
    protein_output = staging_dir / "gffread.protein.faa"
    stdout_path = staging_dir / "gffread.stdout.log"
    stderr_path = staging_dir / "gffread.stderr.log"
    with stdout_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
        process = subprocess.run(
            [
                str(executable),
                str(plain_gff),
                "-g",
                str(plain_genome),
                "-x",
                str(cds_output),
                "-y",
                str(protein_output),
            ],
            stdout=stdout_handle,
            stderr=stderr_handle,
            env=environment,
            check=False,
        )
    if process.returncode != 0:
        raise PrimaryAnnotationError(
            f"gffread sequence extraction returned status {process.returncode}; no output was published"
        )
    if _input_signature(selected_primary_gff_path) != selected_gff_signature:
        raise PrimaryAnnotationError(
            "Selected-primary GFF3 changed during gffread validation; no output was published"
        )
    selected_gff_bytes_after, selected_gff_sha256_after = _checksum(
        selected_primary_gff_path
    )
    if (
        selected_gff_bytes_after != selected_gff_bytes
        or selected_gff_sha256_after != selected_gff_sha256
    ):
        raise PrimaryAnnotationError(
            "Selected-primary GFF3 content changed during gffread validation; "
            "no output was published"
        )
    if not cds_output.is_file() or not protein_output.is_file():
        raise PrimaryAnnotationError("gffread did not create both CDS and protein FASTA outputs")
    gffread_cds = _read_output_fasta(cds_output, "CDS")
    gffread_proteins = _read_output_fasta(protein_output, "protein")
    selected_ids = {item.model.transcript_id for item in selected}
    if len(selected_ids) != len(selected):
        raise PrimaryAnnotationError(
            "Selected transcript IDs are not unique; no output was published"
        )
    for sequence_type, observed_records in (
        ("CDS", gffread_cds),
        ("protein", gffread_proteins),
    ):
        observed_ids = set(observed_records)
        if observed_ids != selected_ids:
            missing = sorted(selected_ids - observed_ids, key=natural_key)
            extra = sorted(observed_ids - selected_ids, key=natural_key)
            raise PrimaryAnnotationError(
                f"gffread {sequence_type} ID set does not exactly equal the selected-primary "
                f"set: missing={len(missing)} ({','.join(missing[:5])}); "
                f"extra={len(extra)} ({','.join(extra[:5])}); no output was published"
            )
    comparison_rows: list[dict[str, object]] = []
    discrepancies: list[str] = []
    for item in selected:
        transcript_id = item.model.transcript_id
        for sequence_type, ours, observed_records in (
            ("CDS", item.cds, gffread_cds),
            ("protein", item.protein.rstrip("*"), gffread_proteins),
        ):
            observed_raw = observed_records.get(transcript_id)
            observed = observed_raw.rstrip("*") if sequence_type == "protein" and observed_raw is not None else observed_raw
            present = observed is not None
            matches = present and observed == ours
            status = "PASS" if matches else ("MISSING_ID" if not present else "SEQUENCE_MISMATCH")
            comparison_rows.append(
                {
                    "transcript_id": transcript_id,
                    "sequence_type": sequence_type,
                    "python_length": len(ours),
                    "gffread_length": len(observed) if observed is not None else "",
                    "gffread_ID_present": str(present).lower(),
                    "sequence_match": str(bool(matches)).lower(),
                    "status": status,
                }
            )
            if not matches:
                discrepancies.append(f"{transcript_id}:{sequence_type}:{status}")
    if discrepancies:
        preview = ", ".join(discrepancies[:5])
        raise PrimaryAnnotationError(
            f"gffread comparison failed for {len(discrepancies)} selected sequences ({preview}); "
            "no output was published"
        )
    stdout_bytes, stdout_digest = _checksum(stdout_path)
    stderr_bytes, stderr_digest = _checksum(stderr_path)
    metadata = {
        "status": "PASS",
        "version": version,
        "comparison_annotation_scope": comparison_annotation_scope,
        "comparison_GFF3_file_name": selected_primary_gff_path.name,
        "comparison_GFF3_bytes": selected_gff_bytes,
        "comparison_GFF3_sha256": selected_gff_sha256,
        "source_full_GFF3_passed_to_gffread": False,
        "selected_transcripts_compared": len(selected),
        "gffread_CDS_records": len(gffread_cds),
        "gffread_protein_records": len(gffread_proteins),
        "exact_selected_CDS_ID_set": True,
        "exact_selected_protein_ID_set": True,
        "stdout_bytes": stdout_bytes,
        "stdout_sha256": stdout_digest,
        "stderr_bytes": stderr_bytes,
        "stderr_sha256": stderr_digest,
    }
    return comparison_rows, metadata


def _validate_inputs(genome: Path, gff: Path, output_dir: Path, sample_id: str) -> None:
    if not SAFE_SAMPLE_ID.fullmatch(sample_id):
        raise PrimaryAnnotationError(
            "sample_id must start with an alphanumeric character and contain only letters, "
            "numbers, periods, underscores, or hyphens"
        )
    for label, path in (("genome", genome), ("GFF3", gff)):
        if not path.is_file() or path.stat().st_size == 0:
            raise PrimaryAnnotationError(f"The {label} input is missing or empty: {path.name}")
    if genome.resolve() == gff.resolve():
        raise PrimaryAnnotationError("Genome and GFF3 inputs must be different files")
    if output_dir.exists() or output_dir.is_symlink():
        raise PrimaryAnnotationError(
            f"Output directory already exists; refusing overwrite: {output_dir.name}"
        )


def standardize_primary_annotation(
    genome_path: str | Path,
    gff_path: str | Path,
    output_dir: str | Path,
    sample_id: str,
    *,
    canonical_rules: Iterable[CanonicalRule] = (),
    transcript_features: Iterable[str] = ("mRNA", "transcript"),
    gene_features: Iterable[str] = ("gene", "pseudogene"),
    missing_phase_policy: str = "fail",
    invalid_coding_gene_policy: str = "fail",
    gene_as_transcript: bool = False,
    gffread: str = "auto",
    require_gffread: bool = False,
) -> StandardizationResult:
    """Create and atomically install one audited primary-annotation bundle."""
    genome = Path(genome_path).expanduser()
    gff = Path(gff_path).expanduser()
    output = Path(output_dir).expanduser()
    canonical = tuple(canonical_rules)
    transcript_types = frozenset(item.strip() for item in transcript_features if item.strip())
    gene_types = frozenset(item.strip() for item in gene_features if item.strip())
    if not transcript_types or not gene_types:
        raise PrimaryAnnotationError("At least one transcript and one gene feature type are required")
    if gene_as_transcript and transcript_types & gene_types:
        raise PrimaryAnnotationError(
            "gene-as-transcript mode requires disjoint transcript and gene feature types"
        )
    if missing_phase_policy not in {"fail", "zero"}:
        raise PrimaryAnnotationError("missing_phase_policy must be 'fail' or 'zero'")
    if invalid_coding_gene_policy not in {"fail", "omit"}:
        raise PrimaryAnnotationError("invalid_coding_gene_policy must be 'fail' or 'omit'")
    _validate_inputs(genome, gff, output, sample_id)
    executable = _resolve_executable(gffread)
    if require_gffread and executable is None:
        raise PrimaryAnnotationError("gffread is required but was not found")

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=str(output.parent))
    )
    try:
        input_records = []
        input_signatures: dict[Path, tuple[int, int, int, int]] = {}
        for role, path in (("chromosome_scope_genome", genome), ("chromosome_scope_GFF3", gff)):
            input_signatures[path] = _input_signature(path)
            size, digest = _checksum(path)
            input_records.append(
                {"role": role, "file_name": path.name, "bytes": size, "sha256": digest}
            )
        genome_records = load_genome(genome)
        annotation = parse_annotation(
            gff,
            {identifier: len(sequence) for identifier, sequence in genome_records.items()},
            transcript_types,
            gene_types,
            gene_as_transcript=gene_as_transcript,
        )
        if annotation.graph_mode == "gene_as_transcript" and canonical:
            raise PrimaryAnnotationError(
                "Canonical transcript-tag rules cannot be applied to synthesized "
                "gene-as-transcript models"
            )
        selected, transcript_audit, gene_audit, invalid_gene_count = select_primary_isoforms(
            annotation, genome_records, canonical, missing_phase_policy
        )
        if invalid_gene_count and invalid_coding_gene_policy == "fail":
            raise PrimaryAnnotationError(
                f"{invalid_gene_count} genes have CDS annotation but no valid coding transcript; "
                "review the source annotation or explicitly use --invalid-coding-gene-policy omit"
            )
        if not selected:
            raise PrimaryAnnotationError("No valid coding transcript remained after validation")

        prefix = sample_id
        protein_path = staging / f"{prefix}.protein.faa"
        cds_path = staging / f"{prefix}.cds.fa"
        primary_gff_path = staging / f"{prefix}.primary.gff3"
        primary_table_path = staging / f"{prefix}.primary_isoforms.tsv"
        transcript_audit_path = staging / f"{prefix}.transcript_audit.tsv"
        gene_audit_path = staging / f"{prefix}.gene_audit.tsv"
        coords_path = staging / f"{prefix}.coords.tsv"
        legacy_coords_path = staging / f"{prefix}.coords"
        comparison_path = staging / f"{prefix}.gffread_comparison.tsv"

        _write_fasta(
            protein_path,
            ((item.model.transcript_id, item.protein) for item in selected),
        )
        _write_fasta(cds_path, ((item.model.transcript_id, item.cds) for item in selected))
        _write_primary_gff(
            primary_gff_path,
            selected,
            annotation.gene_locations,
            annotation.graph_mode,
        )
        primary_rows = []
        coords_rows = []
        gene_selection_rules = {
            str(row["gene_id"]): str(row["selection_rule"]) for row in gene_audit
        }
        for item in selected:
            model = item.model
            primary_rows.append(
                {
                    "sample_id": sample_id,
                    "gene_id": model.gene_id,
                    "transcript_id": model.transcript_id,
                    "seqid": model.seqid,
                    "start": model.start,
                    "end": model.end,
                    "strand": model.strand,
                    "spliced_CDS_length": len(item.cds),
                    "protein_length": len(item.protein),
                    "canonical_rule": item.canonical_rule_label,
                    "selection_rule": gene_selection_rules[model.gene_id],
                    "QC_flags": ";".join(item.flags),
                }
            )
            coords_rows.append(
                {
                    "transcript_id": model.transcript_id,
                    "gene_id": model.gene_id,
                    "chromosome": model.seqid,
                    "start": model.start,
                    "end": model.end,
                    "strand": model.strand,
                }
            )
        _write_tsv(
            primary_table_path,
            primary_rows,
            (
                "sample_id", "gene_id", "transcript_id", "seqid", "start", "end",
                "strand", "spliced_CDS_length", "protein_length", "canonical_rule",
                "selection_rule", "QC_flags",
            ),
        )
        _write_tsv(
            transcript_audit_path,
            transcript_audit,
            (
                "gene_id", "transcript_id", "seqid", "start", "end", "strand",
                "raw_CDS_length", "validated_spliced_CDS_length", "protein_length",
                "canonical_rule", "validation_status", "QC_flags", "disposition",
                "gene_selection_rule",
            ),
        )
        _write_tsv(
            gene_audit_path,
            gene_audit,
            (
                "gene_id", "annotated_transcript_count", "transcripts_with_CDS_count",
                "valid_coding_transcript_count", "invalid_coding_transcript_count",
                "selected_transcript_id", "selection_rule", "status",
            ),
        )
        _write_tsv(
            coords_path,
            coords_rows,
            ("transcript_id", "gene_id", "chromosome", "start", "end", "strand"),
        )
        with legacy_coords_path.open("w", encoding="utf-8", newline="") as handle:
            for row in coords_rows:
                handle.write(
                    f"{row['transcript_id']}\t{row['chromosome']}\t{row['start']}\t"
                    f"{row['end']}\t{row['strand']}\n"
                )

        if executable is not None:
            gffread_work = staging / ".gffread_work"
            gffread_work.mkdir()
            comparison_scope = (
                "selected_primary_GFF3_top_level_mRNA_CDS_only_gene_as_transcript"
                if annotation.graph_mode == "gene_as_transcript"
                else "selected_primary_GFF3_gene_mRNA_CDS_only"
            )
            comparison_rows, comparison_metadata = compare_with_gffread(
                executable,
                genome,
                primary_gff_path,
                selected,
                gffread_work,
                comparison_scope,
            )
            shutil.rmtree(gffread_work)
            gffread_status = "PASS"
        else:
            skipped_status = (
                "NOT_RUN_EXPLICITLY_DISABLED"
                if gffread == "none"
                else "NOT_RUN_GFFREAD_NOT_AVAILABLE"
            )
            skipped_gff_bytes, skipped_gff_sha256 = _checksum(primary_gff_path)
            comparison_rows = [
                {
                    "transcript_id": "",
                    "sequence_type": "",
                    "python_length": "",
                    "gffread_length": "",
                    "gffread_ID_present": "",
                    "sequence_match": "",
                    "status": skipped_status,
                }
            ]
            comparison_metadata = {
                "status": skipped_status,
                "version": "",
                "comparison_annotation_scope": (
                    "selected_primary_GFF3_top_level_mRNA_CDS_only_gene_as_transcript"
                    if annotation.graph_mode == "gene_as_transcript"
                    else "selected_primary_GFF3_gene_mRNA_CDS_only"
                ),
                "comparison_GFF3_file_name": primary_gff_path.name,
                "comparison_GFF3_bytes": skipped_gff_bytes,
                "comparison_GFF3_sha256": skipped_gff_sha256,
                "source_full_GFF3_passed_to_gffread": False,
                "selected_transcripts_compared": 0,
                "exact_selected_CDS_ID_set": None,
                "exact_selected_protein_ID_set": None,
            }
            gffread_status = skipped_status
        _write_tsv(
            comparison_path,
            comparison_rows,
            (
                "transcript_id", "sequence_type", "python_length", "gffread_length",
                "gffread_ID_present", "sequence_match", "status",
            ),
        )

        manifest = {
            "schema_version": 1,
            "workflow": "primary_annotation_standardization",
            "workflow_version": WORKFLOW_VERSION,
            "status": "PASS",
            "publication_gate": (
                "PASS" if gffread_status == "PASS" else "BLOCKED_GFFREAD_NOT_RUN"
            ),
            "sample_id": sample_id,
            "execution": {"processes": 1, "worker_threads": 0},
            "inputs": input_records,
            "policy": {
                "transcript_feature_types": sorted(transcript_types),
                "gene_feature_types": sorted(gene_types),
                "canonical_tag_rules_in_priority_order": [rule.label for rule in canonical],
                "fallback_selection": (
                    "longest_valid_spliced_CDS;then_genomic_span;then_transcript_ID"
                ),
                "missing_phase_policy": missing_phase_policy,
                "invalid_coding_gene_policy": invalid_coding_gene_policy,
                "gene_as_transcript_requested": gene_as_transcript,
                "annotation_graph_mode": annotation.graph_mode,
                "internal_stop_policy": "reject_transcript",
                "terminal_stop_policy": "remove_from_protein_only",
                "ambiguous_codon_policy": "retain_as_X_and_flag",
            },
            "counts": {
                "chromosome_sequences": len(genome_records),
                "source_gene_IDs": len(annotation.gene_ids),
                "source_transcripts": len(annotation.transcripts),
                "source_transcript_rows": annotation.source_transcript_row_count,
                "synthesized_gene_as_transcripts": (
                    len(annotation.transcripts)
                    if annotation.graph_mode == "gene_as_transcript"
                    else 0
                ),
                "selected_genes": len(selected),
                "selected_transcripts": len(selected),
                "invalid_coding_genes": invalid_gene_count,
            },
            "gffread_comparison": comparison_metadata,
        }
        manifest_path = staging / "run_manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        for path, signature in input_signatures.items():
            if _input_signature(path) != signature:
                raise PrimaryAnnotationError(
                    f"Input {path.name} changed while it was being standardized; no output was published"
                )
        checksum_rows = []
        for path in sorted(staging.iterdir(), key=lambda item: item.name):
            if path.is_file() and path.name != "checksums.tsv":
                size, digest = _checksum(path)
                checksum_rows.append({"file": path.name, "bytes": size, "sha256": digest})
        _write_tsv(staging / "checksums.tsv", checksum_rows, ("file", "bytes", "sha256"))
        if output.exists() or output.is_symlink():
            raise PrimaryAnnotationError(
                f"Output directory appeared during the run; refusing overwrite: {output.name}"
            )
        os.replace(staging, output)
        return StandardizationResult(
            output_dir=output,
            source_gene_count=len(annotation.gene_ids),
            selected_gene_count=len(selected),
            invalid_coding_gene_count=invalid_gene_count,
            selected_transcript_count=len(selected),
            gffread_status=gffread_status,
        )
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
