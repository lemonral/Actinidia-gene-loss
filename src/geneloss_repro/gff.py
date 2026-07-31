"""GFF3 and FASTA readers used by the reusable workflow.

The legacy scripts selected the longest *genomic-span* mRNA and silently
dropped CDS features with a phase.  This module uses the spliced CDS length,
honours GFF3 phase, and writes a transparent isoform-selection table.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from .io_utils import SchemaError, natural_key, parse_int


@dataclass(frozen=True)
class GffFeature:
    sequence_id: str
    feature_type: str
    start: int
    end: int
    strand: str
    phase: str
    attributes: dict[str, str]
    line_number: int

    @property
    def identifier(self) -> str | None:
        return self.attributes.get("ID")

    @property
    def parents(self) -> list[str]:
        parent = self.attributes.get("Parent", "")
        return [item.strip() for item in parent.split(",") if item.strip()]


@dataclass(frozen=True)
class Transcript:
    gene_id: str
    transcript_id: str
    chrom: str
    start: int
    end: int
    strand: str
    cds_features: tuple[GffFeature, ...]

    @property
    def spliced_cds_length(self) -> int:
        return sum(feature.end - feature.start + 1 for feature in self.cds_features)

    @property
    def genomic_span(self) -> int:
        return self.end - self.start + 1


def parse_attributes(value: str) -> dict[str, str]:
    """Parse semicolon-delimited GFF3 attributes without assuming field order."""
    attributes: dict[str, str] = {}
    for item in value.strip().split(";"):
        if not item or "=" not in item:
            continue
        key, content = item.split("=", 1)
        attributes[key.strip()] = content.strip()
    return attributes


def iter_gff(path: str | Path) -> Iterator[GffFeature]:
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise SchemaError(f"{path}:{line_number}: expected 9 GFF3 columns, found {len(fields)}")
            start = parse_int(fields[3], "start", f"{path}:{line_number}")
            end = parse_int(fields[4], "end", f"{path}:{line_number}")
            if start < 1 or end < start:
                raise SchemaError(f"{path}:{line_number}: invalid interval {start}-{end}")
            strand = fields[6]
            if strand not in {"+", "-", ".", "?"}:
                raise SchemaError(f"{path}:{line_number}: invalid strand {strand!r}")
            yield GffFeature(
                sequence_id=fields[0],
                feature_type=fields[2],
                start=start,
                end=end,
                strand=strand,
                phase=fields[7],
                attributes=parse_attributes(fields[8]),
                line_number=line_number,
            )


def read_gene_catalog(path: str | Path, gene_feature: str = "gene") -> list[dict[str, object]]:
    """Read one unambiguous GFF gene record per gene ID for spatial denominators."""
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    for feature in iter_gff(path):
        if feature.feature_type != gene_feature:
            continue
        gene_id = feature.identifier
        if not gene_id:
            continue
        if gene_id in seen:
            raise SchemaError(f"{path}:{feature.line_number}: duplicate {gene_feature} ID {gene_id!r}")
        seen.add(gene_id)
        rows.append(
            {
                "target_gene": gene_id,
                "target_chromosome": feature.sequence_id,
                "target_start": feature.start,
                "target_end": feature.end,
                "target_midpoint": (feature.start + feature.end) / 2,
            }
        )
    if not rows:
        raise SchemaError(f"{path}: no {gene_feature!r} features with an ID attribute")
    return sorted(rows, key=lambda row: (natural_key(str(row["target_chromosome"])), row["target_start"], str(row["target_gene"])))


def collect_transcripts(
    path: str | Path,
    transcript_features: set[str] | None = None,
) -> tuple[list[Transcript], list[dict[str, str]]]:
    """Collect transcripts and their CDS records.

    Returns usable transcript objects and an audit list for malformed or
    incomplete records.  Multiple Parent values on a CDS are supported.
    """
    transcript_features = transcript_features or {"mRNA", "transcript"}
    transcript_metadata: dict[str, tuple[str, str, int, int, str]] = {}
    cds_by_transcript: dict[str, list[GffFeature]] = defaultdict(list)
    audit: list[dict[str, str]] = []

    for feature in iter_gff(path):
        if feature.feature_type in transcript_features:
            transcript_id = feature.identifier
            parents = feature.parents
            if not transcript_id or not parents:
                audit.append(
                    {
                        "record_type": "transcript",
                        "record_id": transcript_id or "",
                        "reason": "missing_ID_or_Parent",
                        "line_number": str(feature.line_number),
                    }
                )
                continue
            if len(parents) != 1:
                audit.append(
                    {
                        "record_type": "transcript",
                        "record_id": transcript_id,
                        "reason": "multiple_gene_parents_not_supported",
                        "line_number": str(feature.line_number),
                    }
                )
                continue
            if transcript_id in transcript_metadata:
                raise SchemaError(f"{path}:{feature.line_number}: duplicate transcript ID {transcript_id!r}")
            transcript_metadata[transcript_id] = (
                parents[0], feature.sequence_id, feature.start, feature.end, feature.strand
            )
        elif feature.feature_type == "CDS":
            parents = feature.parents
            if not parents:
                audit.append(
                    {
                        "record_type": "CDS",
                        "record_id": feature.identifier or "",
                        "reason": "missing_Parent",
                        "line_number": str(feature.line_number),
                    }
                )
                continue
            for transcript_id in parents:
                cds_by_transcript[transcript_id].append(feature)

    transcripts: list[Transcript] = []
    for transcript_id, (gene_id, chrom, start, end, strand) in transcript_metadata.items():
        cds = cds_by_transcript.get(transcript_id, [])
        if not cds:
            audit.append(
                {
                    "record_type": "transcript",
                    "record_id": transcript_id,
                    "reason": "no_CDS",
                    "line_number": "",
                }
            )
            continue
        if any(item.sequence_id != chrom or item.strand != strand for item in cds):
            audit.append(
                {
                    "record_type": "transcript",
                    "record_id": transcript_id,
                    "reason": "CDS_chromosome_or_strand_mismatch",
                    "line_number": "",
                }
            )
            continue
        transcripts.append(
            Transcript(gene_id, transcript_id, chrom, start, end, strand, tuple(cds))
        )

    return transcripts, audit


def select_longest_cds_isoform(transcripts: list[Transcript]) -> list[Transcript]:
    """Pick one deterministic primary isoform per gene.

    Tie breaking is: longer spliced CDS, then longer genomic span, then lexical
    transcript ID.  The explicit tie breaker prevents filesystem/GFF ordering
    from changing results between runs.
    """
    by_gene: dict[str, list[Transcript]] = defaultdict(list)
    for transcript in transcripts:
        by_gene[transcript.gene_id].append(transcript)
    chosen: list[Transcript] = []
    for gene_id, choices in by_gene.items():
        selected = sorted(
            choices,
            key=lambda item: (-item.spliced_cds_length, -item.genomic_span, item.transcript_id),
        )[0]
        chosen.append(selected)
    return sorted(chosen, key=lambda item: (natural_key(item.chrom), item.start, item.transcript_id))


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(identifier, uppercase_sequence)`` from a standard FASTA file."""
    identifier: str | None = None
    sequence_parts: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(sequence_parts).upper()
                identifier = line[1:].split()[0]
                if not identifier:
                    raise SchemaError(f"{path}:{line_number}: FASTA header has no identifier")
                sequence_parts = []
            elif identifier is None:
                raise SchemaError(f"{path}:{line_number}: sequence occurs before FASTA header")
            else:
                sequence_parts.append(line.replace(" ", ""))
    if identifier is not None:
        yield identifier, "".join(sequence_parts).upper()


def load_fasta(path: str | Path) -> dict[str, str]:
    result = dict(iter_fasta(path))
    if not result:
        raise SchemaError(f"{path}: no FASTA sequences found")
    return result
