"""Chromosome relabelling and orientation harmonization with sequence closure."""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, TextIO

from .annotation import build_spliced_cds, reverse_complement, translate_standard
from .gff import collect_transcripts


class HarmonizationError(RuntimeError):
    """Raised when a chromosome-scope bundle cannot be transformed exactly."""


@dataclass(frozen=True)
class ChromosomeAction:
    source_chromosome: str
    final_chromosome: str
    orientation: str
    source_length: int


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", value)
    )


def _open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    return path.open("r", encoding="utf-8", errors="strict")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current: str | None = None
    pieces: list[str] = []

    def finish() -> None:
        nonlocal current, pieces
        if current is None:
            return
        sequence = "".join(pieces).upper()
        if not sequence:
            raise HarmonizationError(f"{path.name}: empty FASTA record {current!r}")
        records[current] = sequence

    with _open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                finish()
                header = raw[1:].strip()
                current = header.split()[0] if header else ""
                if not current or current in records:
                    raise HarmonizationError(
                        f"{path.name}:{line_number}: empty or duplicate FASTA identifier"
                    )
                pieces = []
            elif raw.strip():
                if current is None:
                    raise HarmonizationError(
                        f"{path.name}:{line_number}: sequence before FASTA header"
                    )
                pieces.append("".join(raw.split()))
    finish()
    if not records:
        raise HarmonizationError(f"FASTA contains no records: {path}")
    return records


def write_fasta(path: Path, records: Mapping[str, str], width: int = 60) -> None:
    def emit(handle: TextIO) -> None:
        for identifier in sorted(records, key=natural_key):
            sequence = records[identifier]
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")

    if path.suffix.lower() == ".gz":
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as handle:
                    emit(handle)
    else:
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            emit(handle)


def build_actions(
    genome: Mapping[str, str],
    final_by_source: Mapping[str, str],
    orientation_by_source: Mapping[str, str],
) -> dict[str, ChromosomeAction]:
    if set(genome) != set(final_by_source) or set(genome) != set(orientation_by_source):
        raise HarmonizationError("Genome, final map, and orientation map scopes differ")
    if len(set(final_by_source.values())) != len(final_by_source):
        raise HarmonizationError("Final chromosome labels are not a bijection")
    actions: dict[str, ChromosomeAction] = {}
    for source in sorted(genome, key=natural_key):
        orientation = orientation_by_source[source]
        if orientation not in {"+", "-"}:
            raise HarmonizationError(
                f"{source}: orientation is {orientation!r}; manual review is required"
            )
        actions[source] = ChromosomeAction(
            source, final_by_source[source], orientation, len(genome[source])
        )
    return actions


def transform_genome(
    genome: Mapping[str, str], actions: Mapping[str, ChromosomeAction]
) -> dict[str, str]:
    if set(genome) != set(actions):
        raise HarmonizationError("Genome and chromosome action scopes differ")
    transformed: dict[str, str] = {}
    for source, action in actions.items():
        sequence = genome[source]
        if len(sequence) != action.source_length:
            raise HarmonizationError(f"{source}: source length changed")
        output_sequence = sequence if action.orientation == "+" else reverse_complement(sequence)
        if action.final_chromosome in transformed:
            raise HarmonizationError("Final FASTA chromosome label is duplicated")
        transformed[action.final_chromosome] = output_sequence.upper()
    return transformed


def transform_gff(
    source: Path,
    destination: Path,
    actions: Mapping[str, ChromosomeAction],
) -> dict[str, int]:
    feature_rows = 0
    reversed_rows = 0
    sequence_directives = 0
    with source.open("r", encoding="utf-8", errors="strict") as input_handle, destination.open(
        "w", encoding="utf-8", newline=""
    ) as output_handle:
        for line_number, raw in enumerate(input_handle, 1):
            if raw.startswith("##FASTA"):
                raise HarmonizationError("Embedded GFF3 FASTA is not allowed")
            if raw.startswith("##sequence-region"):
                fields = raw.rstrip("\n\r").split()
                if len(fields) != 4 or fields[1] not in actions:
                    raise HarmonizationError(
                        f"{source.name}:{line_number}: invalid sequence-region directive"
                    )
                action = actions[fields[1]]
                if fields[2] != "1" or fields[3] != str(action.source_length):
                    raise HarmonizationError(
                        f"{source.name}:{line_number}: sequence-region length mismatch"
                    )
                output_handle.write(
                    f"##sequence-region {action.final_chromosome} 1 {action.source_length}\n"
                )
                sequence_directives += 1
                continue
            if not raw.strip() or raw.startswith("#"):
                output_handle.write(raw)
                continue
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) != 9:
                raise HarmonizationError(
                    f"{source.name}:{line_number}: expected nine GFF3 columns"
                )
            action = actions.get(fields[0])
            if action is None:
                raise HarmonizationError(
                    f"{source.name}:{line_number}: GFF3 sequence is outside chromosome scope"
                )
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise HarmonizationError(
                    f"{source.name}:{line_number}: non-integer GFF3 coordinate"
                ) from error
            if start < 1 or end < start or end > action.source_length:
                raise HarmonizationError(
                    f"{source.name}:{line_number}: GFF3 interval is outside chromosome"
                )
            fields[0] = action.final_chromosome
            if action.orientation == "-":
                fields[3] = str(action.source_length - end + 1)
                fields[4] = str(action.source_length - start + 1)
                if fields[6] == "+":
                    fields[6] = "-"
                elif fields[6] == "-":
                    fields[6] = "+"
                elif fields[6] not in {".", "?"}:
                    raise HarmonizationError(
                        f"{source.name}:{line_number}: invalid GFF3 strand"
                    )
                reversed_rows += 1
            output_handle.write("\t".join(fields) + "\n")
            feature_rows += 1
    if feature_rows == 0:
        raise HarmonizationError("GFF3 contains no feature rows")
    return {
        "feature_rows": feature_rows,
        "reversed_feature_rows": reversed_rows,
        "sequence_region_directives": sequence_directives,
    }


def validate_sequence_closure(
    *,
    source_genome: Mapping[str, str],
    transformed_genome: Mapping[str, str],
    actions: Mapping[str, ChromosomeAction],
) -> None:
    if set(transformed_genome) != {action.final_chromosome for action in actions.values()}:
        raise HarmonizationError("Transformed FASTA scope differs from final chromosome map")
    for source, action in actions.items():
        expected = (
            source_genome[source]
            if action.orientation == "+"
            else reverse_complement(source_genome[source])
        ).upper()
        if transformed_genome[action.final_chromosome].upper() != expected:
            raise HarmonizationError(f"{source}: transformed FASTA sequence mismatch")


def validate_cds_protein_closure(
    *,
    source_genome: Mapping[str, str],
    source_gff: Path,
    transformed_genome: Mapping[str, str],
    transformed_gff: Path,
    expected_cds: Mapping[str, str],
    expected_proteins: Mapping[str, str],
) -> dict[str, int]:
    source_transcripts, source_audit = collect_transcripts(source_gff)
    target_transcripts, target_audit = collect_transcripts(transformed_gff)
    if source_audit or target_audit:
        raise HarmonizationError(
            "Source or transformed GFF3 contains transcript-model audit failures"
        )
    source_by_id = {item.transcript_id: item for item in source_transcripts}
    target_by_id = {item.transcript_id: item for item in target_transcripts}
    if set(source_by_id) != set(target_by_id):
        raise HarmonizationError("Source and transformed transcript ID sets differ")
    if set(source_by_id) != set(expected_cds) or set(source_by_id) != set(expected_proteins):
        raise HarmonizationError("GFF3, CDS FASTA, and protein FASTA ID sets differ")
    source_genome_dict = dict(source_genome)
    transformed_genome_dict = dict(transformed_genome)
    for transcript_id in sorted(source_by_id):
        source_cds = build_spliced_cds(source_by_id[transcript_id], source_genome_dict)
        transformed_cds = build_spliced_cds(
            target_by_id[transcript_id], transformed_genome_dict
        )
        if source_cds != transformed_cds or source_cds != expected_cds[transcript_id].upper():
            raise HarmonizationError(
                f"{transcript_id}: CDS sequence changed during harmonization"
            )
        protein = translate_standard(source_cds).rstrip("*")
        if protein != expected_proteins[transcript_id].upper().rstrip("*"):
            raise HarmonizationError(
                f"{transcript_id}: translated protein differs from frozen protein FASTA"
            )
    return {
        "transcripts": len(source_by_id),
        "exact_cds_matches": len(source_by_id),
        "exact_protein_matches": len(source_by_id),
    }
