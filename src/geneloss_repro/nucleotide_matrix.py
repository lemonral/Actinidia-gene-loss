"""Build validated chromosome nucleotide-homology matrices from paired PAFs.

The matrix represents one target assembly against one frozen reference.  It
combines target-to-reference and reference-to-target primary minimap2 rows,
uses interval unions for both chromosome denominators, and retains a complete
query-by-reference Cartesian product so an absent alignment remains explicit.
"""

from __future__ import annotations

import gzip
import math
from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .chromosome_assignment import NUCLEOTIDE_COLUMNS
from .io_utils import natural_key


class NucleotideMatrixError(RuntimeError):
    """Raised when FASTA or PAF evidence violates the frozen contract."""


@dataclass(frozen=True)
class PafRecord:
    """One retained, role-normalized target/reference PAF alignment."""

    target_chromosome: str
    target_length: int
    target_start: int
    target_end: int
    reference_chromosome: str
    reference_length: int
    reference_start: int
    reference_end: int
    orientation: str
    matching_bases: int
    alignment_block_length: int
    divergence: Decimal


@dataclass(frozen=True)
class PafAudit:
    """Counts proving how the fixed PAF filters were applied."""

    input_rows: int
    retained_rows: int
    nonprimary_rows: int
    low_mapq_rows: int
    short_block_rows: int
    high_divergence_rows: int


def fasta_lengths(path: str | Path) -> dict[str, int]:
    """Return strict unique FASTA identifiers and sequence lengths."""

    source = Path(path)
    opener = gzip.open if source.suffix.lower() == ".gz" else open
    lengths: dict[str, int] = {}
    current: str | None = None
    try:
        with opener(source, "rt", encoding="utf-8", errors="strict") as handle:
            for line_number, line in enumerate(handle, start=1):
                if line.startswith(">"):
                    identifier = line[1:].strip().split()[0] if line[1:].strip() else ""
                    if not identifier:
                        raise NucleotideMatrixError(
                            f"{source.name}:{line_number}: empty FASTA identifier"
                        )
                    if identifier in lengths:
                        raise NucleotideMatrixError(
                            f"{source.name}:{line_number}: duplicate FASTA identifier {identifier!r}"
                        )
                    current = identifier
                    lengths[current] = 0
                    continue
                if current is None:
                    if line.strip():
                        raise NucleotideMatrixError(
                            f"{source.name}:{line_number}: sequence before first FASTA header"
                        )
                    continue
                sequence = "".join(line.split())
                if sequence:
                    lengths[current] += len(sequence)
    except (OSError, UnicodeError) as error:
        raise NucleotideMatrixError(f"Cannot read FASTA {source.name}: {error}") from error
    if not lengths or any(length <= 0 for length in lengths.values()):
        raise NucleotideMatrixError(f"{source.name}: every FASTA record must be non-empty")
    return lengths


def _exact_int(value: str, *, label: str, line: int, path: Path) -> int:
    if not value.isdigit():
        raise NucleotideMatrixError(f"{path.name}:{line}: {label} is not an integer")
    return int(value)


def _tags(fields: Sequence[str], *, line: int, path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) != 3 or not parts[0]:
            raise NucleotideMatrixError(f"{path.name}:{line}: malformed PAF tag {field!r}")
        key = parts[0]
        if key in result:
            raise NucleotideMatrixError(f"{path.name}:{line}: duplicate PAF tag {key!r}")
        result[key] = field
    return result


def read_role_normalized_paf(
    path: str | Path,
    *,
    query_role: str,
    target_lengths: Mapping[str, int],
    reference_lengths: Mapping[str, int],
    minimum_mapq: int = 20,
    minimum_alignment_block_bp: int = 10_000,
    maximum_de: Decimal = Decimal("0.15"),
) -> tuple[list[PafRecord], PafAudit]:
    """Read one PAF and normalize coordinates to target/reference roles."""

    if query_role not in {"target", "reference"}:
        raise NucleotideMatrixError("query_role must be 'target' or 'reference'")
    source = Path(path)
    records: list[PafRecord] = []
    counts = defaultdict(int)
    try:
        with source.open("rt", encoding="utf-8", errors="strict") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                if not raw_line.strip():
                    continue
                counts["input"] += 1
                fields = raw_line.rstrip("\n\r").split("\t")
                if len(fields) < 12:
                    raise NucleotideMatrixError(
                        f"{source.name}:{line_number}: PAF row has fewer than 12 fields"
                    )
                q_name, t_name = fields[0], fields[5]
                q_length = _exact_int(fields[1], label="query length", line=line_number, path=source)
                q_start = _exact_int(fields[2], label="query start", line=line_number, path=source)
                q_end = _exact_int(fields[3], label="query end", line=line_number, path=source)
                strand = fields[4]
                t_length = _exact_int(fields[6], label="target length", line=line_number, path=source)
                t_start = _exact_int(fields[7], label="target start", line=line_number, path=source)
                t_end = _exact_int(fields[8], label="target end", line=line_number, path=source)
                matching = _exact_int(fields[9], label="matching bases", line=line_number, path=source)
                block = _exact_int(fields[10], label="alignment block length", line=line_number, path=source)
                mapq = _exact_int(fields[11], label="MAPQ", line=line_number, path=source)
                if strand not in {"+", "-"}:
                    raise NucleotideMatrixError(
                        f"{source.name}:{line_number}: invalid PAF strand {strand!r}"
                    )
                if not (0 <= q_start < q_end <= q_length and 0 <= t_start < t_end <= t_length):
                    raise NucleotideMatrixError(
                        f"{source.name}:{line_number}: PAF coordinates exceed declared lengths"
                    )
                if matching > block or block <= 0:
                    raise NucleotideMatrixError(
                        f"{source.name}:{line_number}: invalid matching/block lengths"
                    )
                tags = _tags(fields[12:], line=line_number, path=source)
                required = {"tp", "de", "cg", "cs"}
                if not required.issubset(tags):
                    missing = ",".join(sorted(required.difference(tags)))
                    raise NucleotideMatrixError(
                        f"{source.name}:{line_number}: missing required PAF tags: {missing}"
                    )
                if tags["tp"] != "tp:A:P":
                    counts["nonprimary"] += 1
                    continue
                try:
                    divergence = Decimal(tags["de"].split(":", 2)[2])
                except Exception as error:
                    raise NucleotideMatrixError(
                        f"{source.name}:{line_number}: invalid de tag"
                    ) from error
                if not divergence.is_finite() or not Decimal(0) <= divergence <= Decimal(1):
                    raise NucleotideMatrixError(
                        f"{source.name}:{line_number}: de must lie in [0,1]"
                    )
                if mapq < minimum_mapq:
                    counts["low_mapq"] += 1
                    continue
                if block < minimum_alignment_block_bp:
                    counts["short_block"] += 1
                    continue
                if divergence > maximum_de:
                    counts["high_divergence"] += 1
                    continue

                if query_role == "target":
                    target_name, target_length, target_start, target_end = (
                        q_name,
                        q_length,
                        q_start,
                        q_end,
                    )
                    reference_name, reference_length, reference_start, reference_end = (
                        t_name,
                        t_length,
                        t_start,
                        t_end,
                    )
                else:
                    target_name, target_length, target_start, target_end = (
                        t_name,
                        t_length,
                        t_start,
                        t_end,
                    )
                    reference_name, reference_length, reference_start, reference_end = (
                        q_name,
                        q_length,
                        q_start,
                        q_end,
                    )
                if target_lengths.get(target_name) != target_length:
                    raise NucleotideMatrixError(
                        f"{source.name}:{line_number}: target chromosome/length is not in FASTA"
                    )
                if reference_lengths.get(reference_name) != reference_length:
                    raise NucleotideMatrixError(
                        f"{source.name}:{line_number}: reference chromosome/length is not in FASTA"
                    )
                records.append(
                    PafRecord(
                        target_chromosome=target_name,
                        target_length=target_length,
                        target_start=target_start,
                        target_end=target_end,
                        reference_chromosome=reference_name,
                        reference_length=reference_length,
                        reference_start=reference_start,
                        reference_end=reference_end,
                        orientation=strand,
                        matching_bases=matching,
                        alignment_block_length=block,
                        divergence=divergence,
                    )
                )
                counts["retained"] += 1
    except (OSError, UnicodeError) as error:
        raise NucleotideMatrixError(f"Cannot read PAF {source.name}: {error}") from error
    if counts["input"] == 0:
        raise NucleotideMatrixError(f"{source.name}: PAF contains no rows")
    return records, PafAudit(
        input_rows=counts["input"],
        retained_rows=counts["retained"],
        nonprimary_rows=counts["nonprimary"],
        low_mapq_rows=counts["low_mapq"],
        short_block_rows=counts["short_block"],
        high_divergence_rows=counts["high_divergence"],
    )


def _union_length(intervals: Iterable[tuple[int, int]]) -> int:
    ordered = sorted(intervals)
    if not ordered:
        return 0
    total = 0
    start, end = ordered[0]
    for next_start, next_end in ordered[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end
    return total + end - start


def _fraction(numerator: int, denominator: int) -> Decimal:
    return Decimal(numerator) / Decimal(denominator)


def _harmonic(first: Decimal, second: Decimal) -> Decimal:
    if first == 0 or second == 0:
        return Decimal(0)
    return Decimal(2) * first * second / (first + second)


def _number(value: Decimal) -> str:
    if value == 0:
        return "0"
    if value == 1:
        return "1"
    return format(value, ".15g")


def build_nucleotide_rows(
    *,
    target_lengths: Mapping[str, int],
    reference_lengths: Mapping[str, int],
    canonical_by_reference: Mapping[str, str],
    forward_records: Sequence[PafRecord],
    reverse_records: Sequence[PafRecord],
    orientation_minimum_fraction: Decimal = Decimal("0.80"),
    orientation_minimum_matching_bp: int = 1_000_000,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Create exact matrix rows plus a companion orientation-support audit."""

    if set(reference_lengths) != set(canonical_by_reference):
        raise NucleotideMatrixError("Reference FASTA IDs do not equal the frozen map IDs")
    grouped: dict[tuple[str, str], list[PafRecord]] = defaultdict(list)
    for record in (*forward_records, *reverse_records):
        if target_lengths.get(record.target_chromosome) != record.target_length:
            raise NucleotideMatrixError("PAF target length changed after parsing")
        if reference_lengths.get(record.reference_chromosome) != record.reference_length:
            raise NucleotideMatrixError("PAF reference length changed after parsing")
        grouped[(record.target_chromosome, record.reference_chromosome)].append(record)

    rows: list[dict[str, str]] = []
    orientation_rows: list[dict[str, str]] = []
    target_ids = sorted(target_lengths, key=natural_key)
    reference_ids = sorted(reference_lengths, key=natural_key)
    for target in target_ids:
        for reference in reference_ids:
            evidence = grouped.get((target, reference), [])
            target_covered = _union_length(
                (record.target_start, record.target_end) for record in evidence
            )
            reference_covered = _union_length(
                (record.reference_start, record.reference_end) for record in evidence
            )
            matching = sum(record.matching_bases for record in evidence)
            aligned = sum(record.alignment_block_length for record in evidence)
            plus = sum(
                record.matching_bases for record in evidence if record.orientation == "+"
            )
            minus = matching - plus
            if not evidence:
                divergence = Decimal(1)
                orientation = "none"
                dominant_fraction = Decimal(0)
                dominant = "none"
            else:
                divergence = sum(
                    record.divergence * record.alignment_block_length for record in evidence
                ) / Decimal(aligned)
                dominant = "+" if plus > minus else "-" if minus > plus else "mixed"
                dominant_fraction = Decimal(max(plus, minus)) / Decimal(matching)
                orientation = (
                    dominant
                    if dominant in {"+", "-"}
                    and dominant_fraction >= orientation_minimum_fraction
                    and matching >= orientation_minimum_matching_bp
                    else "mixed"
                )
            target_coverage = _fraction(target_covered, target_lengths[target])
            reference_coverage = _fraction(reference_covered, reference_lengths[reference])
            reciprocal = min(target_coverage, reference_coverage)
            score = _harmonic(target_coverage, reference_coverage) * (Decimal(1) - divergence)
            rows.append(
                {
                    "query_chromosome": target,
                    "reference_chromosome": reference,
                    "canonical_chromosome": canonical_by_reference[reference],
                    "score": _number(score),
                    "query_covered_bp": str(target_covered),
                    "query_length_bp": str(target_lengths[target]),
                    "query_coverage": _number(target_coverage),
                    "reference_covered_bp": str(reference_covered),
                    "reference_length_bp": str(reference_lengths[reference]),
                    "reference_coverage": _number(reference_coverage),
                    "reciprocal_coverage": _number(reciprocal),
                    "matching_bases": str(matching),
                    "weighted_divergence": _number(divergence),
                    "orientation": orientation,
                }
            )
            orientation_rows.append(
                {
                    "query_chromosome": target,
                    "reference_chromosome": reference,
                    "plus_matching_bases": str(plus),
                    "minus_matching_bases": str(minus),
                    "total_matching_bases": str(matching),
                    "dominant_orientation": dominant,
                    "dominant_fraction": _number(dominant_fraction),
                    "automatic_orientation_gate": str(orientation in {"+", "-"}).lower(),
                }
            )
    expected = len(target_lengths) * len(reference_lengths)
    if len(rows) != expected or tuple(rows[0]) != NUCLEOTIDE_COLUMNS:
        raise NucleotideMatrixError("Internal matrix schema or Cartesian closure failure")
    return rows, orientation_rows


__all__ = [
    "NucleotideMatrixError",
    "PafAudit",
    "PafRecord",
    "build_nucleotide_rows",
    "fasta_lengths",
    "read_role_normalized_paf",
]
