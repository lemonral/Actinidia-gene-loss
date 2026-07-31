#!/usr/bin/env python3
"""Dependency-free helpers for reconstructing JCVI gene-depth coverage.

The legacy ``jcvi.compara.synteny depth`` metric is based on gene indices,
not covered genome bases.  For each anchor block, JCVI uses the minimum and
maximum BED indices as a half-open interval ``[minimum, maximum)``.  This
module deliberately preserves that endpoint convention so that new and
historical comparisons use the same denominator and coverage definition.

All parsers fail closed on malformed coordinates, duplicate BED accessions,
missing anchor accessions, empty inputs, and blocks spanning sequence IDs.
They do not encode any project-specific species or sample expectations.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


class ReconstructionError(RuntimeError):
    """Raised when JCVI BED or anchor inputs violate the data contract."""


@dataclass(frozen=True)
class BedIndex:
    """A deterministically ordered BED and its accession lookup table."""

    path: Path
    rows: tuple[tuple[str, int, int, str], ...]
    order: dict[str, tuple[int, str]]


def natural_key(value: str) -> tuple[object, ...]:
    """Return a natural-sort key suitable for chromosome/scaffold names."""

    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", value)
    )


def load_bed(path: Path) -> BedIndex:
    """Load BED features and assign deterministic zero-based gene indices."""

    rows: list[tuple[str, int, int, str]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                fields = raw.split()
            if len(fields) < 4:
                raise ReconstructionError(
                    f"{path}:{line_number}: expected at least four BED fields"
                )
            seqid, start_text, end_text, accession = fields[:4]
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as exc:
                raise ReconstructionError(
                    f"{path}:{line_number}: non-integer BED coordinate"
                ) from exc
            if not seqid or not accession or start < 0 or end < start:
                raise ReconstructionError(f"{path}:{line_number}: invalid BED feature")
            rows.append((seqid, start, end, accession))
    if not rows:
        raise ReconstructionError(f"BED contains no features: {path}")

    rows.sort(key=lambda row: (natural_key(row[0]), row[1], row[3]))
    order: dict[str, tuple[int, str]] = {}
    for index, row in enumerate(rows):
        accession = row[3]
        if accession in order:
            raise ReconstructionError(
                f"BED contains duplicate accession {accession!r}: {path}"
            )
        order[accession] = (index, row[0])
    return BedIndex(path=path, rows=tuple(rows), order=order)


def load_anchor_blocks(path: Path) -> tuple[list[list[tuple[str, str]]], int]:
    """Load raw JCVI anchors, returning blocks and the total pair-row count."""

    blocks: list[list[tuple[str, str]]] = []
    current: list[tuple[str, str]] = []
    pair_rows = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            if raw.startswith("#"):
                if current:
                    blocks.append(current)
                    current = []
                continue
            fields = raw.split()
            if len(fields) < 2:
                raise ReconstructionError(
                    f"{path}:{line_number}: expected reference and query accession columns"
                )
            current.append((fields[0], fields[1]))
            pair_rows += 1
    if current:
        blocks.append(current)
    if not blocks or pair_rows == 0:
        raise ReconstructionError(f"Anchor file contains no blocks: {path}")
    return blocks, pair_rows


def _merge_half_open(intervals: Iterable[tuple[int, int]]) -> tuple[int, int]:
    ordered = sorted((start, end) for start, end in intervals if end > start)
    if not ordered:
        return 0, 0
    covered = 0
    merged_count = 0
    start, end = ordered[0]
    for next_start, next_end in ordered[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            covered += end - start
            merged_count += 1
            start, end = next_start, next_end
    covered += end - start
    merged_count += 1
    return covered, merged_count


def summarize_depth(
    blocks: list[list[tuple[str, str]]], bed: BedIndex, side: int
) -> dict[str, int | float]:
    """Reconstruct non-zero gene-depth counts for one side of raw anchors.

    ``side`` is zero for reference accessions and one for query accessions.
    The right endpoint of every block is intentionally excluded to reproduce
    JCVI's legacy range-depth behavior.
    """

    if side not in {0, 1}:
        raise ReconstructionError(f"Anchor side must be 0 or 1, found {side}")
    intervals: list[tuple[int, int]] = []
    missing: set[str] = set()
    cross_sequence_blocks = 0
    for block in blocks:
        if not block:
            raise ReconstructionError("Anchor block is empty")
        located: list[tuple[int, str]] = []
        for pair in block:
            record = bed.order.get(pair[side])
            if record is None:
                missing.add(pair[side])
            else:
                located.append(record)
        if not located:
            continue
        seqids = {seqid for _, seqid in located}
        if len(seqids) != 1:
            cross_sequence_blocks += 1
        indices = [index for index, _ in located]
        intervals.append((min(indices), max(indices)))

    if missing:
        examples = ", ".join(sorted(missing)[:5])
        raise ReconstructionError(
            f"{len(missing)} anchor accessions are absent from {bed.path}; "
            f"examples: {examples}"
        )
    if cross_sequence_blocks:
        raise ReconstructionError(
            f"{cross_sequence_blocks} anchor blocks cross BED sequence IDs in {bed.path}"
        )

    nonzero, merged_count = _merge_half_open(intervals)
    total = len(bed.rows)
    if nonzero < 0 or nonzero > total:
        raise ReconstructionError(
            f"Invalid non-zero depth count {nonzero}/{total} for {bed.path}"
        )
    zero = total - nonzero
    return {
        "total": total,
        "zero": zero,
        "nonzero": nonzero,
        "coverage": nonzero * 100.0 / total,
        "merged_intervals": merged_count,
        "cross_chromosome": cross_sequence_blocks,
    }
