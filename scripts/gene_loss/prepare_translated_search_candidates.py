#!/usr/bin/env python3
"""Prepare exact SynOrths-flanked candidates for translated genome search.

The historical workflow considered every unobserved reference gene within 20
reference-gene positions of a SynOrths anchor.  This script reproduces that
candidate universe, but declares a locus callable only when it has bilateral
anchors within the same 20-gene window, both anchors map unambiguously to the
same target chromosome, the inferred target interval is not excessive, and
the interval contains enough unambiguous A/C/G/T sequence.
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
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import TextIO


class CandidateError(RuntimeError):
    pass


@dataclass(frozen=True)
class Coordinate:
    gene: str
    chromosome: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class TargetAnchor:
    chromosome: str
    start: int
    end: int


REQUIRED_MANIFEST = ("unit", "target_genome", "synorth_pairs", "output_dir")
CANDIDATE_COLUMNS = (
    "unit",
    "reference_gene",
    "reference_chromosome",
    "reference_gene_index_1based",
    "has_reference_cds",
    "left_anchor",
    "left_anchor_distance_genes",
    "right_anchor",
    "right_anchor_distance_genes",
    "target_chromosome",
    "target_interval_start_1based",
    "target_interval_end_1based",
    "target_interval_span_bp",
    "target_interval_acgt_fraction",
    "callable",
    "callability_reason",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path) -> dict[str, object]:
    source = path.resolve()
    if not source.is_file() or source.stat().st_size <= 0:
        raise CandidateError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def open_text(path: Path) -> TextIO:
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
            raise CandidateError(f"{path.name}: empty FASTA record {current!r}")
        records[current] = sequence

    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                finish()
                header = raw[1:].strip()
                current = header.split()[0] if header else ""
                if not current or current in records:
                    raise CandidateError(
                        f"{path.name}:{line_number}: empty or duplicate FASTA identifier"
                    )
                pieces = []
            elif raw.strip():
                if current is None:
                    raise CandidateError(f"{path.name}:{line_number}: sequence before header")
                pieces.append("".join(raw.split()))
    finish()
    if not records:
        raise CandidateError(f"{path.name}: FASTA contains no records")
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for identifier in sorted(records):
            handle.write(f">{identifier}\n")
            sequence = records[identifier]
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def read_coords(path: Path) -> tuple[dict[str, Coordinate], dict[str, list[str]]]:
    by_gene: dict[str, Coordinate] = {}
    by_chromosome: dict[str, list[Coordinate]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if len(fields) != 5:
                raise CandidateError(f"{path.name}:{line_number}: expected five columns")
            gene, chromosome, raw_start, raw_end, strand = fields
            if not gene or gene in by_gene or not chromosome or strand not in {"+", "-"}:
                raise CandidateError(f"{path.name}:{line_number}: invalid or duplicate coordinate")
            try:
                start, end = int(raw_start), int(raw_end)
            except ValueError as error:
                raise CandidateError(f"{path.name}:{line_number}: non-integer coordinate") from error
            if start < 1 or end < start:
                raise CandidateError(f"{path.name}:{line_number}: invalid interval")
            item = Coordinate(gene, chromosome, start, end, strand)
            by_gene[gene] = item
            by_chromosome[chromosome].append(item)
    ordered = {
        chromosome: [item.gene for item in sorted(items, key=lambda x: (x.start, x.end, x.gene))]
        for chromosome, items in by_chromosome.items()
    }
    if not by_gene:
        raise CandidateError(f"{path.name}: no coordinates")
    return by_gene, ordered


def read_synorth(
    path: Path,
    reference: dict[str, Coordinate],
    reference_column_1based: int = 5,
    target_gene_column_1based: int = 1,
) -> dict[str, list[TargetAnchor]]:
    if reference_column_1based < 1 or target_gene_column_1based < 1:
        raise CandidateError("SynOrths reference/target columns must be >=1")
    anchors: dict[str, list[TargetAnchor]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if len(fields) < 8:
                raise CandidateError(f"{path.name}:{line_number}: fewer than eight columns")
            if len(fields) < reference_column_1based:
                raise CandidateError(
                    f"{path.name}:{line_number}: missing audited reference column"
                )
            reference_gene = fields[reference_column_1based - 1]
            if reference_gene not in reference:
                raise CandidateError(f"{path.name}:{line_number}: unknown reference gene")
            reference_chromosome_column = reference_column_1based
            if (
                len(fields) <= reference_chromosome_column
                or fields[reference_chromosome_column] != reference[reference_gene].chromosome
            ):
                raise CandidateError(f"{path.name}:{line_number}: reference chromosome mismatch")
            target_chromosome_column = target_gene_column_1based
            try:
                query_chr = fields[target_chromosome_column]
                start = int(fields[target_chromosome_column + 1])
                end = int(fields[target_chromosome_column + 2])
            except ValueError as error:
                raise CandidateError(f"{path.name}:{line_number}: non-integer target coordinate") from error
            except IndexError as error:
                raise CandidateError(f"{path.name}:{line_number}: missing target coordinate fields") from error
            if not query_chr or start < 1 or end < start:
                raise CandidateError(f"{path.name}:{line_number}: invalid target anchor")
            anchors[reference_gene].append(TargetAnchor(query_chr, start, end))
    if not anchors:
        raise CandidateError(f"{path.name}: no SynOrths anchors")
    return anchors


def nearest_anchor(
    genes: list[str], observed: set[str], index: int, step: int, window: int
) -> tuple[str, int] | None:
    for distance in range(1, window + 1):
        candidate = index + step * distance
        if candidate < 0 or candidate >= len(genes):
            return None
        if genes[candidate] in observed:
            return genes[candidate], distance
    return None


def anchor_scope(items: list[TargetAnchor]) -> tuple[str, int, int] | None:
    chromosomes = {item.chromosome for item in items}
    if len(chromosomes) != 1:
        return None
    return next(iter(chromosomes)), min(item.start for item in items), max(item.end for item in items)


def build_unit(
    *,
    unit: str,
    genome_path: Path,
    synorth_path: Path,
    reference_coords_path: Path,
    reference_cds_path: Path,
    output_dir: Path,
    window: int,
    padding_bp: int,
    maximum_locus_bp: int,
    minimum_acgt_fraction: float,
    synorth_reference_column_1based: int,
    synorth_target_gene_column_1based: int,
) -> dict[str, object]:
    if output_dir.exists():
        raise CandidateError(f"refusing to overwrite output: {output_dir}")
    reference, ordered = read_coords(reference_coords_path)
    reference_cds = read_fasta(reference_cds_path)
    if not set(reference_cds).issubset(reference):
        raise CandidateError("reference CDS contains IDs absent from reference coordinates")
    genome = read_fasta(genome_path)
    anchors = read_synorth(
        synorth_path,
        reference,
        reference_column_1based=synorth_reference_column_1based,
        target_gene_column_1based=synorth_target_gene_column_1based,
    )
    observed = set(anchors)
    rows: list[dict[str, str]] = []
    query_records: dict[str, str] = {}

    for chromosome in sorted(ordered):
        genes = ordered[chromosome]
        legacy_indices: set[int] = set()
        for index, gene in enumerate(genes):
            if gene not in observed:
                continue
            legacy_indices.update(
                item for item in range(max(0, index - window), min(len(genes), index + window + 1))
                if genes[item] not in observed
            )
        for index in sorted(legacy_indices):
            gene = genes[index]
            left = nearest_anchor(genes, observed, index, -1, window)
            right = nearest_anchor(genes, observed, index, +1, window)
            reason = "callable"
            target_chromosome = ""
            interval_start = interval_end = span = 0
            acgt_fraction: float | None = None
            if gene not in reference_cds:
                reason = "missing_reference_cds"
            elif left is None or right is None:
                reason = "missing_bilateral_anchor"
            else:
                left_scope = anchor_scope(anchors[left[0]])
                right_scope = anchor_scope(anchors[right[0]])
                if left_scope is None or right_scope is None:
                    reason = "ambiguous_anchor_target"
                elif left_scope[0] != right_scope[0]:
                    reason = "flanks_map_to_different_target_chromosomes"
                else:
                    target_chromosome = left_scope[0]
                    if target_chromosome not in genome:
                        raise CandidateError(
                            f"{unit}: SynOrths chromosome {target_chromosome!r} absent from genome"
                        )
                    interval_start = max(1, min(left_scope[1], right_scope[1]) - padding_bp)
                    interval_end = min(
                        len(genome[target_chromosome]),
                        max(left_scope[2], right_scope[2]) + padding_bp,
                    )
                    span = interval_end - interval_start + 1
                    sequence = genome[target_chromosome][interval_start - 1 : interval_end]
                    acgt_fraction = sum(base in "ACGT" for base in sequence) / len(sequence)
                    if span > maximum_locus_bp:
                        reason = "excessive_target_interval"
                    elif acgt_fraction < minimum_acgt_fraction:
                        reason = "low_unambiguous_sequence_fraction"
            callable_locus = reason == "callable"
            if gene in reference_cds:
                query_records[gene] = reference_cds[gene]
            rows.append(
                {
                    "unit": unit,
                    "reference_gene": gene,
                    "reference_chromosome": chromosome,
                    "reference_gene_index_1based": str(index + 1),
                    "has_reference_cds": str(gene in reference_cds).lower(),
                    "left_anchor": "" if left is None else left[0],
                    "left_anchor_distance_genes": "" if left is None else str(left[1]),
                    "right_anchor": "" if right is None else right[0],
                    "right_anchor_distance_genes": "" if right is None else str(right[1]),
                    "target_chromosome": target_chromosome,
                    "target_interval_start_1based": "" if not interval_start else str(interval_start),
                    "target_interval_end_1based": "" if not interval_end else str(interval_end),
                    "target_interval_span_bp": "" if not span else str(span),
                    "target_interval_acgt_fraction": (
                        "" if acgt_fraction is None else f"{acgt_fraction:.12g}"
                    ),
                    "callable": str(callable_locus).lower(),
                    "callability_reason": reason,
                }
            )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.staging.", dir=output_dir.parent))
    try:
        candidate_table = staging / "candidates.tsv"
        with candidate_table.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=CANDIDATE_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        query_fasta = staging / "candidate_reference_cds.fasta"
        write_fasta(query_fasta, query_records)
        reason_counts: dict[str, int] = defaultdict(int)
        for row in rows:
            reason_counts[row["callability_reason"]] += 1
        report = {
            "schema_version": 1,
            "workflow": "synorth_flanked_translated_search_candidates",
            "status": "PASS",
            "unit": unit,
            "created_at_utc": utc_now(),
            "parameters": {
                "reference_gene_window": window,
                "target_interval_padding_bp": padding_bp,
                "maximum_target_interval_bp": maximum_locus_bp,
                "minimum_target_interval_acgt_fraction": minimum_acgt_fraction,
                "synorth_reference_column_1based": synorth_reference_column_1based,
                "synorth_target_gene_column_1based": synorth_target_gene_column_1based,
                "candidate_definition": "unobserved reference gene within the window of any SynOrths anchor",
                "callable_definition": "bilateral unambiguous same-chromosome anchors plus interval quality gates",
            },
            "inputs": {
                "target_genome": binding(genome_path),
                "synorth_pairs": binding(synorth_path),
                "reference_coords": binding(reference_coords_path),
                "reference_cds": binding(reference_cds_path),
            },
            "metrics": {
                "reference_coordinate_genes": len(reference),
                "reference_cds_records": len(reference_cds),
                "target_chromosomes": len(genome),
                "observed_reference_anchors": len(observed),
                "candidate_rows": len(rows),
                "candidate_query_cds_records": len(query_records),
                "callable_rows": reason_counts.get("callable", 0),
                "callability_reasons": dict(sorted(reason_counts.items())),
            },
            "outputs": {
                "candidates": binding(candidate_table),
                "candidate_reference_cds": binding(query_fasta),
            },
        }
        report_path = staging / "run_manifest.json"
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        checksums = staging / "checksums.tsv"
        with checksums.open("w", encoding="utf-8", newline="") as handle:
            handle.write("file\tbytes\tsha256\n")
            for path in sorted(staging.iterdir()):
                if path.is_file() and path.name != checksums.name:
                    item = binding(path)
                    handle.write(f"{path.name}\t{item['bytes']}\t{item['sha256']}\n")
        os.rename(staging, output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return report


def resolve(root: Path, value: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise CandidateError(f"unsafe data-root-relative path: {value!r}")
    # Validate the lexical data-root-relative location before following a
    # registered compatibility symlink.  Legacy assets are intentionally
    # exposed below data_root through symlinks whose real targets may live in
    # the frozen manuscript data tree; resolving first incorrectly rejects
    # those otherwise safe manifest entries.
    result = (root / relative).absolute()
    if not result.is_relative_to(root) or not result.is_file() or result.stat().st_size <= 0:
        raise CandidateError(f"missing or unsafe input: {value!r}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--reference-coords", required=True, type=Path)
    parser.add_argument("--reference-cds", required=True, type=Path)
    parser.add_argument("--window", type=int, default=20)
    parser.add_argument("--padding-bp", type=int, default=10000)
    parser.add_argument("--maximum-locus-bp", type=int, default=5000000)
    parser.add_argument("--minimum-acgt-fraction", type=float, default=0.80)
    parser.add_argument(
        "--synorth-reference-column",
        type=int,
        default=5,
        help=(
            "1-based reference-gene column in SynOrths pairs. New-unit pairs use 5; "
            "the exact-bound legacy pairs use audited column 1."
        ),
    )
    parser.add_argument(
        "--synorth-target-gene-column",
        type=int,
        default=1,
        help=(
            "1-based target-genome gene column in SynOrths pairs. New-unit pairs use 1; "
            "the exact-bound legacy pairs use audited column 5."
        ),
    )
    args = parser.parse_args()
    try:
        root = args.data_root.resolve()
        if (
            args.window < 1
            or args.padding_bp < 0
            or args.maximum_locus_bp < 1
            or args.synorth_reference_column < 1
            or args.synorth_target_gene_column < 1
        ):
            raise CandidateError("invalid positive integer parameter")
        if not 0 <= args.minimum_acgt_fraction <= 1:
            raise CandidateError("minimum ACGT fraction must lie in [0,1]")
        with args.manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != REQUIRED_MANIFEST:
                raise CandidateError("manifest columns differ from the exact schema")
            rows = list(reader)
        if not rows or len({row["unit"] for row in rows}) != len(rows):
            raise CandidateError("manifest is empty or has duplicate units")
        reference_coords = args.reference_coords.resolve()
        reference_cds = args.reference_cds.resolve()
        summaries = []
        for row in rows:
            output = root / row["output_dir"]
            if not output.resolve().is_relative_to(root):
                raise CandidateError("output escapes data root")
            summaries.append(
                build_unit(
                    unit=row["unit"],
                    genome_path=resolve(root, row["target_genome"]),
                    synorth_path=resolve(root, row["synorth_pairs"]),
                    reference_coords_path=reference_coords,
                    reference_cds_path=reference_cds,
                    output_dir=output,
                    window=args.window,
                    padding_bp=args.padding_bp,
                    maximum_locus_bp=args.maximum_locus_bp,
                    minimum_acgt_fraction=args.minimum_acgt_fraction,
                    synorth_reference_column_1based=args.synorth_reference_column,
                    synorth_target_gene_column_1based=args.synorth_target_gene_column,
                )
            )
        print(f"PASS\t{len(summaries)} units")
        return 0
    except (OSError, UnicodeError, csv.Error, json.JSONDecodeError, CandidateError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
