"""Fail-closed helpers for bidirectional JCVI chromosome-homology matrices."""

from __future__ import annotations

import re
from dataclasses import dataclass
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Mapping


getcontext().prec = 50


class JcviMatrixError(RuntimeError):
    """Raised when BED, GFF, protein, or raw-anchor evidence is inconsistent."""


@dataclass(frozen=True)
class BedCatalog:
    chromosomes: tuple[str, ...]
    gene_to_chromosome: Mapping[str, str]
    eligible_by_chromosome: Mapping[str, int]
    rows_by_gene: Mapping[str, tuple[str, int, int, str]]


@dataclass(frozen=True)
class AnchorAudit:
    pair_rows: int
    unique_pairs: int
    block_headers: int


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", value)
    )


def read_bed(path: Path) -> BedCatalog:
    gene_to_chromosome: dict[str, str] = {}
    rows_by_gene: dict[str, tuple[str, int, int, str]] = {}
    counts: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 4:
                fields = raw.split()
            if len(fields) < 4:
                raise JcviMatrixError(
                    f"{path.name}:{line_number}: expected at least four BED columns"
                )
            chromosome, start_text, end_text, gene = fields[:4]
            strand = fields[5] if len(fields) >= 6 else "."
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as error:
                raise JcviMatrixError(
                    f"{path.name}:{line_number}: non-integer BED coordinate"
                ) from error
            if (
                not chromosome
                or not gene
                or start < 0
                or end <= start
                or strand not in {"+", "-", "."}
            ):
                raise JcviMatrixError(f"{path.name}:{line_number}: invalid BED row")
            if gene in gene_to_chromosome:
                raise JcviMatrixError(
                    f"{path.name}:{line_number}: duplicate BED accession {gene!r}"
                )
            gene_to_chromosome[gene] = chromosome
            rows_by_gene[gene] = (chromosome, start, end, strand)
            counts[chromosome] = counts.get(chromosome, 0) + 1
    if not gene_to_chromosome:
        raise JcviMatrixError(f"BED contains no features: {path}")
    chromosomes = tuple(sorted(counts, key=natural_key))
    return BedCatalog(chromosomes, gene_to_chromosome, counts, rows_by_gene)


def read_fasta_ids(path: Path) -> set[str]:
    identifiers: set[str] = set()
    sequence_seen = False
    current: str | None = None
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if current is not None and not sequence_seen:
                    raise JcviMatrixError(f"{path.name}: empty FASTA record {current!r}")
                header = raw[1:].strip()
                current = header.split()[0] if header else ""
                if not current or current in identifiers:
                    raise JcviMatrixError(
                        f"{path.name}:{line_number}: empty or duplicate FASTA identifier"
                    )
                identifiers.add(current)
                sequence_seen = False
            elif raw.strip():
                if current is None:
                    raise JcviMatrixError(
                        f"{path.name}:{line_number}: sequence precedes FASTA header"
                    )
                sequence_seen = True
    if current is not None and not sequence_seen:
        raise JcviMatrixError(f"{path.name}: empty FASTA record {current!r}")
    if not identifiers:
        raise JcviMatrixError(f"Protein FASTA contains no records: {path}")
    return identifiers


def require_bed_protein_identity(bed: BedCatalog, protein: Path, *, label: str) -> None:
    protein_ids = read_fasta_ids(protein)
    bed_ids = set(bed.gene_to_chromosome)
    if protein_ids != bed_ids:
        missing = sorted(protein_ids.difference(bed_ids))[:5]
        extra = sorted(bed_ids.difference(protein_ids))[:5]
        raise JcviMatrixError(
            f"{label} BED/protein IDs differ: protein_only={missing}; bed_only={extra}"
        )


def require_bed_gff_identity(bed: BedCatalog, gff: Path) -> None:
    observed: dict[str, tuple[str, int, int, str]] = {}
    with gff.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) != 9:
                raise JcviMatrixError(f"{gff.name}:{line_number}: expected nine GFF3 columns")
            chromosome, _, feature, start_text, end_text, _, strand, _, attributes = fields
            if feature not in {"mRNA", "transcript"}:
                continue
            values: dict[str, str] = {}
            for item in attributes.split(";"):
                if not item:
                    continue
                if "=" not in item:
                    raise JcviMatrixError(
                        f"{gff.name}:{line_number}: malformed GFF3 attributes"
                    )
                key, value = item.split("=", 1)
                if key in values or not key or not value:
                    raise JcviMatrixError(
                        f"{gff.name}:{line_number}: invalid GFF3 attributes"
                    )
                values[key] = value
            transcript = values.get("ID")
            if transcript not in bed.gene_to_chromosome:
                continue
            if transcript in observed:
                raise JcviMatrixError(
                    f"{gff.name}:{line_number}: duplicate selected transcript {transcript!r}"
                )
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as error:
                raise JcviMatrixError(
                    f"{gff.name}:{line_number}: non-integer GFF3 coordinate"
                ) from error
            observed[transcript] = (chromosome, start - 1, end, strand)
    if set(observed) != set(bed.gene_to_chromosome):
        missing = sorted(set(bed.gene_to_chromosome).difference(observed))[:5]
        raise JcviMatrixError(f"Target GFF3 lacks BED transcripts: {missing}")
    mismatches = [gene for gene, row in bed.rows_by_gene.items() if observed[gene] != row]
    if mismatches:
        raise JcviMatrixError(
            f"Target BED coordinates/strands differ from GFF3: {mismatches[:5]}"
        )


def relabel_reference_bed_from_canonical_truth(
    bed: BedCatalog, canonical_by_reference: Mapping[str, str]
) -> tuple[BedCatalog, dict[str, str]]:
    """Map canonical GFF/BED labels to frozen haplome-specific reference IDs.

    Some Hongyang publisher GFF3 files use ``Chr01`` while the matching genome
    uses ``Chr01A`` or ``Chr01P``.  Relabelling is allowed only when the BED
    chromosome set is *exactly* the unique canonical-label set in the frozen
    truth registry; no prefix stripping or fuzzy name matching is performed.
    """

    inverse: dict[str, str] = {}
    for reference, canonical in canonical_by_reference.items():
        if canonical in inverse:
            raise JcviMatrixError("Frozen canonical chromosome labels are not unique")
        inverse[canonical] = reference
    if set(bed.chromosomes) != set(inverse):
        raise JcviMatrixError(
            "Reference BED chromosome set is neither frozen reference IDs nor the exact "
            "canonical-label set"
        )
    gene_to_chromosome = {
        gene: inverse[chromosome] for gene, chromosome in bed.gene_to_chromosome.items()
    }
    eligible = {inverse[chromosome]: count for chromosome, count in bed.eligible_by_chromosome.items()}
    rows = {
        gene: (inverse[chromosome], start, end, strand)
        for gene, (chromosome, start, end, strand) in bed.rows_by_gene.items()
    }
    relabelled = BedCatalog(
        tuple(sorted(eligible, key=natural_key)), gene_to_chromosome, eligible, rows
    )
    return relabelled, {canonical: inverse[canonical] for canonical in sorted(inverse, key=natural_key)}


def read_normalized_anchor_pairs(
    path: Path,
    *,
    first_bed: BedCatalog,
    second_bed: BedCatalog,
    first_role: str,
) -> tuple[set[tuple[str, str]], AnchorAudit]:
    """Return role-normalized ``(target_gene, reference_gene)`` raw-anchor pairs."""

    if first_role not in {"target", "reference"}:
        raise JcviMatrixError("first_role must be target or reference")
    pairs: set[tuple[str, str]] = set()
    pair_rows = 0
    block_headers = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            if raw.startswith("#"):
                block_headers += 1
                continue
            fields = raw.split()
            if len(fields) < 2:
                raise JcviMatrixError(
                    f"{path.name}:{line_number}: expected two anchor accessions"
                )
            first, second = fields[:2]
            if first not in first_bed.gene_to_chromosome:
                raise JcviMatrixError(
                    f"{path.name}:{line_number}: first accession absent from its BED"
                )
            if second not in second_bed.gene_to_chromosome:
                raise JcviMatrixError(
                    f"{path.name}:{line_number}: second accession absent from its BED"
                )
            pair = (first, second) if first_role == "target" else (second, first)
            pair_rows += 1
            if pair in pairs:
                raise JcviMatrixError(
                    f"{path.name}:{line_number}: duplicate raw anchor pair {pair!r}"
                )
            pairs.add(pair)
    if not pairs:
        raise JcviMatrixError(f"Raw anchor file contains no pairs: {path}")
    return pairs, AnchorAudit(pair_rows, len(pairs), block_headers)


def _decimal_text(value: Decimal) -> str:
    if value == 0:
        return "0"
    return format(value, ".15g")


def _harmonic(first: Decimal, second: Decimal) -> Decimal:
    if first == 0 or second == 0:
        return Decimal(0)
    return Decimal(2) * first * second / (first + second)


def build_jcvi_rows(
    *,
    target_bed: BedCatalog,
    reference_bed: BedCatalog,
    canonical_by_reference: Mapping[str, str],
    forward_pairs: set[tuple[str, str]],
    reverse_pairs: set[tuple[str, str]],
) -> tuple[list[dict[str, str]], dict[str, int]]:
    if set(reference_bed.chromosomes) != set(canonical_by_reference):
        raise JcviMatrixError("Reference BED chromosome set differs from truth registry")
    pairs = forward_pairs.union(reverse_pairs)
    by_cell: dict[tuple[str, str], set[tuple[str, str]]] = {}
    for target_gene, reference_gene in pairs:
        try:
            target_chromosome = target_bed.gene_to_chromosome[target_gene]
            reference_chromosome = reference_bed.gene_to_chromosome[reference_gene]
        except KeyError as error:
            raise JcviMatrixError("Normalized anchor pair is absent from a BED") from error
        by_cell.setdefault((target_chromosome, reference_chromosome), set()).add(
            (target_gene, reference_gene)
        )

    rows: list[dict[str, str]] = []
    for target_chromosome in target_bed.chromosomes:
        query_eligible = target_bed.eligible_by_chromosome[target_chromosome]
        for reference_chromosome in reference_bed.chromosomes:
            reference_eligible = reference_bed.eligible_by_chromosome[reference_chromosome]
            cell_pairs = by_cell.get((target_chromosome, reference_chromosome), set())
            query_anchored = len({pair[0] for pair in cell_pairs})
            reference_anchored = len({pair[1] for pair in cell_pairs})
            query_coverage = Decimal(query_anchored) / Decimal(query_eligible)
            reference_coverage = Decimal(reference_anchored) / Decimal(reference_eligible)
            score = _harmonic(query_coverage, reference_coverage)
            rows.append(
                {
                    "query_chromosome": target_chromosome,
                    "reference_chromosome": reference_chromosome,
                    "canonical_chromosome": canonical_by_reference[reference_chromosome],
                    "score": _decimal_text(score),
                    "query_anchored_genes": str(query_anchored),
                    "query_eligible_genes": str(query_eligible),
                    "query_gene_coverage": _decimal_text(query_coverage),
                    "reference_anchored_genes": str(reference_anchored),
                    "reference_eligible_genes": str(reference_eligible),
                    "reference_gene_coverage": _decimal_text(reference_coverage),
                    "unique_anchor_pairs": str(len(cell_pairs)),
                }
            )
    return rows, {
        "forward_unique_pairs": len(forward_pairs),
        "reverse_unique_pairs": len(reverse_pairs),
        "bidirectional_union_pairs": len(pairs),
        "bidirectional_intersection_pairs": len(forward_pairs.intersection(reverse_pairs)),
    }
