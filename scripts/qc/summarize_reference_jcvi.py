#!/usr/bin/env python3
"""Summarize fixed-reference JCVI anchors across the 23 analysis units."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import statistics
import tempfile


HEADER = (
    "plot_order",
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "source_class",
    "anchors",
    "target_bed",
)


class SummaryError(RuntimeError):
    """Raised for incomplete or inconsistent JCVI evidence."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise SummaryError(f"{label} is missing or empty: {resolved}")
    return resolved


def beneath(root: Path, raw: str, label: str) -> Path:
    relative = Path(raw.strip())
    if not raw.strip() or relative.is_absolute() or ".." in relative.parts:
        raise SummaryError(f"{label} must be a safe relative path")
    resolved = (root / relative).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise SummaryError(f"{label} escapes data root") from error
    return require_file(resolved, label)


def read_bed(path: Path) -> tuple[dict[str, str], dict[str, int]]:
    gene_to_chromosome: dict[str, str] = {}
    chromosome_counts: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise SummaryError(f"{path}:{line_number}: BED requires at least four fields")
            chromosome, gene_id = fields[0], fields[3]
            if not chromosome or not gene_id or gene_id in gene_to_chromosome:
                raise SummaryError(f"{path}:{line_number}: duplicate or empty BED identifier")
            gene_to_chromosome[gene_id] = chromosome
            chromosome_counts[chromosome] = chromosome_counts.get(chromosome, 0) + 1
    if not gene_to_chromosome:
        raise SummaryError(f"empty BED: {path}")
    return gene_to_chromosome, chromosome_counts


def parse_anchors(
    path: Path,
    reference_genes: dict[str, str],
    target_genes: dict[str, str],
) -> dict[str, object]:
    reference_unique: set[str] = set()
    target_unique: set[str] = set()
    reference_to_targets: dict[str, set[str]] = {}
    chromosome_reference: dict[str, set[str]] = {}
    chromosome_target: dict[str, set[str]] = {}
    chromosome_rows: dict[str, int] = {}
    anchor_rows = 0
    blocks = 0
    block_has_rows = False
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            stripped = raw.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if block_has_rows:
                    blocks += 1
                    block_has_rows = False
                continue
            fields = stripped.split("\t")
            if len(fields) < 2:
                raise SummaryError(f"{path}:{line_number}: malformed anchor row")
            reference_gene, target_gene = fields[:2]
            if reference_gene not in reference_genes:
                raise SummaryError(
                    f"{path}:{line_number}: reference ID missing from BED: {reference_gene}"
                )
            if target_gene not in target_genes:
                raise SummaryError(
                    f"{path}:{line_number}: target ID missing from BED: {target_gene}"
                )
            block_has_rows = True
            anchor_rows += 1
            reference_unique.add(reference_gene)
            target_unique.add(target_gene)
            reference_to_targets.setdefault(reference_gene, set()).add(target_gene)
            reference_chromosome = reference_genes[reference_gene]
            chromosome_reference.setdefault(reference_chromosome, set()).add(reference_gene)
            chromosome_target.setdefault(reference_chromosome, set()).add(target_gene)
            chromosome_rows[reference_chromosome] = (
                chromosome_rows.get(reference_chromosome, 0) + 1
            )
    if block_has_rows:
        blocks += 1
    if not anchor_rows:
        raise SummaryError(f"no anchor rows: {path}")
    copy_depths = [len(targets) for targets in reference_to_targets.values()]
    return {
        "anchor_blocks": blocks,
        "anchor_pair_rows": anchor_rows,
        "reference_unique": reference_unique,
        "target_unique": target_unique,
        "mean_target_genes_per_anchored_reference_gene": statistics.fmean(copy_depths),
        "median_target_genes_per_anchored_reference_gene": statistics.median(copy_depths),
        "reference_genes_by_chromosome": chromosome_reference,
        "target_genes_by_reference_chromosome": chromosome_target,
        "anchor_rows_by_reference_chromosome": chromosome_rows,
    }


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--reference-bed", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        root = args.data_root.expanduser().resolve()
        sources = require_file(args.sources, "source manifest")
        reference_bed = require_file(args.reference_bed, "reference BED")
        output = args.output_dir.expanduser().resolve()
        if output.exists():
            raise SummaryError(f"refusing existing output directory: {output}")
        reference_genes, reference_chromosome_counts = read_bed(reference_bed)
        reference_chromosomes = sorted(reference_chromosome_counts)
        if len(reference_chromosomes) != 24:
            raise SummaryError(
                f"expected 24 reference assembly sequences, found {len(reference_chromosomes)}"
            )
        with sources.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != HEADER:
                raise SummaryError(f"source header must be exactly {HEADER}")
            source_rows = list(reader)
        if len(source_rows) != 23:
            raise SummaryError(f"expected 23 analysis units, found {len(source_rows)}")

        temporary = Path(
            tempfile.mkdtemp(prefix=f".{output.name}.tmp.", dir=output.parent)
        )
        summary_rows: list[list[object]] = []
        chromosome_rows: list[list[object]] = []
        input_checksums: list[dict[str, object]] = []
        observed_units: set[str] = set()
        reference_coverages: list[float] = []
        target_coverages: list[float] = []

        for row in source_rows:
            unit = row["assembly_unit_id"].strip()
            if not unit or unit in observed_units:
                raise SummaryError(f"duplicate or empty unit: {unit!r}")
            observed_units.add(unit)
            anchors = beneath(root, row["anchors"], f"{unit} anchors")
            target_bed = beneath(root, row["target_bed"], f"{unit} BED")
            target_genes, target_chromosome_counts = read_bed(target_bed)
            metrics = parse_anchors(anchors, reference_genes, target_genes)
            reference_unique = metrics["reference_unique"]
            target_unique = metrics["target_unique"]
            assert isinstance(reference_unique, set)
            assert isinstance(target_unique, set)
            anchor_rows = int(metrics["anchor_pair_rows"])
            reference_coverage = 100.0 * len(reference_unique) / len(reference_genes)
            target_coverage = 100.0 * len(target_unique) / len(target_genes)
            reference_coverages.append(reference_coverage)
            target_coverages.append(target_coverage)
            summary_rows.append(
                [
                    int(row["plot_order"]),
                    unit,
                    row["biological_species"],
                    row["haplotype_or_subgenome"],
                    row["source_class"],
                    int(metrics["anchor_blocks"]),
                    anchor_rows,
                    len(reference_genes),
                    len(reference_unique),
                    reference_coverage,
                    len(target_genes),
                    len(target_unique),
                    target_coverage,
                    len(
                        {
                            reference_genes[gene]
                            for gene in reference_unique
                        }
                    ),
                    len(
                        {
                            target_genes[gene]
                            for gene in target_unique
                        }
                    ),
                    float(metrics["mean_target_genes_per_anchored_reference_gene"]),
                    float(metrics["median_target_genes_per_anchored_reference_gene"]),
                ]
            )
            reference_by_chromosome = metrics["reference_genes_by_chromosome"]
            target_by_chromosome = metrics["target_genes_by_reference_chromosome"]
            rows_by_chromosome = metrics["anchor_rows_by_reference_chromosome"]
            assert isinstance(reference_by_chromosome, dict)
            assert isinstance(target_by_chromosome, dict)
            assert isinstance(rows_by_chromosome, dict)
            for chromosome in reference_chromosomes:
                anchored_reference = len(reference_by_chromosome.get(chromosome, set()))
                denominator = reference_chromosome_counts[chromosome]
                chromosome_rows.append(
                    [
                        int(row["plot_order"]),
                        unit,
                        chromosome,
                        denominator,
                        anchored_reference,
                        100.0 * anchored_reference / denominator,
                        len(target_by_chromosome.get(chromosome, set())),
                        int(rows_by_chromosome.get(chromosome, 0)),
                    ]
                )
            for role, path in (("anchors", anchors), ("target_bed", target_bed)):
                input_checksums.append(
                    {
                        "assembly_unit_id": unit,
                        "role": role,
                        "basename": path.name,
                        "bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )

        summary_rows.sort(key=lambda values: int(values[0]))
        chromosome_rows.sort(key=lambda values: (int(values[0]), str(values[2])))
        summary_path = temporary / "clematoclethra_actinidia_jcvi_summary.tsv"
        chromosome_path = temporary / "clematoclethra_actinidia_jcvi_chromosome_summary.tsv"
        write_tsv(
            summary_path,
            [
                "plot_order",
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "source_class",
                "anchor_blocks",
                "anchor_pair_rows",
                "reference_bed_genes",
                "reference_anchored_genes",
                "reference_anchored_gene_percent",
                "target_bed_genes",
                "target_anchored_genes",
                "target_anchored_gene_percent",
                "reference_chromosomes_with_anchors",
                "target_chromosomes_with_anchors",
                "mean_target_genes_per_anchored_reference_gene",
                "median_target_genes_per_anchored_reference_gene",
            ],
            summary_rows,
        )
        write_tsv(
            chromosome_path,
            [
                "plot_order",
                "assembly_unit_id",
                "reference_chromosome",
                "reference_chromosome_genes",
                "anchored_reference_genes",
                "anchored_reference_gene_percent",
                "anchored_target_genes",
                "anchor_pair_rows",
            ],
            chromosome_rows,
        )
        validation = {
            "schema_version": 1,
            "status": "PASS_CLEMATOCLETHRA_ACTINIDIA_JCVI_REFERENCE_SUPPORT",
            "unit_count": len(summary_rows),
            "reference_gene_count": len(reference_genes),
            "reference_chromosome_count": len(reference_chromosomes),
            "reference_bed_sha256": sha256(reference_bed),
            "metric_definition": (
                "unique genes occurring in raw JCVI collinear anchors divided by "
                "the corresponding JCVI BED gene count"
            ),
            "minimum_reference_anchored_gene_percent": min(reference_coverages),
            "maximum_reference_anchored_gene_percent": max(reference_coverages),
            "minimum_target_anchored_gene_percent": min(target_coverages),
            "maximum_target_anchored_gene_percent": max(target_coverages),
            "inputs": input_checksums,
            "outputs": [],
        }
        validation_path = temporary / "validation.json"
        validation_path.write_text(
            json.dumps(validation, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        for path in (summary_path, chromosome_path):
            validation["outputs"].append(
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
        validation_path.write_text(
            json.dumps(validation, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, output)
        print(f"PASS: {output}")
        return 0
    except (OSError, UnicodeError, ValueError, SummaryError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
