#!/usr/bin/env python3
"""Localize conservative legacy deletions between exact SynOrths flanks.

Only historical ``deleted`` calls are used.  A spatial coordinate is emitted
when the deleted reference gene has bilateral accepted reference anchors
within the frozen 20-gene window and both anchors map to one target
chromosome.  The coordinate is the expected interval midpoint, not an
observed deleted fragment or a pseudogene hit.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


class SpatialInputError(RuntimeError):
    pass


MANIFEST_COLUMNS = (
    "assembly_unit_id", "legacy_sample", "biological_species", "haplotype_or_subgenome",
    "assembly_scope", "synorth_audit", "synorth_pairs", "deleted_genes", "genome", "gff",
)


@dataclass(frozen=True)
class Anchor:
    chromosome: str
    start: int
    end: int


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.resolve().open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    source = path.resolve()
    if not source.is_file() or (not allow_empty and source.stat().st_size == 0):
        raise SpatialInputError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def resolve(root: Path, value: str, *, allow_empty: bool = False) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise SpatialInputError(f"unsafe data-root-relative path: {value!r}")
    path = (root / relative).absolute()
    if not path.is_relative_to(root) or not path.is_file() or (
        not allow_empty and path.stat().st_size == 0
    ):
        raise SpatialInputError(f"missing or unsafe input: {value!r}")
    return path


def read_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise SpatialInputError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise SpatialInputError(f"{path}: JSON root is not an object")
    return value


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MANIFEST_COLUMNS:
            raise SpatialInputError("legacy spatial manifest columns differ from exact schema")
        rows = list(reader)
    units = [row["assembly_unit_id"] for row in rows]
    if not rows or any(not value for row in rows for value in row.values()) or len(units) != len(set(units)):
        raise SpatialInputError("legacy spatial manifest is empty or has empty/duplicate units")
    return rows


def read_id_list(path: Path) -> list[str]:
    values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(values) != len(set(values)):
        raise SpatialInputError(f"{path.name}: duplicate identifier rows")
    return values


def read_reference_order(path: Path) -> tuple[dict[str, tuple[str, int]], dict[str, list[str]]]:
    records: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    seen: set[str] = set()
    with path.resolve().open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if len(fields) != 5:
                raise SpatialInputError(f"{path.name}:{line_number}: expected five coordinate fields")
            gene, chromosome = fields[0], fields[1]
            try:
                start, end = int(fields[2]), int(fields[3])
            except ValueError as error:
                raise SpatialInputError(f"{path.name}:{line_number}: invalid coordinate") from error
            if not gene or gene in seen or not chromosome or start < 1 or end < start:
                raise SpatialInputError(f"{path.name}:{line_number}: invalid/duplicate reference row")
            seen.add(gene)
            records[chromosome].append((start, end, gene))
    order: dict[str, list[str]] = {}
    lookup: dict[str, tuple[str, int]] = {}
    for chromosome, items in records.items():
        genes = [gene for _, _, gene in sorted(items)]
        order[chromosome] = genes
        for index, gene in enumerate(genes):
            lookup[gene] = (chromosome, index)
    if not lookup:
        raise SpatialInputError(f"{path.name}: no reference coordinates")
    return lookup, order


def open_text(path: Path):
    return gzip.open(path.resolve(), "rt", encoding="utf-8") if path.suffix == ".gz" else path.resolve().open(
        "r", encoding="utf-8"
    )


def read_genome(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = ""
    pieces: list[str] = []
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if current:
                    records[current] = "".join(pieces).upper()
                current = raw[1:].strip().split()[0]
                if not current or current in records:
                    raise SpatialInputError(f"{path.name}:{line_number}: empty/duplicate FASTA ID")
                pieces = []
            elif raw.strip():
                if not current:
                    raise SpatialInputError(f"{path.name}:{line_number}: sequence before header")
                pieces.append("".join(raw.split()))
    if current:
        records[current] = "".join(pieces).upper()
    if not records or any(not sequence for sequence in records.values()):
        raise SpatialInputError(f"{path.name}: empty genome or record")
    return records


def read_anchors(path: Path, reference_column: int, query_column: int) -> dict[str, list[Anchor]]:
    anchors: dict[str, list[Anchor]] = defaultdict(list)
    seen: set[tuple[str, ...]] = set()
    with path.resolve().open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if len(fields) < 8:
                raise SpatialInputError(f"{path.name}:{line_number}: fewer than eight fields")
            key = tuple(fields)
            if key in seen:
                raise SpatialInputError(f"{path.name}:{line_number}: duplicate pair row")
            seen.add(key)
            reference_gene = fields[reference_column - 1]
            query_gene = fields[query_column - 1]
            # The frozen legacy SynOrths schema places query chromosome/start/end
            # immediately after the audited query-gene column.
            try:
                chromosome = fields[query_column]
                start, end = int(fields[query_column + 1]), int(fields[query_column + 2])
            except (IndexError, ValueError) as error:
                raise SpatialInputError(f"{path.name}:{line_number}: invalid query anchor fields") from error
            start, end = min(start, end), max(start, end)
            if not reference_gene or not query_gene or not chromosome or start < 1:
                raise SpatialInputError(f"{path.name}:{line_number}: invalid anchor")
            anchors[reference_gene].append(Anchor(chromosome, start, end))
    if not anchors:
        raise SpatialInputError(f"{path.name}: no anchors")
    return anchors


def anchor_scope(values: list[Anchor]) -> tuple[str, int, int] | None:
    chromosomes = {value.chromosome for value in values}
    if len(chromosomes) != 1:
        return None
    return next(iter(chromosomes)), min(value.start for value in values), max(value.end for value in values)


def nearest(genes: list[str], anchored: set[str], index: int, step: int, window: int) -> str | None:
    for distance in range(1, window + 1):
        candidate = index + step * distance
        if candidate < 0 or candidate >= len(genes):
            return None
        if genes[candidate] in anchored:
            return genes[candidate]
    return None


def write_tsv(path: Path, rows: list[dict[str, object]], columns: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--reference-coords", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--window", type=int, default=20)
    parser.add_argument("--padding-bp", type=int, default=10000)
    parser.add_argument("--maximum-locus-bp", type=int, default=5000000)
    parser.add_argument("--minimum-acgt-fraction", type=float, default=0.80)
    args = parser.parse_args()
    staging: Path | None = None
    try:
        if args.window < 1 or args.padding_bp < 0 or args.maximum_locus_bp < 1:
            raise SpatialInputError("invalid interval parameters")
        if not 0 <= args.minimum_acgt_fraction <= 1:
            raise SpatialInputError("minimum-acgt-fraction must lie in [0,1]")
        root = args.data_root.resolve()
        rows = read_manifest(args.manifest)
        reference_coords = resolve(root, args.reference_coords)
        lookup, order = read_reference_order(reference_coords)
        output = args.output_dir.resolve()
        if output.exists():
            raise SpatialInputError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        calls: list[dict[str, object]] = []
        coordinates: list[dict[str, object]] = []
        assemblies: list[dict[str, object]] = []
        audits: list[dict[str, object]] = []
        for row in rows:
            unit = row["assembly_unit_id"]
            sample = row["legacy_sample"]
            audit_path = resolve(root, row["synorth_audit"])
            pairs_path = resolve(root, row["synorth_pairs"])
            deleted_path = resolve(root, row["deleted_genes"], allow_empty=True)
            genome_path = resolve(root, row["genome"])
            gff_path = resolve(root, row["gff"])
            audit = read_json(audit_path)
            inputs = audit.get("inputs")
            metrics = audit.get("metrics")
            if audit.get("schema_version") != 2 or audit.get("sample") != sample or not isinstance(
                inputs, dict
            ) or not isinstance(metrics, dict):
                raise SpatialInputError(f"{unit}: incompatible SynOrths audit")
            if Path(str(inputs.get("pairs", ""))).resolve() != pairs_path.resolve():
                raise SpatialInputError(f"{unit}: pair/audit realpath mismatch")
            reference_column = int(inputs["reference_column_1_based"])
            query_column = int(inputs["query_column_1_based"])
            anchors = read_anchors(pairs_path, reference_column, query_column)
            if int(metrics.get("unique_reference_anchors", -1)) != len(anchors):
                raise SpatialInputError(f"{unit}: anchor count does not close to audit")
            genome = read_genome(genome_path)
            deleted = read_id_list(deleted_path)
            reasons: Counter[str] = Counter()
            localized = 0
            for gene in deleted:
                if gene not in lookup:
                    reasons["missing_reference_coordinate"] += 1
                    continue
                chromosome, index = lookup[gene]
                genes = order[chromosome]
                left = nearest(genes, set(anchors), index, -1, args.window)
                right = nearest(genes, set(anchors), index, 1, args.window)
                if left is None or right is None:
                    reasons["missing_bilateral_anchor"] += 1
                    continue
                left_scope = anchor_scope(anchors[left])
                right_scope = anchor_scope(anchors[right])
                if left_scope is None or right_scope is None:
                    reasons["ambiguous_anchor_target"] += 1
                    continue
                if left_scope[0] != right_scope[0]:
                    reasons["flanks_map_to_different_target_chromosomes"] += 1
                    continue
                target_chromosome = left_scope[0]
                if target_chromosome not in genome:
                    raise SpatialInputError(f"{unit}/{gene}: target chromosome absent from genome")
                start = max(1, min(left_scope[1], right_scope[1]) - args.padding_bp)
                end = min(len(genome[target_chromosome]), max(left_scope[2], right_scope[2]) + args.padding_bp)
                sequence = genome[target_chromosome][start - 1 : end]
                span = end - start + 1
                acgt = sum(base in "ACGT" for base in sequence) / len(sequence)
                if span > args.maximum_locus_bp:
                    reasons["excessive_target_interval"] += 1
                    continue
                if acgt < args.minimum_acgt_fraction:
                    reasons["low_unambiguous_sequence_fraction"] += 1
                    continue
                localized += 1
                reasons["localized"] += 1
                calls.append(
                    {
                        "assembly_unit_id": unit, "reference_gene_id": gene,
                        "classification": "positive_deleted", "callable": "true",
                    }
                )
                coordinates.append(
                    {
                        "assembly_unit_id": unit, "reference_gene_id": gene,
                        "classification": "positive_deleted", "chromosome": target_chromosome,
                        "expected_locus_start_1based": start, "expected_locus_end_1based": end,
                        "coordinate_semantics": "midpoint_of_bilateral_legacy_synorth_bounded_expected_interval",
                    }
                )
            assemblies.append(
                {
                    "assembly_unit_id": unit,
                    "biological_species": row["biological_species"],
                    "haplotype_or_subgenome": row["haplotype_or_subgenome"],
                    "assembly_scope": row["assembly_scope"],
                    "genome": os.path.relpath(genome_path.resolve(), staging),
                    "gff": os.path.relpath(gff_path.resolve(), staging),
                    "genome_local_sha256": sha256(genome_path),
                    "gff_local_sha256": sha256(gff_path),
                }
            )
            audits.append(
                {
                    "assembly_unit_id": unit, "legacy_sample": sample,
                    "historical_deleted_count": len(deleted), "spatially_localized_count": localized,
                    "localization_reasons": dict(sorted(reasons.items())),
                    "synorth_audit": binding(audit_path), "synorth_pairs": binding(pairs_path),
                    "deleted_genes": binding(deleted_path, allow_empty=True),
                    "genome": binding(genome_path), "gff": binding(gff_path),
                }
            )
        calls.sort(key=lambda row: (str(row["assembly_unit_id"]), str(row["reference_gene_id"])))
        coordinates.sort(key=lambda row: (str(row["assembly_unit_id"]), str(row["reference_gene_id"])))
        assemblies.sort(key=lambda row: str(row["assembly_unit_id"]))
        call_path = staging / "positive_deleted_calls.tsv"
        coordinate_path = staging / "expected_deleted_locus_coordinates.tsv"
        assembly_path = staging / "assembly_manifest.tsv"
        write_tsv(call_path, calls, ("assembly_unit_id", "reference_gene_id", "classification", "callable"))
        write_tsv(
            coordinate_path, coordinates,
            (
                "assembly_unit_id", "reference_gene_id", "classification", "chromosome",
                "expected_locus_start_1based", "expected_locus_end_1based", "coordinate_semantics",
            ),
        )
        write_tsv(
            assembly_path, assemblies,
            (
                "assembly_unit_id", "biological_species", "haplotype_or_subgenome", "assembly_scope",
                "genome", "gff", "genome_local_sha256", "gff_local_sha256",
            ),
        )
        report = {
            "schema_version": 1,
            "workflow": "legacy_conservative_deleted_expected_locus_spatial_inputs",
            "status": "PASS",
            "coordinate_semantics": "expected interval bounded by bilateral exact-bound legacy SynOrths anchors; midpoint is not an observed remnant",
            "chromosome_labels": "legacy analyzed assembly labels; primary spatial comparison uses normalized within-chromosome distance",
            "parameters": {
                "reference_gene_window": args.window,
                "target_interval_padding_bp": args.padding_bp,
                "maximum_target_interval_bp": args.maximum_locus_bp,
                "minimum_target_interval_acgt_fraction": args.minimum_acgt_fraction,
            },
            "unit_count": len(rows),
            "historical_deleted_count": sum(item["historical_deleted_count"] for item in audits),
            "spatially_localized_positive_deleted_count": len(calls),
            "source_manifest": binding(args.manifest),
            "reference_coords": binding(reference_coords),
            "units": audits,
            "outputs": {
                "positive_calls": binding(call_path, allow_empty=True),
                "expected_locus_coordinates": binding(coordinate_path, allow_empty=True),
                "assembly_manifest": binding(assembly_path),
            },
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        with (staging / "checksums.tsv").open("w", encoding="utf-8", newline="") as handle:
            handle.write("file\tbytes\tsha256\n")
            for path in sorted(staging.iterdir()):
                if path.is_file() and path.name != "checksums.tsv":
                    item = binding(path, allow_empty=True)
                    handle.write(f"{path.name}\t{item['bytes']}\t{item['sha256']}\n")
        os.rename(staging, output)
        staging = None
        print(json.dumps({"status": "PASS", "units": len(rows), "localized": len(calls)}))
        return 0
    except (SpatialInputError, OSError, UnicodeError, csv.Error, ValueError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    finally:
        if staging is not None and staging.exists():
            shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
