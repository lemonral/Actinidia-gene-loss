#!/usr/bin/env python3
"""Summarize complete/partial loss and callable copy-opportunity rates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


class SummaryError(RuntimeError):
    pass


UNIT_COLUMNS = (
    "reference_gene_id", "assembly_unit_id", "biological_species", "aggregation_rule",
    "classification", "callable", "evidence_state", "positive_call", "confident_negative",
)
SPECIES_REQUIRED = {
    "reference_gene_id", "biological_species", "species_gene_status", "assembly_unit_count",
}
ALLOWED_SPECIES_STATES = {"positive_complete", "positive_partial", "not_positive", "uncertain"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    source = path.resolve()
    if not source.is_file() or (not allow_empty and source.stat().st_size == 0):
        raise SummaryError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def boolean(value: str, context: str) -> bool:
    normalized = value.strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise SummaryError(f"{context}: expected true/false, found {value!r}")


def rate(numerator: int, denominator: int) -> str:
    return "NA" if denominator == 0 else format(numerator / denominator, ".12g")


def read_shared_ids(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "reference_gene_id" not in (reader.fieldnames or []):
            raise SummaryError(f"{path.name}: missing reference_gene_id")
        values = [row["reference_gene_id"] for row in reader]
    if any(not value for value in values) or len(values) != len(set(values)):
        raise SummaryError(f"{path.name}: empty/duplicate shared reference IDs")
    return set(values)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-calls", required=True, type=Path)
    parser.add_argument("--species-matrix", required=True, type=Path)
    parser.add_argument("--shared-genes", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    staging: Path | None = None
    try:
        for path in (args.unit_calls, args.species_matrix, args.shared_genes):
            binding(path, allow_empty=False)
        shared = read_shared_ids(args.shared_genes)
        unit_counts: dict[str, Counter[str]] = defaultdict(Counter)
        species_copy_counts: dict[str, Counter[str]] = defaultdict(Counter)
        unit_species: dict[str, str] = {}
        reference_ids: set[str] = set()
        unit_pairs: set[tuple[str, str]] = set()
        with args.unit_calls.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != UNIT_COLUMNS:
                raise SummaryError("unit-call columns differ from exact aggregation schema")
            for line_number, row in enumerate(reader, 2):
                gene = row["reference_gene_id"]
                unit = row["assembly_unit_id"]
                species = row["biological_species"]
                if not gene or not unit or not species:
                    raise SummaryError(f"unit calls:{line_number}: empty gene/unit/species")
                pair = (gene, unit)
                if pair in unit_pairs:
                    raise SummaryError(f"unit calls:{line_number}: duplicate unit/gene pair")
                unit_pairs.add(pair)
                if unit in unit_species and unit_species[unit] != species:
                    raise SummaryError(f"unit calls:{line_number}: unit maps to multiple species")
                unit_species[unit] = species
                callable_value = boolean(row["callable"], f"unit calls:{line_number}:callable")
                positive = boolean(row["positive_call"], f"unit calls:{line_number}:positive_call")
                if positive and not callable_value:
                    raise SummaryError(f"unit calls:{line_number}: positive opportunity is non-callable")
                reference_ids.add(gene)
                for counts in (unit_counts[unit], species_copy_counts[species]):
                    counts["total"] += 1
                    counts["callable"] += int(callable_value)
                    counts["positive"] += int(positive)
                    if gene not in shared:
                        counts["nonshared_total"] += 1
                        counts["nonshared_callable"] += int(callable_value)
                        counts["nonshared_positive"] += int(positive)

        units = sorted(unit_counts)
        if not units or not reference_ids or len(unit_pairs) != len(units) * len(reference_ids):
            raise SummaryError("unit call grid is incomplete")
        species_state_counts: dict[str, Counter[str]] = defaultdict(Counter)
        species_unit_counts: dict[str, int] = {}
        species_pairs: set[tuple[str, str]] = set()
        with args.species_matrix.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not SPECIES_REQUIRED.issubset(reader.fieldnames or []):
                raise SummaryError("species matrix is missing required columns")
            for line_number, row in enumerate(reader, 2):
                gene = row["reference_gene_id"]
                species = row["biological_species"]
                status = row["species_gene_status"]
                if gene not in reference_ids or not species or status not in ALLOWED_SPECIES_STATES:
                    raise SummaryError(f"species matrix:{line_number}: invalid gene/species/status")
                pair = (gene, species)
                if pair in species_pairs:
                    raise SummaryError(f"species matrix:{line_number}: duplicate species/gene pair")
                species_pairs.add(pair)
                try:
                    count = int(row["assembly_unit_count"])
                except ValueError as error:
                    raise SummaryError(f"species matrix:{line_number}: invalid unit count") from error
                if species in species_unit_counts and species_unit_counts[species] != count:
                    raise SummaryError(f"species matrix:{line_number}: inconsistent unit count")
                species_unit_counts[species] = count
                species_state_counts[species][status] += 1
                if gene not in shared:
                    species_state_counts[species][f"nonshared_{status}"] += 1
        species_names = sorted(species_state_counts)
        if (
            set(species_names) != set(species_copy_counts)
            or len(species_pairs) != len(species_names) * len(reference_ids)
        ):
            raise SummaryError("species matrix/grid does not close to unit calls")
        if not shared.issubset(reference_ids):
            raise SummaryError("shared genes lie outside reference universe")
        for species in species_names:
            observed_units = sum(1 for value in unit_species.values() if value == species)
            if species_unit_counts[species] != observed_units:
                raise SummaryError(f"{species}: species/unit count closure failed")

        output = args.output_dir.absolute()
        if output.exists():
            raise SummaryError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        unit_path = staging / "unit_copy_opportunity_summary.tsv"
        unit_fields = (
            "assembly_unit_id", "biological_species", "reference_gene_count",
            "callable_copy_opportunities", "positive_loss_copy_opportunities",
            "callable_copy_opportunity_loss_rate", "nonshared_reference_gene_count",
            "nonshared_callable_copy_opportunities", "nonshared_positive_loss_copy_opportunities",
            "nonshared_callable_copy_opportunity_loss_rate",
        )
        with unit_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=unit_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for unit in units:
                counts = unit_counts[unit]
                writer.writerow(
                    {
                        "assembly_unit_id": unit,
                        "biological_species": unit_species[unit],
                        "reference_gene_count": len(reference_ids),
                        "callable_copy_opportunities": counts["callable"],
                        "positive_loss_copy_opportunities": counts["positive"],
                        "callable_copy_opportunity_loss_rate": rate(counts["positive"], counts["callable"]),
                        "nonshared_reference_gene_count": len(reference_ids) - len(shared),
                        "nonshared_callable_copy_opportunities": counts["nonshared_callable"],
                        "nonshared_positive_loss_copy_opportunities": counts["nonshared_positive"],
                        "nonshared_callable_copy_opportunity_loss_rate": rate(
                            counts["nonshared_positive"], counts["nonshared_callable"]
                        ),
                    }
                )

        species_path = staging / "species_loss_mode_and_copy_opportunity_summary.tsv"
        species_fields = (
            "biological_species", "assembly_unit_count", "reference_gene_count",
            "total_copy_opportunities", "callable_copy_opportunities",
            "positive_loss_copy_opportunities", "callable_copy_opportunity_loss_rate",
            "complete_loss_gene_count", "partial_homeolog_loss_gene_count",
            "not_positive_gene_count", "uncertain_gene_count", "shared_complete_loss_gene_count",
            "nonshared_reference_gene_count", "nonshared_callable_copy_opportunities",
            "nonshared_positive_loss_copy_opportunities", "nonshared_callable_copy_opportunity_loss_rate",
            "nonshared_complete_loss_gene_count", "nonshared_partial_homeolog_loss_gene_count",
        )
        with species_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=species_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for species in species_names:
                copy = species_copy_counts[species]
                states = species_state_counts[species]
                writer.writerow(
                    {
                        "biological_species": species,
                        "assembly_unit_count": species_unit_counts[species],
                        "reference_gene_count": len(reference_ids),
                        "total_copy_opportunities": copy["total"],
                        "callable_copy_opportunities": copy["callable"],
                        "positive_loss_copy_opportunities": copy["positive"],
                        "callable_copy_opportunity_loss_rate": rate(copy["positive"], copy["callable"]),
                        "complete_loss_gene_count": states["positive_complete"],
                        "partial_homeolog_loss_gene_count": states["positive_partial"],
                        "not_positive_gene_count": states["not_positive"],
                        "uncertain_gene_count": states["uncertain"],
                        "shared_complete_loss_gene_count": len(shared),
                        "nonshared_reference_gene_count": len(reference_ids) - len(shared),
                        "nonshared_callable_copy_opportunities": copy["nonshared_callable"],
                        "nonshared_positive_loss_copy_opportunities": copy["nonshared_positive"],
                        "nonshared_callable_copy_opportunity_loss_rate": rate(
                            copy["nonshared_positive"], copy["nonshared_callable"]
                        ),
                        "nonshared_complete_loss_gene_count": states["nonshared_positive_complete"],
                        "nonshared_partial_homeolog_loss_gene_count": states["nonshared_positive_partial"],
                    }
                )
        report = {
            "schema_version": 1,
            "workflow": "callable_copy_opportunity_and_loss_mode_summary",
            "status": "PASS",
            "definitions": {
                "complete_loss": "positive loss in every included callable assembly unit for the lineage",
                "partial_homeolog_loss": "positive loss in at least one but not every included assembly unit",
                "callable_copy_opportunity_loss_rate": "positive callable unit-gene opportunities divided by all callable unit-gene opportunities",
                "nonshared": "reference genes shared-positive-complete across every lineage are excluded from numerator and denominator",
            },
            "assembly_unit_count": len(units),
            "biological_species_count": len(species_names),
            "reference_gene_count": len(reference_ids),
            "shared_positive_complete_gene_count": len(shared),
            "inputs": {
                "unit_calls": binding(args.unit_calls),
                "species_matrix": binding(args.species_matrix),
                "shared_genes": binding(args.shared_genes),
            },
            "outputs": {
                "unit_summary": binding(unit_path),
                "species_summary": binding(species_path),
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
        print(json.dumps({"status": "PASS", "units": len(units), "lineages": len(species_names)}))
        return 0
    except (SummaryError, OSError, UnicodeError, csv.Error, ValueError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    finally:
        if staging is not None and staging.exists():
            shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
