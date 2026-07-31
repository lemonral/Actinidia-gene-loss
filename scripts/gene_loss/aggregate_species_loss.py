#!/usr/bin/env python3
"""Aggregate complete assembly-unit calls into biological-species evidence.

This program deliberately separates technical assembly units (haplotypes or
subgenomes) from biological species.  A species aggregation rule must be
declared in metadata; no rule is inferred from the number or names of units.
The conservative ``positive_complete`` state means every included unit has a
positive call.  ``positive_partial`` records mixed unit evidence and is never,
by itself, labelled a lineage-restricted species loss.
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
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from geneloss_repro.io_utils import SchemaError  # noqa: E402
from geneloss_repro.pgls import (  # noqa: E402
    _capture_snapshot,
    _rename_directory_no_replace,
    _require_unchanged,
)


POSITIVE_CLASSES = frozenset({"pseudogenized", "deleted"})
NEGATIVE_CLASSES = frozenset({"retained"})
ALLOWED_CLASSES = POSITIVE_CLASSES | NEGATIVE_CLASSES | {"not_called_loss", "uncertain"}
ALLOWED_RULES = frozenset({"all_units_positive", "any_unit_positive"})

UNIT_COLUMNS = [
    "reference_gene_id",
    "assembly_unit_id",
    "biological_species",
    "aggregation_rule",
    "classification",
    "callable",
    "evidence_state",
    "positive_call",
    "confident_negative",
]
SPECIES_MATRIX_COLUMNS = [
    "reference_gene_id",
    "biological_species",
    "aggregation_rule",
    "species_gene_status",
    "species_positive_by_rule",
    "assembly_unit_count",
    "callable_unit_count",
    "positive_unit_count",
    "pseudogenized_unit_count",
    "deleted_unit_count",
    "confident_negative_unit_count",
    "uncertain_unit_count",
    "positive_units",
    "confident_negative_units",
    "uncertain_units",
]
PREVALENCE_COLUMNS = [
    "reference_gene_id",
    "biological_species_count",
    "positive_complete_species_count",
    "positive_partial_species_count",
    "positive_by_rule_species_count",
    "not_positive_species_count",
    "uncertain_species_count",
    "positive_complete_prevalence",
    "positive_by_rule_prevalence",
    "shared_positive_complete",
    "confident_lineage_restricted_species_loss",
    "positive_complete_species",
    "positive_partial_species",
    "positive_by_rule_species",
    "not_positive_species",
    "uncertain_species",
]
NON_SHARED_COLUMNS = SPECIES_MATRIX_COLUMNS + [
    "positive_complete_species_count",
    "positive_partial_species_count",
    "positive_by_rule_species_count",
    "confident_lineage_restricted_species_loss",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_boolean(value: str, *, context: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "t", "yes", "y"}:
        return True
    if normalized in {"0", "false", "f", "no", "n"}:
        return False
    raise ValueError(f"{context}: expected an explicit boolean, found {value!r}")


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            yield handle
    else:
        with path.open(encoding="utf-8", newline="") as handle:
            yield handle


def read_metadata(
    path: Path, include_column: str
) -> tuple[dict[str, dict[str, str]], dict[str, list[str]], set[str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"assembly_unit_id", "biological_species", "aggregation_rule", include_column}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path.name}: missing metadata columns: {', '.join(sorted(missing))}")
        rows = list(reader)

    if not rows:
        raise ValueError(f"{path.name}: metadata has no data rows")

    seen_units: set[str] = set()
    selected: dict[str, dict[str, str]] = {}
    units_by_species: dict[str, list[str]] = defaultdict(list)
    rules_by_species: dict[str, set[str]] = defaultdict(set)
    for line_number, row in enumerate(rows, start=2):
        unit = row["assembly_unit_id"].strip()
        species = row["biological_species"].strip()
        if not unit:
            raise ValueError(f"{path.name}:{line_number}: empty assembly_unit_id")
        if not species:
            raise ValueError(f"{path.name}:{line_number}: empty biological_species")
        if unit in seen_units:
            raise ValueError(f"{path.name}:{line_number}: duplicate assembly_unit_id {unit!r}")
        seen_units.add(unit)
        included = parse_boolean(
            row[include_column], context=f"{path.name}:{line_number}:{include_column}"
        )
        if not included:
            continue
        rule = row["aggregation_rule"].strip()
        if not rule:
            raise ValueError(
                f"{path.name}:{line_number}: included unit {unit!r} has no aggregation_rule"
            )
        if rule not in ALLOWED_RULES:
            raise ValueError(
                f"{path.name}:{line_number}: unsupported aggregation_rule {rule!r}; "
                f"choose one of {', '.join(sorted(ALLOWED_RULES))}"
            )
        selected[unit] = {
            "assembly_unit_id": unit,
            "biological_species": species,
            "aggregation_rule": rule,
        }
        units_by_species[species].append(unit)
        rules_by_species[species].add(rule)

    if not selected:
        raise ValueError(f"{path.name}: no assembly units selected by {include_column!r}")
    inconsistent = {species: rules for species, rules in rules_by_species.items() if len(rules) != 1}
    if inconsistent:
        species = sorted(inconsistent)[0]
        rules = ", ".join(sorted(inconsistent[species]))
        raise ValueError(
            f"{path.name}: biological species {species!r} has inconsistent aggregation rules: {rules}"
        )
    for species in units_by_species:
        units_by_species[species].sort()
    return selected, dict(units_by_species), seen_units


def evidence_state(classification: str, callable_value: bool) -> str:
    if classification in POSITIVE_CLASSES:
        return "positive"
    if classification in NEGATIVE_CLASSES and callable_value:
        return "confident_negative"
    return "uncertain"


def read_complete_matrix(
    path: Path,
    selected_metadata: dict[str, dict[str, str]],
    known_metadata_units: set[str],
) -> tuple[list[str], dict[str, dict[str, tuple[str, bool, str]]]]:
    """Read and validate the exact selected unit-by-gene grid."""
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"reference_gene_id", "assembly_unit_id", "classification", "callable"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path.name}: missing matrix columns: {', '.join(sorted(missing))}")

        calls_by_gene: dict[str, dict[str, tuple[str, bool, str]]] = defaultdict(dict)
        all_pairs: set[tuple[str, str]] = set()
        selected_units_seen: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            gene = row["reference_gene_id"].strip()
            unit = row["assembly_unit_id"].strip()
            classification = row["classification"].strip()
            if not gene:
                raise ValueError(f"{path.name}:{line_number}: empty reference_gene_id")
            if not unit:
                raise ValueError(f"{path.name}:{line_number}: empty assembly_unit_id")
            if unit not in known_metadata_units:
                raise ValueError(
                    f"{path.name}:{line_number}: assembly unit {unit!r} is absent from metadata"
                )
            if classification not in ALLOWED_CLASSES:
                raise ValueError(
                    f"{path.name}:{line_number}: unsupported classification {classification!r}"
                )
            callable_value = parse_boolean(
                row["callable"], context=f"{path.name}:{line_number}:callable"
            )
            if classification in POSITIVE_CLASSES and not callable_value:
                raise ValueError(
                    f"{path.name}:{line_number}: positive classification {classification!r} "
                    "requires callable=true"
                )
            pair = (gene, unit)
            if pair in all_pairs:
                raise ValueError(
                    f"{path.name}:{line_number}: duplicate unit-gene row for {gene!r}, {unit!r}"
                )
            all_pairs.add(pair)
            # Infer one reference universe from the complete input, including
            # rows belonging to metadata units excluded from this cohort.
            calls_by_gene.setdefault(gene, {})
            if unit not in selected_metadata:
                continue
            selected_units_seen.add(unit)
            calls_by_gene[gene][unit] = (
                classification,
                callable_value,
                evidence_state(classification, callable_value),
            )

    if not calls_by_gene:
        raise ValueError(f"{path.name}: selected call matrix has no data rows")
    missing_units = set(selected_metadata).difference(selected_units_seen)
    if missing_units:
        raise ValueError(
            f"{path.name}: selected assembly units have no matrix rows: {', '.join(sorted(missing_units))}"
        )

    selected_units = set(selected_metadata)
    missing_pairs: list[tuple[str, str]] = []
    for gene, unit_calls in calls_by_gene.items():
        for unit in selected_units.difference(unit_calls):
            missing_pairs.append((gene, unit))
    if missing_pairs:
        examples = ", ".join(f"{gene}/{unit}" for gene, unit in sorted(missing_pairs)[:5])
        raise ValueError(
            f"{path.name}: incomplete selected unit-by-gene grid; "
            f"{len(missing_pairs)} rows are missing ({examples})"
        )
    expected = len(calls_by_gene) * len(selected_units)
    observed = sum(len(unit_calls) for unit_calls in calls_by_gene.values())
    if observed != expected:
        raise ValueError(f"{path.name}: observed {observed} selected rows; expected {expected}")
    return sorted(calls_by_gene), dict(calls_by_gene)


def species_record(
    gene: str,
    species: str,
    units: list[str],
    rule: str,
    calls: dict[str, tuple[str, bool, str]],
) -> dict[str, object]:
    positive_units = [unit for unit in units if calls[unit][2] == "positive"]
    negative_units = [unit for unit in units if calls[unit][2] == "confident_negative"]
    uncertain_units = [unit for unit in units if calls[unit][2] == "uncertain"]
    positive_count = len(positive_units)
    if positive_count == len(units):
        status = "positive_complete"
    elif positive_count:
        status = "positive_partial"
    elif len(negative_units) == len(units):
        status = "not_positive"
    else:
        status = "uncertain"
    species_positive = status == "positive_complete" or (
        rule == "any_unit_positive" and status == "positive_partial"
    )
    return {
        "reference_gene_id": gene,
        "biological_species": species,
        "aggregation_rule": rule,
        "species_gene_status": status,
        "species_positive_by_rule": str(species_positive).lower(),
        "assembly_unit_count": len(units),
        "callable_unit_count": sum(1 for unit in units if calls[unit][1]),
        "positive_unit_count": positive_count,
        "pseudogenized_unit_count": sum(1 for unit in units if calls[unit][0] == "pseudogenized"),
        "deleted_unit_count": sum(1 for unit in units if calls[unit][0] == "deleted"),
        "confident_negative_unit_count": len(negative_units),
        "uncertain_unit_count": len(uncertain_units),
        "positive_units": ";".join(positive_units),
        "confident_negative_units": ";".join(negative_units),
        "uncertain_units": ";".join(uncertain_units),
    }


def write_header(handle: TextIO, columns: list[str]) -> csv.DictWriter:
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    return writer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-call-matrix", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--include-column", default="include")
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if os.path.lexists(args.output_dir):
        raise SystemExit(f"ERROR: output path already exists: {args.output_dir}")

    try:
        input_snapshots = {
            "unit_call_matrix": _capture_snapshot(args.unit_call_matrix),
            "unit_metadata": _capture_snapshot(args.unit_metadata),
        }
        with tempfile.TemporaryDirectory(prefix="species-loss-input-snapshots.") as frozen:
            frozen_root = Path(frozen)
            frozen_matrix = frozen_root / f"unit_call_matrix.{args.unit_call_matrix.name}"
            frozen_metadata = frozen_root / f"unit_metadata.{args.unit_metadata.name}"
            frozen_matrix.write_bytes(input_snapshots["unit_call_matrix"].payload)
            frozen_metadata.write_bytes(input_snapshots["unit_metadata"].payload)
            selected, units_by_species, known_units = read_metadata(
                frozen_metadata, args.include_column
            )
            genes, calls_by_gene = read_complete_matrix(
                frozen_matrix, selected, known_units
            )
    except (OSError, ValueError, SchemaError, csv.Error) as exc:
        raise SystemExit(f"ERROR: {exc}")

    species_names = sorted(units_by_species)
    rule_by_species = {
        species: selected[units_by_species[species][0]]["aggregation_rule"]
        for species in species_names
    }
    unit_names = sorted(selected)

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = args.output_dir.parent / f".{args.output_dir.name}.tmp.{os.getpid()}"
    if staging.exists():
        raise SystemExit(f"ERROR: staging directory already exists: {staging}")
    staging.mkdir()

    status_counts: Counter[str] = Counter()
    shared_count = 0
    non_shared_positive_call_count = 0
    confident_lineage_count = 0
    try:
        unit_path = staging / "unit_calls_long.tsv"
        with unit_path.open("w", encoding="utf-8", newline="") as handle:
            writer = write_header(handle, UNIT_COLUMNS)
            for gene in genes:
                for unit in unit_names:
                    classification, callable_value, state = calls_by_gene[gene][unit]
                    writer.writerow(
                        {
                            "reference_gene_id": gene,
                            "assembly_unit_id": unit,
                            "biological_species": selected[unit]["biological_species"],
                            "aggregation_rule": selected[unit]["aggregation_rule"],
                            "classification": classification,
                            "callable": str(callable_value).lower(),
                            "evidence_state": state,
                            "positive_call": str(state == "positive").lower(),
                            "confident_negative": str(state == "confident_negative").lower(),
                        }
                    )

        matrix_path = staging / "species_gene_matrix.tsv"
        prevalence_path = staging / "species_prevalence.tsv"
        shared_path = staging / "shared_positive_complete_genes.tsv"
        non_shared_path = staging / "non_shared_positive_calls.tsv"
        with (
            matrix_path.open("w", encoding="utf-8", newline="") as matrix_handle,
            prevalence_path.open("w", encoding="utf-8", newline="") as prevalence_handle,
            shared_path.open("w", encoding="utf-8", newline="") as shared_handle,
            non_shared_path.open("w", encoding="utf-8", newline="") as non_shared_handle,
        ):
            matrix_writer = write_header(matrix_handle, SPECIES_MATRIX_COLUMNS)
            prevalence_writer = write_header(prevalence_handle, PREVALENCE_COLUMNS)
            shared_writer = write_header(shared_handle, PREVALENCE_COLUMNS)
            non_shared_writer = write_header(non_shared_handle, NON_SHARED_COLUMNS)

            for gene in genes:
                species_rows = [
                    species_record(
                        gene,
                        species,
                        units_by_species[species],
                        rule_by_species[species],
                        calls_by_gene[gene],
                    )
                    for species in species_names
                ]
                for row in species_rows:
                    matrix_writer.writerow(row)
                    status_counts[str(row["species_gene_status"])] += 1

                by_status = {
                    status: [
                        str(row["biological_species"])
                        for row in species_rows
                        if row["species_gene_status"] == status
                    ]
                    for status in ("positive_complete", "positive_partial", "not_positive", "uncertain")
                }
                by_rule = [
                    str(row["biological_species"])
                    for row in species_rows
                    if row["species_positive_by_rule"] == "true"
                ]
                shared = len(by_status["positive_complete"]) == len(species_names)
                confident_lineage = (
                    not shared
                    and len(by_status["positive_complete"]) == 1
                    and not by_status["positive_partial"]
                    and not by_status["uncertain"]
                    and len(by_status["not_positive"]) == len(species_names) - 1
                )
                prevalence = {
                    "reference_gene_id": gene,
                    "biological_species_count": len(species_names),
                    "positive_complete_species_count": len(by_status["positive_complete"]),
                    "positive_partial_species_count": len(by_status["positive_partial"]),
                    "positive_by_rule_species_count": len(by_rule),
                    "not_positive_species_count": len(by_status["not_positive"]),
                    "uncertain_species_count": len(by_status["uncertain"]),
                    "positive_complete_prevalence": format(
                        len(by_status["positive_complete"]) / len(species_names), ".12g"
                    ),
                    "positive_by_rule_prevalence": format(len(by_rule) / len(species_names), ".12g"),
                    "shared_positive_complete": str(shared).lower(),
                    "confident_lineage_restricted_species_loss": str(confident_lineage).lower(),
                    "positive_complete_species": ";".join(by_status["positive_complete"]),
                    "positive_partial_species": ";".join(by_status["positive_partial"]),
                    "positive_by_rule_species": ";".join(by_rule),
                    "not_positive_species": ";".join(by_status["not_positive"]),
                    "uncertain_species": ";".join(by_status["uncertain"]),
                }
                prevalence_writer.writerow(prevalence)
                if shared:
                    shared_writer.writerow(prevalence)
                    shared_count += 1
                if confident_lineage:
                    confident_lineage_count += 1
                if not shared:
                    for row in species_rows:
                        if row["species_positive_by_rule"] != "true":
                            continue
                        non_shared_writer.writerow(
                            {
                                **row,
                                "positive_complete_species_count": len(by_status["positive_complete"]),
                                "positive_partial_species_count": len(by_status["positive_partial"]),
                                "positive_by_rule_species_count": len(by_rule),
                                "confident_lineage_restricted_species_loss": str(
                                    confident_lineage
                                    and row["species_gene_status"] == "positive_complete"
                                ).lower(),
                            }
                        )
                        non_shared_positive_call_count += 1

        output_files = [unit_path, matrix_path, prevalence_path, shared_path, non_shared_path]
        summary = {
            "schema_version": "2.0",
            "status": "PASS",
            "definitions": {
                "shared_positive_complete": (
                    "positive_complete in every included biological species; technical assembly "
                    "units are never counted as species"
                ),
                "positive_complete": "every included assembly unit for the species has a positive call",
                "positive_partial": "at least one but not every included assembly unit has a positive call",
                "not_positive": "all included units are callable and confidently negative",
                "uncertain": "no positive unit and at least one uncertain or non-callable unit",
                "not_called_loss": (
                    "absence from a positive list is not retained evidence and is treated as uncertain"
                ),
                "lineage_restricted": (
                    "one positive_complete species, every other species not_positive, and no partial "
                    "or uncertain species"
                ),
            },
            "inputs": [
                {
                    "role": "unit_call_matrix",
                    "basename": args.unit_call_matrix.name,
                    "sha256": input_snapshots["unit_call_matrix"].sha256,
                },
                {
                    "role": "unit_metadata",
                    "basename": args.unit_metadata.name,
                    "sha256": input_snapshots["unit_metadata"].sha256,
                },
            ],
            "include_column": args.include_column,
            "assembly_unit_count": len(unit_names),
            "biological_species_count": len(species_names),
            "reference_gene_count": len(genes),
            "expected_unit_matrix_rows": len(unit_names) * len(genes),
            "species_gene_matrix_rows": len(species_names) * len(genes),
            "species_gene_status_counts": dict(sorted(status_counts.items())),
            "shared_positive_complete_gene_count": shared_count,
            "non_shared_positive_call_count": non_shared_positive_call_count,
            "confident_lineage_restricted_gene_count": confident_lineage_count,
            "aggregation_rule_species_counts": dict(sorted(Counter(rule_by_species.values()).items())),
            "species_aggregation": [
                {
                    "biological_species": species,
                    "aggregation_rule": rule_by_species[species],
                    "assembly_unit_count": len(units_by_species[species]),
                    "assembly_units": units_by_species[species],
                }
                for species in species_names
            ],
            "checks": {
                "complete_selected_unit_gene_grid": True,
                "positive_classification_requires_callable": True,
                "not_called_loss_treated_as_uncertain": True,
                "species_status_counts_reconciled": (
                    sum(status_counts.values()) == len(species_names) * len(genes)
                ),
                "shared_set_reconciled": shared_count
                == sum(
                    1
                    for gene in genes
                    if all(
                        species_record(
                            gene,
                            species,
                            units_by_species[species],
                            rule_by_species[species],
                            calls_by_gene[gene],
                        )["species_gene_status"]
                        == "positive_complete"
                        for species in species_names
                    )
                ),
                "output_checksums_reconciled": True,
            },
            "outputs": [
                {"basename": path.name, "sha256": sha256(path)} for path in output_files
            ],
        }
        (staging / "species_loss_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

        for snapshot in input_snapshots.values():
            _require_unchanged(snapshot)
        _rename_directory_no_replace(staging, args.output_dir)
    except (OSError, ValueError, SchemaError, csv.Error) as exc:
        shutil.rmtree(staging, ignore_errors=True)
        raise SystemExit(f"ERROR: {exc}")
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    print(
        f"Aggregated {len(genes)} reference genes across {len(unit_names)} assembly units "
        f"and {len(species_names)} biological species into {args.output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
