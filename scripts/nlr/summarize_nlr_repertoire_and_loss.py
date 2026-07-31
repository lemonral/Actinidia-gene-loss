#!/usr/bin/env python3
"""Summarize complete NLR repertoires and positive reference-NLR loss calls.

The analysis unit is an explicitly named assembly unit, not a biological
species and not an anonymous code.  One invocation publishes exactly one
declared cohort.  Primary and A. rufa sensitivity cohorts must therefore be
run into separate output directories and can never share a denominator.

Required input schemas (tab-delimited)
--------------------------------------
metadata
    assembly_unit_id, biological_species, haplotype_or_subgenome,
    assembly_scope, include, analysis_cohort
repertoire counts
    assembly_unit_id, assembly_scope, total_nlr_count,
    repertoire_source_basename, repertoire_source_sha256
positive calls
    assembly_unit_id, assembly_scope, reference_nlr_id,
    reference_nlr_universe_id
callable denominators
    Either one count row per unit with callable_reference_nlr_denominator, or
    one catalog row per callable reference with reference_nlr_id.  Both forms
    also require assembly_unit_id, assembly_scope, reference_nlr_universe_id,
    denominator_source_basename, and denominator_source_sha256.

An empty positive-call table is valid when its header is present.  A zero
callable denominator is valid only with zero positive calls; its percentage is
reported as an empty field with an explicit undefined_zero_denominator status.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


COHORT_ROLES = frozenset({"primary", "a_rufa_sensitivity"})
SHA256_PATTERN = re.compile(r"^[0-9a-fA-F]{64}$")
UNIT_COLUMNS = (
    "analysis_cohort",
    "cohort_role",
    "assembly_unit_id",
    "biological_species",
    "haplotype_or_subgenome",
    "assembly_scope",
    "total_nlr_count",
    "positive_reference_nlr_loss_count",
    "callable_reference_nlr_denominator",
    "positive_reference_nlr_loss_percentage",
    "percentage_status",
    "reference_nlr_universe_id",
    "denominator_source_basename",
    "denominator_source_sha256",
    "repertoire_source_basename",
    "repertoire_source_sha256",
    "positive_loss_calls_source_basename",
    "positive_loss_calls_source_sha256",
)
SPECIES_COLUMNS = (
    "analysis_cohort",
    "cohort_role",
    "biological_species",
    "assembly_unit_count",
    "assembly_unit_ids",
    "total_nlr_count_sum_across_units",
    "positive_reference_nlr_loss_count_sum_across_units",
    "callable_reference_nlr_denominator_sum_across_units",
    "positive_reference_nlr_loss_percentage_across_unit_comparisons",
    "percentage_status",
)


@dataclass(frozen=True)
class Denominator:
    count: int
    assembly_scope: str
    universe_id: str
    source_basename: str
    source_sha256: str
    callable_ids: frozenset[str] | None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--repertoire-counts", type=Path, required=True)
    parser.add_argument("--positive-loss-calls", type=Path, required=True)
    parser.add_argument("--callable-denominators", type=Path, required=True)
    parser.add_argument("--analysis-cohort", required=True)
    parser.add_argument("--cohort-role", required=True, choices=sorted(COHORT_ROLES))
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, *, required: Iterable[str], label: str) -> tuple[list[str], list[dict[str, str]]]:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{label} must be a regular non-symlink file: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not fields:
            raise ValueError(f"{label} has no header")
        if len(set(fields)) != len(fields):
            raise ValueError(f"{label} has duplicate column names")
        missing = set(required).difference(fields)
        if missing:
            raise ValueError(f"{label} is missing columns: {', '.join(sorted(missing))}")
        rows = []
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise ValueError(f"{label}:{line_number} has extra tab-delimited fields")
            rows.append({field: (raw.get(field) or "").strip() for field in fields})
    return fields, rows


def parse_boolean(value: str, *, context: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "t", "yes", "y", "1"}:
        return True
    if normalized in {"false", "f", "no", "n", "0"}:
        return False
    raise ValueError(f"{context}: expected an explicit boolean, found {value!r}")


def parse_nonnegative_integer(value: str, *, context: str) -> int:
    if not re.fullmatch(r"0|[1-9][0-9]*", value):
        raise ValueError(f"{context}: expected a non-negative integer, found {value!r}")
    return int(value)


def require_text(value: str, *, context: str) -> str:
    if not value:
        raise ValueError(f"{context}: value must not be empty")
    if any(character in value for character in "\t\r\n"):
        raise ValueError(f"{context}: value contains a tab or newline")
    return value


def require_basename(value: str, *, context: str) -> str:
    require_text(value, context=context)
    if Path(value).name != value or "/" in value or "\\" in value or value in {".", ".."}:
        raise ValueError(f"{context}: expected a basename without a directory path")
    return value


def require_sha256(value: str, *, context: str) -> str:
    if SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{context}: expected a 64-character SHA-256 digest")
    return value.lower()


def read_selected_metadata(path: Path, cohort: str) -> tuple[list[str], dict[str, dict[str, str]]]:
    required = {
        "assembly_unit_id",
        "biological_species",
        "haplotype_or_subgenome",
        "assembly_scope",
        "include",
        "analysis_cohort",
    }
    _, rows = read_tsv(path, required=required, label="metadata")
    if not rows:
        raise ValueError("metadata has no data rows")
    seen: set[str] = set()
    order: list[str] = []
    selected: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(rows, start=2):
        unit = require_text(row["assembly_unit_id"], context=f"metadata:{line_number}:assembly_unit_id")
        if unit in seen:
            raise ValueError(f"metadata:{line_number}: duplicate assembly_unit_id {unit!r}")
        seen.add(unit)
        included = parse_boolean(row["include"], context=f"metadata:{line_number}:include")
        if not included or row["analysis_cohort"] != cohort:
            continue
        require_text(row["biological_species"], context=f"metadata:{line_number}:biological_species")
        require_text(row["assembly_scope"], context=f"metadata:{line_number}:assembly_scope")
        selected[unit] = row
        order.append(unit)
    if not selected:
        raise ValueError(f"metadata selects no included assembly units for cohort {cohort!r}")
    return order, selected


def require_exact_units(observed: set[str], expected: set[str], *, label: str) -> None:
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        raise ValueError(
            f"{label} does not exactly match the included cohort; missing={missing}, extra={extra}"
        )


def read_repertoire_counts(path: Path, metadata: Mapping[str, Mapping[str, str]]) -> dict[str, dict[str, object]]:
    required = {
        "assembly_unit_id",
        "assembly_scope",
        "total_nlr_count",
        "repertoire_source_basename",
        "repertoire_source_sha256",
    }
    _, rows = read_tsv(path, required=required, label="repertoire counts")
    parsed: dict[str, dict[str, object]] = {}
    for line_number, row in enumerate(rows, start=2):
        unit = require_text(row["assembly_unit_id"], context=f"repertoire counts:{line_number}:assembly_unit_id")
        if unit in parsed:
            raise ValueError(f"repertoire counts:{line_number}: duplicate assembly_unit_id {unit!r}")
        parsed[unit] = {
            "assembly_scope": require_text(row["assembly_scope"], context=f"repertoire counts:{line_number}:assembly_scope"),
            "total_nlr_count": parse_nonnegative_integer(row["total_nlr_count"], context=f"repertoire counts:{line_number}:total_nlr_count"),
            "source_basename": require_basename(row["repertoire_source_basename"], context=f"repertoire counts:{line_number}:repertoire_source_basename"),
            "source_sha256": require_sha256(row["repertoire_source_sha256"], context=f"repertoire counts:{line_number}:repertoire_source_sha256"),
        }
    require_exact_units(set(parsed), set(metadata), label="repertoire counts")
    for unit, values in parsed.items():
        expected_scope = metadata[unit]["assembly_scope"]
        if values["assembly_scope"] != expected_scope:
            raise ValueError(
                f"repertoire counts: assembly scope mismatch for {unit!r}: "
                f"{values['assembly_scope']!r} != {expected_scope!r}"
            )
    return parsed


def _consistent_denominator_metadata(
    rows: Sequence[dict[str, str]], *, unit: str, label: str
) -> tuple[str, str, str, str]:
    fields = (
        "assembly_scope",
        "reference_nlr_universe_id",
        "denominator_source_basename",
        "denominator_source_sha256",
    )
    distinct = {tuple(row[field] for field in fields) for row in rows}
    if len(distinct) != 1:
        raise ValueError(f"{label}: inconsistent scope or denominator provenance for {unit!r}")
    scope, universe, basename, digest = next(iter(distinct))
    return (
        require_text(scope, context=f"{label}:{unit}:assembly_scope"),
        require_text(universe, context=f"{label}:{unit}:reference_nlr_universe_id"),
        require_basename(basename, context=f"{label}:{unit}:denominator_source_basename"),
        require_sha256(digest, context=f"{label}:{unit}:denominator_source_sha256"),
    )


def read_denominators(path: Path, metadata: Mapping[str, Mapping[str, str]]) -> tuple[dict[str, Denominator], str]:
    base_required = {
        "assembly_unit_id",
        "assembly_scope",
        "reference_nlr_universe_id",
        "denominator_source_basename",
        "denominator_source_sha256",
    }
    fields, rows = read_tsv(path, required=base_required, label="callable denominators")
    count_mode = "callable_reference_nlr_denominator" in fields
    catalog_mode = "reference_nlr_id" in fields
    if count_mode == catalog_mode:
        raise ValueError(
            "callable denominators must contain exactly one of "
            "callable_reference_nlr_denominator (count mode) or reference_nlr_id (catalog mode)"
        )
    parsed: dict[str, Denominator] = {}
    if count_mode:
        for line_number, row in enumerate(rows, start=2):
            unit = require_text(row["assembly_unit_id"], context=f"callable denominators:{line_number}:assembly_unit_id")
            if unit in parsed:
                raise ValueError(f"callable denominators:{line_number}: duplicate assembly_unit_id {unit!r}")
            scope, universe, basename, digest = _consistent_denominator_metadata(
                [row], unit=unit, label="callable denominators"
            )
            parsed[unit] = Denominator(
                count=parse_nonnegative_integer(
                    row["callable_reference_nlr_denominator"],
                    context=f"callable denominators:{line_number}:callable_reference_nlr_denominator",
                ),
                assembly_scope=scope,
                universe_id=universe,
                source_basename=basename,
                source_sha256=digest,
                callable_ids=None,
            )
        mode = "count"
    else:
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        seen_pairs: set[tuple[str, str]] = set()
        for line_number, row in enumerate(rows, start=2):
            unit = require_text(row["assembly_unit_id"], context=f"callable denominators:{line_number}:assembly_unit_id")
            reference = require_text(row["reference_nlr_id"], context=f"callable denominators:{line_number}:reference_nlr_id")
            pair = (unit, reference)
            if pair in seen_pairs:
                raise ValueError(f"callable denominators:{line_number}: duplicate callable pair {unit!r}/{reference!r}")
            seen_pairs.add(pair)
            grouped[unit].append(row)
        for unit, unit_rows in grouped.items():
            scope, universe, basename, digest = _consistent_denominator_metadata(
                unit_rows, unit=unit, label="callable denominators"
            )
            parsed[unit] = Denominator(
                count=len(unit_rows),
                assembly_scope=scope,
                universe_id=universe,
                source_basename=basename,
                source_sha256=digest,
                callable_ids=frozenset(row["reference_nlr_id"] for row in unit_rows),
            )
        mode = "catalog"
    require_exact_units(set(parsed), set(metadata), label="callable denominators")
    for unit, denominator in parsed.items():
        expected_scope = metadata[unit]["assembly_scope"]
        if denominator.assembly_scope != expected_scope:
            raise ValueError(
                f"callable denominators: assembly scope mismatch for {unit!r}: "
                f"{denominator.assembly_scope!r} != {expected_scope!r}"
            )
    return parsed, mode


def read_positive_calls(
    path: Path,
    metadata: Mapping[str, Mapping[str, str]],
    denominators: Mapping[str, Denominator],
) -> dict[str, int]:
    required = {
        "assembly_unit_id",
        "assembly_scope",
        "reference_nlr_id",
        "reference_nlr_universe_id",
    }
    _, rows = read_tsv(path, required=required, label="positive loss calls")
    counts = {unit: 0 for unit in metadata}
    seen_pairs: set[tuple[str, str]] = set()
    for line_number, row in enumerate(rows, start=2):
        unit = require_text(row["assembly_unit_id"], context=f"positive loss calls:{line_number}:assembly_unit_id")
        reference = require_text(row["reference_nlr_id"], context=f"positive loss calls:{line_number}:reference_nlr_id")
        if unit not in metadata:
            raise ValueError(
                f"positive loss calls:{line_number}: assembly unit {unit!r} is outside the exact included cohort"
            )
        pair = (unit, reference)
        if pair in seen_pairs:
            raise ValueError(f"positive loss calls:{line_number}: duplicate positive call {unit!r}/{reference!r}")
        seen_pairs.add(pair)
        expected_scope = metadata[unit]["assembly_scope"]
        if row["assembly_scope"] != expected_scope:
            raise ValueError(
                f"positive loss calls:{line_number}: assembly scope mismatch for {unit!r}: "
                f"{row['assembly_scope']!r} != {expected_scope!r}"
            )
        denominator = denominators[unit]
        if row["reference_nlr_universe_id"] != denominator.universe_id:
            raise ValueError(
                f"positive loss calls:{line_number}: reference universe mismatch for {unit!r}"
            )
        if denominator.callable_ids is not None and reference not in denominator.callable_ids:
            raise ValueError(
                f"positive loss calls:{line_number}: {reference!r} is absent from the callable denominator for {unit!r}"
            )
        counts[unit] += 1
    for unit, count in counts.items():
        denominator = denominators[unit].count
        if count > denominator:
            raise ValueError(
                f"positive reference-NLR loss count exceeds callable denominator for {unit!r}: "
                f"{count} > {denominator}"
            )
    return counts


def percentage(numerator: int, denominator: int) -> tuple[float | str, str]:
    if denominator == 0:
        return "", "undefined_zero_denominator"
    return round(100.0 * numerator / denominator, 6), "defined"


def build_unit_rows(
    *,
    order: Sequence[str],
    metadata: Mapping[str, Mapping[str, str]],
    repertoire: Mapping[str, Mapping[str, object]],
    denominators: Mapping[str, Denominator],
    positive_counts: Mapping[str, int],
    analysis_cohort: str,
    cohort_role: str,
    calls_basename: str,
    calls_sha256: str,
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for unit in order:
        denominator = denominators[unit]
        rate, rate_status = percentage(positive_counts[unit], denominator.count)
        output.append(
            {
                "analysis_cohort": analysis_cohort,
                "cohort_role": cohort_role,
                "assembly_unit_id": unit,
                "biological_species": metadata[unit]["biological_species"],
                "haplotype_or_subgenome": metadata[unit]["haplotype_or_subgenome"],
                "assembly_scope": metadata[unit]["assembly_scope"],
                "total_nlr_count": repertoire[unit]["total_nlr_count"],
                "positive_reference_nlr_loss_count": positive_counts[unit],
                "callable_reference_nlr_denominator": denominator.count,
                "positive_reference_nlr_loss_percentage": rate,
                "percentage_status": rate_status,
                "reference_nlr_universe_id": denominator.universe_id,
                "denominator_source_basename": denominator.source_basename,
                "denominator_source_sha256": denominator.source_sha256,
                "repertoire_source_basename": repertoire[unit]["source_basename"],
                "repertoire_source_sha256": repertoire[unit]["source_sha256"],
                "positive_loss_calls_source_basename": calls_basename,
                "positive_loss_calls_source_sha256": calls_sha256,
            }
        )
    return output


def build_species_rows(unit_rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    species_order: list[str] = []
    for row in unit_rows:
        species = str(row["biological_species"])
        if species not in grouped:
            species_order.append(species)
        grouped[species].append(row)
    output: list[dict[str, object]] = []
    for species in species_order:
        rows = grouped[species]
        total = sum(int(row["total_nlr_count"]) for row in rows)
        positive = sum(int(row["positive_reference_nlr_loss_count"]) for row in rows)
        denominator = sum(int(row["callable_reference_nlr_denominator"]) for row in rows)
        rate, rate_status = percentage(positive, denominator)
        output.append(
            {
                "analysis_cohort": rows[0]["analysis_cohort"],
                "cohort_role": rows[0]["cohort_role"],
                "biological_species": species,
                "assembly_unit_count": len(rows),
                "assembly_unit_ids": ";".join(str(row["assembly_unit_id"]) for row in rows),
                "total_nlr_count_sum_across_units": total,
                "positive_reference_nlr_loss_count_sum_across_units": positive,
                "callable_reference_nlr_denominator_sum_across_units": denominator,
                "positive_reference_nlr_loss_percentage_across_unit_comparisons": rate,
                "percentage_status": rate_status,
            }
        )
    return output


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def publish(
    *,
    output_dir: Path,
    unit_rows: Sequence[Mapping[str, object]],
    species_rows: Sequence[Mapping[str, object]],
    validation: Mapping[str, object],
    inputs: Sequence[Path],
) -> None:
    if output_dir.is_symlink():
        raise ValueError(f"output directory must not be a symlink: {output_dir}")
    if output_dir.exists():
        raise FileExistsError(f"output directory already exists: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", suffix=".tmp", dir=output_dir.parent))
    try:
        write_tsv(staging / "nlr_unit_summary.tsv", UNIT_COLUMNS, unit_rows)
        write_tsv(staging / "nlr_species_aggregate.tsv", SPECIES_COLUMNS, species_rows)
        (staging / "validation.json").write_text(
            json.dumps(validation, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        checksums = [
            {
                "role": role,
                "basename": path.name,
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for role, path in zip(
                ("metadata", "repertoire_counts", "positive_loss_calls", "callable_denominators"),
                inputs,
                strict=True,
            )
        ]
        write_tsv(
            staging / "input_checksums.tsv",
            ("role", "basename", "sha256", "bytes"),
            checksums,
        )
        for path in staging.iterdir():
            os.chmod(path, 0o644)
        os.chmod(staging, 0o755)
        os.replace(staging, output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def run(args: argparse.Namespace) -> dict[str, object]:
    cohort = require_text(args.analysis_cohort.strip(), context="--analysis-cohort")
    order, metadata = read_selected_metadata(args.metadata, cohort)
    repertoire = read_repertoire_counts(args.repertoire_counts, metadata)
    denominators, denominator_mode = read_denominators(args.callable_denominators, metadata)
    positive_counts = read_positive_calls(args.positive_loss_calls, metadata, denominators)
    calls_sha = sha256_file(args.positive_loss_calls)
    unit_rows = build_unit_rows(
        order=order,
        metadata=metadata,
        repertoire=repertoire,
        denominators=denominators,
        positive_counts=positive_counts,
        analysis_cohort=cohort,
        cohort_role=args.cohort_role,
        calls_basename=args.positive_loss_calls.name,
        calls_sha256=calls_sha,
    )
    species_rows = build_species_rows(unit_rows)
    validation: dict[str, object] = {
        "status": "pass",
        "analysis_cohort": cohort,
        "cohort_role": args.cohort_role,
        "assembly_unit_count": len(unit_rows),
        "biological_species_count": len(species_rows),
        "denominator_input_mode": denominator_mode,
        "positive_reference_nlr_loss_call_count": sum(positive_counts.values()),
        "callable_reference_nlr_denominator_sum": sum(item.count for item in denominators.values()),
        "undefined_percentage_unit_count": sum(
            1 for row in unit_rows if row["percentage_status"] == "undefined_zero_denominator"
        ),
        "checks": {
            "unique_unit_rows": "pass",
            "exact_included_cohort": "pass",
            "assembly_scope_reconciliation": "pass",
            "denominator_provenance": "pass",
            "positive_calls_not_greater_than_denominator": "pass",
            "single_cohort_only": "pass",
        },
        "species_aggregate_definition": (
            "Descriptive sums across included assembly units; the percentage denominator is "
            "the sum of callable reference-NLR unit comparisons. This table is not a set of "
            "independent biological-species replicates."
        ),
    }
    publish(
        output_dir=args.output_dir,
        unit_rows=unit_rows,
        species_rows=species_rows,
        validation=validation,
        inputs=(
            args.metadata,
            args.repertoire_counts,
            args.positive_loss_calls,
            args.callable_denominators,
        ),
    )
    return validation


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validation = run(args)
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        "Published NLR repertoire/loss summary for "
        f"{validation['assembly_unit_count']} assembly units in cohort "
        f"{validation['analysis_cohort']!r}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
