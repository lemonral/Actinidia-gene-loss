"""Fail-closed species-level phylogenetic generalized least squares.

This module implements the deliberately narrow PGLS model used for the
revision analysis: one row per biological species, a continuity-corrected
logit of a lineage-specific/non-shared positive loss count, and one numeric
predictor.  Technical assembly units, haplotypes, and subgenomes must be
aggregated before this module is called.
"""

from __future__ import annotations

import csv
import ctypes
import errno
import hashlib
import json
import math
import os
import platform
import re
import shutil
import stat
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .io_utils import SchemaError, write_tsv


LOSS_SCOPE = "lineage_specific_nonshared"
ANALYSIS_LEVEL = "biological_species"
MIN_FIT_SPECIES = 5
MIN_PRIMARY_SPECIES = MIN_FIT_SPECIES + 1  # leave-one-species-out is mandatory
PRIMARY_PREDICTOR = "log2_ploidy"
MAX_EXACT_COUNT = 2**63 - 1
MAX_LOG2_PLOIDY = 10.0
MAX_COVARIANCE_CONDITION_NUMBER = 1e12
INPUT_PASS_SCHEMA = "species_pgls_input_pass_v1"
TREE_PASS_SCHEMA = "species_time_tree_pass_v1"
PLOIDY_PASS_SCHEMA = "species_ploidy_ledger_pass_v1"
PGLS_SCHEMA = "species_pgls_v2"
INFERENCE_STATUS = "EXPLORATORY_BLOCKED_DENOMINATOR_AWARE_MODEL_REQUIRED"
DIGEST_PATTERN = re.compile(r"[0-9a-f]{64}")
INPUT_PASS_CHECKS = frozenset(
    {
        "complete_species_gene_grid",
        "all_multiunit_species_use_all_units_positive",
        "positive_calls_are_callable",
        "not_called_loss_is_uncertain",
        "shared_set_is_exact_positive_complete_intersection",
        "shared_removed_from_numerator_and_denominator",
        "uncertain_and_noncallable_excluded_from_denominator",
        "pgls_rows_reconciled_to_species_matrix",
        "ploidy_values_reconciled_to_passed_ledger",
        "exact_biological_species_set",
    }
)
TREE_PASS_CHECKS = frozenset(
    {
        "dated_tree_manifest_pass",
        "rooted_by_accepted_topology",
        "exact_biological_species_tip_set",
        "strictly_bifurcating",
        "ultrametric",
        "finite_nonnegative_branch_lengths",
    }
)
PLOIDY_PASS_CHECKS = frozenset(
    {
        "one_row_per_biological_species",
        "positive_integer_ploidy",
        "log2_ploidy_recalculated_exactly",
        "exact_biological_species_set",
        "source_provenance_complete",
    }
)
SPECIES_LOSS_PASS_CHECKS = frozenset(
    {
        "complete_selected_unit_gene_grid",
        "positive_classification_requires_callable",
        "not_called_loss_treated_as_uncertain",
        "species_status_counts_reconciled",
        "shared_set_reconciled",
        "output_checksums_reconciled",
    }
)
AGGREGATION_POLICY = {
    "species_positive_definition": "positive_complete",
    "multiunit_species_rule": "all_units_positive",
    "positive_requires_every_selected_unit_callable": True,
    "partial_positive_is_not_species_loss": True,
    "partial_positive_denominator_policy": "exclude",
    "not_called_loss_state": "uncertain",
    "uncertain_and_noncallable_denominator_policy": "exclude",
    "shared_definition": "positive_complete_in_every_included_biological_species",
    "shared_removal_policy": "remove_from_numerator_and_denominator",
}
TECHNICAL_UNIT_COLUMNS = frozenset(
    {
        "assembly_unit_id",
        "haplotype_id",
        "subgenome_id",
        "sample_id",
        "terminal_id",
        "technical_unit_id",
    }
)


def _scientific_dependencies():
    try:
        import numpy as np
        from Bio import Phylo
        from scipy import linalg, optimize, stats
    except ImportError as exc:  # pragma: no cover - exercised in minimal installations
        raise RuntimeError(
            "Species PGLS requires NumPy, Biopython, and SciPy. Install the project "
            "with: python -m pip install -e '.[statistics]'"
        ) from exc
    return np, Phylo, linalg, optimize, stats


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class InputSnapshot:
    """Immutable bytes plus identity metadata for one external input."""

    path: Path
    basename: str
    size: int
    sha256: str
    device: int
    inode: int
    mtime_ns: int
    payload: bytes

    def public_binding(self) -> dict[str, str | int]:
        return {"basename": self.basename, "bytes": self.size, "sha256": self.sha256}


def _capture_snapshot(path: str | Path) -> InputSnapshot:
    source = Path(path)
    try:
        with source.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise SchemaError(f"{source.name}: input must resolve to a regular file")
            payload = handle.read()
            after = os.fstat(handle.fileno())
    except OSError:
        raise
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after or len(payload) != before.st_size:
        raise SchemaError(f"{source.name}: input changed while it was being snapshotted")
    try:
        current = source.stat()
    except OSError as exc:
        raise SchemaError(f"{source.name}: input disappeared after snapshot") from exc
    if (current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns) != identity_after:
        raise SchemaError(f"{source.name}: input path changed while it was being snapshotted")
    return InputSnapshot(
        path=source,
        basename=source.name,
        size=len(payload),
        sha256=hashlib.sha256(payload).hexdigest(),
        device=before.st_dev,
        inode=before.st_ino,
        mtime_ns=before.st_mtime_ns,
        payload=payload,
    )


def _require_unchanged(snapshot: InputSnapshot) -> None:
    current = _capture_snapshot(snapshot.path)
    expected = (
        snapshot.device,
        snapshot.inode,
        snapshot.size,
        snapshot.mtime_ns,
        snapshot.sha256,
    )
    observed = (
        current.device,
        current.inode,
        current.size,
        current.mtime_ns,
        current.sha256,
    )
    if observed != expected:
        raise SchemaError(
            f"{snapshot.basename}: input changed after validation; refusing publication"
        )


def _write_snapshot(snapshot: InputSnapshot, destination: Path) -> Path:
    destination.write_bytes(snapshot.payload)
    if _sha256(destination) != snapshot.sha256:
        raise SchemaError(f"{snapshot.basename}: immutable working snapshot checksum mismatch")
    return destination


def _finite_float(value: str, *, column: str, line: int, path: Path) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise SchemaError(f"{path.name}:{line}: {column} is not numeric: {value!r}") from exc
    if not math.isfinite(parsed):
        raise SchemaError(f"{path.name}:{line}: {column} must be finite, found {value!r}")
    return parsed


def _nonnegative_integer(value: str, *, column: str, line: int, path: Path) -> int:
    if not re.fullmatch(r"[0-9]+", value):
        raise SchemaError(
            f"{path.name}:{line}: {column} must be a non-negative integer, found {value!r}"
        )
    canonical = value.lstrip("0") or "0"
    maximum_text = str(MAX_EXACT_COUNT)
    if len(canonical) > len(maximum_text) or (
        len(canonical) == len(maximum_text) and canonical > maximum_text
    ):
        raise SchemaError(
            f"{path.name}:{line}: {column} exceeds the exact signed-64-bit limit "
            f"{MAX_EXACT_COUNT}"
        )
    return int(canonical, 10)


def _read_exact_tsv(path: Path) -> tuple[list[str], list[tuple[dict[str, str], int]]]:
    """Read a rectangular TSV while rejecting ambiguous headers and ragged rows."""

    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", strict=True)
            try:
                raw_fields = next(reader)
            except StopIteration as exc:
                raise SchemaError(f"{path.name}: missing TSV header") from exc
            fields = [field.strip() for field in raw_fields]
            if not fields or any(not field for field in fields):
                raise SchemaError(f"{path.name}: TSV header contains an empty column name")
            if fields != raw_fields:
                raise SchemaError(f"{path.name}: TSV header names may not have surrounding whitespace")
            duplicates = sorted({field for field in fields if fields.count(field) > 1})
            if duplicates:
                raise SchemaError(
                    f"{path.name}: duplicate TSV header column(s): {', '.join(duplicates)}"
                )
            rows: list[tuple[dict[str, str], int]] = []
            for line, values in enumerate(reader, start=2):
                if not values or not any(value.strip() for value in values):
                    continue
                if len(values) != len(fields):
                    raise SchemaError(
                        f"{path.name}:{line}: ragged TSV row has {len(values)} fields; "
                        f"expected exactly {len(fields)}"
                    )
                rows.append(
                    ({field: value.strip() for field, value in zip(fields, values)}, line)
                )
    except csv.Error as exc:
        raise SchemaError(f"{path.name}: malformed TSV: {exc}") from exc
    return fields, rows


@dataclass(frozen=True)
class SpeciesObservation:
    biological_species: str
    positive_loss_count: int
    callable_denominator: int
    predictor: float
    response: float

    @property
    def observed_loss_rate(self) -> float:
        return self.positive_loss_count / self.callable_denominator


@dataclass(frozen=True)
class TreeCovariance:
    species: tuple[str, ...]
    covariance: Any
    root_to_tip: tuple[float, ...]
    tree_height: float
    maximum_tip_height_difference: float
    condition_number: float


@dataclass(frozen=True)
class PGLSFit:
    model_id: str
    species: tuple[str, ...]
    beta: Any
    standard_errors: Any
    confidence_lower: Any
    confidence_upper: Any
    statistics: Any
    p_values: Any
    fitted: Any
    residuals: Any
    marginal_standardized_residuals: Any
    lambda_ml: float
    log_likelihood_ml: float
    residual_variance_ml: float
    residual_variance_unbiased: float
    residual_df: int
    aic: float
    covariance_lambda: Any


def _parse_json_snapshot(snapshot: InputSnapshot) -> dict[str, Any]:
    try:
        decoded = snapshot.payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise SchemaError(f"{snapshot.basename}: JSON input is not UTF-8") from exc
    try:
        value = json.loads(decoded)
    except (json.JSONDecodeError, ValueError) as exc:
        raise SchemaError(f"{snapshot.basename}: invalid JSON: {exc}") from exc
    if not isinstance(value, dict):
        raise SchemaError(f"{snapshot.basename}: JSON top level must be an object")
    return value


def _require_exact_keys(value: Any, expected: set[str] | frozenset[str], *, label: str) -> None:
    if not isinstance(value, dict) or set(value) != set(expected):
        raise SchemaError(f"{label}: unexpected JSON key set")


def _require_binding(value: Any, *, label: str) -> dict[str, str | int]:
    _require_exact_keys(value, {"basename", "bytes", "sha256"}, label=label)
    basename = value["basename"]
    size = value["bytes"]
    digest = value["sha256"]
    if not isinstance(basename, str) or not basename or Path(basename).name != basename:
        raise SchemaError(f"{label}: basename must be one path-free filename")
    if type(size) is not int or size < 0:
        raise SchemaError(f"{label}: bytes must be a non-negative integer")
    if not isinstance(digest, str) or DIGEST_PATTERN.fullmatch(digest) is None:
        raise SchemaError(f"{label}: sha256 must be 64 lowercase hexadecimal characters")
    return {"basename": basename, "bytes": size, "sha256": digest}


def _require_snapshot_binding(value: Any, snapshot: InputSnapshot, *, label: str) -> None:
    if _require_binding(value, label=label) != snapshot.public_binding():
        raise SchemaError(f"{label}: checksum/size/basename binding is not exact")


def _require_all_true_checks(value: Any, expected: frozenset[str], *, label: str) -> None:
    _require_exact_keys(value, expected, label=label)
    if any(type(item) is not bool or not item for item in value.values()):
        raise SchemaError(f"{label}: every required validation check must be true")


def _require_species_list(value: Any, *, label: str) -> list[str]:
    if not isinstance(value, list) or not value:
        raise SchemaError(f"{label}: biological_species must be a non-empty JSON list")
    if any(not isinstance(item, str) or not item.strip() or item != item.strip() for item in value):
        raise SchemaError(f"{label}: biological_species contains an invalid name")
    if len(set(value)) != len(value):
        raise SchemaError(f"{label}: biological_species contains duplicates")
    return value


def _validate_species_loss_manifest(
    snapshot: InputSnapshot, *, expected_species: Sequence[str]
) -> dict[str, Any]:
    report = _parse_json_snapshot(snapshot)
    required = {
        "schema_version",
        "status",
        "definitions",
        "inputs",
        "include_column",
        "assembly_unit_count",
        "biological_species_count",
        "reference_gene_count",
        "expected_unit_matrix_rows",
        "species_gene_matrix_rows",
        "species_gene_status_counts",
        "shared_positive_complete_gene_count",
        "non_shared_positive_call_count",
        "confident_lineage_restricted_gene_count",
        "aggregation_rule_species_counts",
        "species_aggregation",
        "checks",
        "outputs",
    }
    _require_exact_keys(report, required, label=snapshot.basename)
    if report["schema_version"] != "2.0" or report["status"] != "PASS":
        raise SchemaError(f"{snapshot.basename}: species-loss aggregation is not schema 2.0 PASS")
    integer_fields = (
        "assembly_unit_count",
        "biological_species_count",
        "reference_gene_count",
        "expected_unit_matrix_rows",
        "species_gene_matrix_rows",
        "shared_positive_complete_gene_count",
        "non_shared_positive_call_count",
        "confident_lineage_restricted_gene_count",
    )
    if any(type(report[field]) is not int or report[field] < 0 for field in integer_fields):
        raise SchemaError(f"{snapshot.basename}: count fields must be non-negative JSON integers")
    if report["assembly_unit_count"] <= 0 or report["reference_gene_count"] <= 0:
        raise SchemaError(f"{snapshot.basename}: unit and reference-gene counts must be positive")
    if type(report["biological_species_count"]) is not int or report[
        "biological_species_count"
    ] != len(expected_species):
        raise SchemaError(f"{snapshot.basename}: biological species count does not reconcile")
    if report["expected_unit_matrix_rows"] != (
        report["assembly_unit_count"] * report["reference_gene_count"]
    ) or report["species_gene_matrix_rows"] != (
        report["biological_species_count"] * report["reference_gene_count"]
    ):
        raise SchemaError(f"{snapshot.basename}: matrix row arithmetic does not reconcile")
    if not (
        report["shared_positive_complete_gene_count"] <= report["reference_gene_count"]
        and report["confident_lineage_restricted_gene_count"] <= report["reference_gene_count"]
        and report["non_shared_positive_call_count"] <= report["species_gene_matrix_rows"]
    ):
        raise SchemaError(f"{snapshot.basename}: reported loss counts exceed their universes")
    status_counts = report["species_gene_status_counts"]
    if (
        not isinstance(status_counts, dict)
        or any(
            not isinstance(key, str) or type(value) is not int or value < 0
            for key, value in status_counts.items()
        )
        or sum(status_counts.values()) != report["species_gene_matrix_rows"]
    ):
        raise SchemaError(f"{snapshot.basename}: species status counts do not reconcile")
    _require_all_true_checks(
        report["checks"], SPECIES_LOSS_PASS_CHECKS, label=f"{snapshot.basename}.checks"
    )
    aggregation = report["species_aggregation"]
    if not isinstance(aggregation, list) or len(aggregation) != len(expected_species):
        raise SchemaError(f"{snapshot.basename}: species_aggregation has the wrong row count")
    observed: list[str] = []
    observed_unit_total = 0
    for index, row in enumerate(aggregation):
        label = f"{snapshot.basename}.species_aggregation[{index}]"
        _require_exact_keys(
            row,
            {"biological_species", "aggregation_rule", "assembly_unit_count", "assembly_units"},
            label=label,
        )
        species = row["biological_species"]
        units = row["assembly_units"]
        if not isinstance(species, str) or not species:
            raise SchemaError(f"{label}: biological_species is invalid")
        if row["aggregation_rule"] != "all_units_positive":
            raise SchemaError(
                f"{label}: PGLS requires all_units_positive; partial/any-unit loss is invalid"
            )
        if not isinstance(units, list) or not units or any(
            not isinstance(unit, str) or not unit for unit in units
        ):
            raise SchemaError(f"{label}: assembly_units is invalid")
        if len(set(units)) != len(units) or row["assembly_unit_count"] != len(units):
            raise SchemaError(f"{label}: assembly-unit count/list does not reconcile")
        observed.append(species)
        observed_unit_total += len(units)
    if observed != list(expected_species):
        raise SchemaError(
            f"{snapshot.basename}: species aggregation order/set differs from PGLS data"
        )
    if observed_unit_total != report["assembly_unit_count"]:
        raise SchemaError(f"{snapshot.basename}: total assembly-unit count does not reconcile")
    counts = report["aggregation_rule_species_counts"]
    if counts != {"all_units_positive": len(expected_species)}:
        raise SchemaError(
            f"{snapshot.basename}: every PGLS species must use all_units_positive"
        )
    return report


def _validate_ploidy_pass_report(
    snapshot: InputSnapshot, *, expected_species: Sequence[str]
) -> dict[str, Any]:
    report = _parse_json_snapshot(snapshot)
    required = {
        "schema_version",
        "workflow",
        "workflow_version",
        "status",
        "analysis_level",
        "predictor",
        "ploidy_ledger",
        "biological_species",
        "checks",
    }
    _require_exact_keys(report, required, label=snapshot.basename)
    expected_scalars = {
        "schema_version": PLOIDY_PASS_SCHEMA,
        "workflow": "species_ploidy_ledger_validation",
        "workflow_version": "1.0.0",
        "status": "PASS",
        "analysis_level": ANALYSIS_LEVEL,
        "predictor": PRIMARY_PREDICTOR,
    }
    for key, expected in expected_scalars.items():
        if type(report[key]) is not type(expected) or report[key] != expected:
            raise SchemaError(f"{snapshot.basename}: invalid {key}")
    _require_binding(report["ploidy_ledger"], label=f"{snapshot.basename}.ploidy_ledger")
    species = _require_species_list(report["biological_species"], label=snapshot.basename)
    if species != list(expected_species):
        raise SchemaError(f"{snapshot.basename}: ploidy-ledger species order/set mismatch")
    _require_all_true_checks(
        report["checks"], PLOIDY_PASS_CHECKS, label=f"{snapshot.basename}.checks"
    )
    return report


def _validate_input_pass_report(
    snapshot: InputSnapshot,
    *,
    data_snapshot: InputSnapshot,
    species_loss_snapshot: InputSnapshot,
    ploidy_pass_snapshot: InputSnapshot,
    expected_species: Sequence[str],
    predictor_column: str,
) -> dict[str, Any]:
    report = _parse_json_snapshot(snapshot)
    required = {
        "schema_version",
        "workflow",
        "workflow_version",
        "status",
        "analysis_level",
        "loss_scope",
        "predictor",
        "input_data",
        "species_count",
        "biological_species",
        "aggregation_policy",
        "upstream_bindings",
        "checks",
    }
    _require_exact_keys(report, required, label=snapshot.basename)
    expected_scalars = {
        "schema_version": INPUT_PASS_SCHEMA,
        "workflow": "species_pgls_input_builder",
        "workflow_version": "1.0.0",
        "status": "PASS",
        "analysis_level": ANALYSIS_LEVEL,
        "loss_scope": LOSS_SCOPE,
        "predictor": PRIMARY_PREDICTOR,
        "species_count": len(expected_species),
    }
    for key, expected in expected_scalars.items():
        if type(report[key]) is not type(expected) or report[key] != expected:
            raise SchemaError(f"{snapshot.basename}: invalid {key}")
    if predictor_column != PRIMARY_PREDICTOR:
        raise SchemaError(
            f"species PGLS primary predictor must be exactly {PRIMARY_PREDICTOR!r}"
        )
    _require_snapshot_binding(
        report["input_data"], data_snapshot, label=f"{snapshot.basename}.input_data"
    )
    species = _require_species_list(report["biological_species"], label=snapshot.basename)
    if species != list(expected_species):
        raise SchemaError(f"{snapshot.basename}: PGLS input species order/set mismatch")
    if report["aggregation_policy"] != AGGREGATION_POLICY:
        raise SchemaError(f"{snapshot.basename}: aggregation policy is not the frozen PGLS policy")
    bindings = report["upstream_bindings"]
    _require_exact_keys(
        bindings,
        {
            "species_loss_manifest",
            "species_gene_matrix",
            "shared_positive_complete_gene_set",
            "ploidy_ledger",
            "ploidy_ledger_pass_report",
        },
        label=f"{snapshot.basename}.upstream_bindings",
    )
    _require_snapshot_binding(
        bindings["species_loss_manifest"],
        species_loss_snapshot,
        label=f"{snapshot.basename}.upstream_bindings.species_loss_manifest",
    )
    _require_snapshot_binding(
        bindings["ploidy_ledger_pass_report"],
        ploidy_pass_snapshot,
        label=f"{snapshot.basename}.upstream_bindings.ploidy_ledger_pass_report",
    )
    for role in ("species_gene_matrix", "shared_positive_complete_gene_set", "ploidy_ledger"):
        _require_binding(bindings[role], label=f"{snapshot.basename}.upstream_bindings.{role}")
    _require_all_true_checks(
        report["checks"], INPUT_PASS_CHECKS, label=f"{snapshot.basename}.checks"
    )
    return report


def _validate_tree_pass_report(
    snapshot: InputSnapshot,
    *,
    tree_snapshot: InputSnapshot,
    expected_species: Sequence[str],
) -> dict[str, Any]:
    report = _parse_json_snapshot(snapshot)
    required = {
        "schema_version",
        "workflow",
        "workflow_version",
        "status",
        "analysis_level",
        "tree",
        "source_dating_manifest",
        "biological_species",
        "root_semantics",
        "branch_length_units",
        "checks",
    }
    _require_exact_keys(report, required, label=snapshot.basename)
    expected_scalars = {
        "schema_version": TREE_PASS_SCHEMA,
        "workflow": "species_time_tree_validation",
        "workflow_version": "1.0.0",
        "status": "PASS",
        "analysis_level": ANALYSIS_LEVEL,
        "root_semantics": "accepted_biological_species_mrca",
        "branch_length_units": "million_years",
    }
    for key, expected in expected_scalars.items():
        if type(report[key]) is not type(expected) or report[key] != expected:
            raise SchemaError(f"{snapshot.basename}: invalid {key}")
    _require_snapshot_binding(report["tree"], tree_snapshot, label=f"{snapshot.basename}.tree")
    _require_binding(
        report["source_dating_manifest"], label=f"{snapshot.basename}.source_dating_manifest"
    )
    species = _require_species_list(report["biological_species"], label=snapshot.basename)
    if species != list(expected_species):
        raise SchemaError(f"{snapshot.basename}: time-tree species order/set mismatch")
    _require_all_true_checks(
        report["checks"], TREE_PASS_CHECKS, label=f"{snapshot.basename}.checks"
    )
    return report


def read_species_data(
    path: str | Path,
    *,
    predictor_column: str,
    species_column: str = "biological_species",
    count_column: str = "lineage_specific_nonshared_positive_loss_count",
    denominator_column: str = "callable_denominator",
    scope_column: str = "loss_scope",
    level_column: str = "analysis_level",
) -> list[SpeciesObservation]:
    """Read and validate the exact biological-species analysis table."""

    source = Path(path)
    required = {
        species_column,
        count_column,
        denominator_column,
        predictor_column,
        scope_column,
        level_column,
    }
    fields, raw_rows = _read_exact_tsv(source)
    missing = sorted(required.difference(fields))
    if missing:
        raise SchemaError(f"{source.name}: missing required columns: {', '.join(missing)}")
    forbidden = sorted(TECHNICAL_UNIT_COLUMNS.intersection(fields))
    if forbidden:
        raise SchemaError(
            f"{source.name}: technical-unit columns are forbidden in species PGLS: "
            f"{', '.join(forbidden)}; aggregate haplotypes/subgenomes first"
        )

    if len(raw_rows) < MIN_PRIMARY_SPECIES:
        raise SchemaError(
            f"{source.name}: species PGLS requires at least {MIN_PRIMARY_SPECIES} biological "
            "species so every leave-one-species-out fit retains at least "
            f"{MIN_FIT_SPECIES}; found {len(raw_rows)}"
        )

    observations: list[SpeciesObservation] = []
    seen: dict[str, int] = {}
    for row, line in raw_rows:
        species = row[species_column]
        if not species:
            raise SchemaError(f"{source.name}:{line}: empty {species_column}")
        if any(character in species for character in ("\x00", "\r", "\n", "\t")):
            raise SchemaError(f"{source.name}:{line}: {species_column} contains a control character")
        if species in seen:
            raise SchemaError(
                f"{source.name}:{line}: duplicate biological species {species!r}; first used at "
                f"line {seen[species]}. Haplotype/subgenome rows are pseudoreplicates here"
            )
        seen[species] = line
        if row[scope_column] != LOSS_SCOPE:
            raise SchemaError(
                f"{source.name}:{line}: {scope_column} must be exactly {LOSS_SCOPE!r}; "
                f"found {row[scope_column]!r}"
            )
        if row[level_column] != ANALYSIS_LEVEL:
            raise SchemaError(
                f"{source.name}:{line}: {level_column} must be exactly {ANALYSIS_LEVEL!r}; "
                f"found {row[level_column]!r}"
            )
        count = _nonnegative_integer(
            row[count_column], column=count_column, line=line, path=source
        )
        denominator = _nonnegative_integer(
            row[denominator_column], column=denominator_column, line=line, path=source
        )
        if denominator <= 0:
            raise SchemaError(f"{source.name}:{line}: {denominator_column} must be positive")
        if count > denominator:
            raise SchemaError(
                f"{source.name}:{line}: {count_column} ({count}) exceeds "
                f"{denominator_column} ({denominator})"
            )
        if predictor_column == PRIMARY_PREDICTOR and re.fullmatch(
            r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?", row[predictor_column]
        ) is None:
            raise SchemaError(
                f"{source.name}:{line}: {PRIMARY_PREDICTOR} must be a canonical "
                f"non-negative decimal, found {row[predictor_column]!r}"
            )
        predictor = _finite_float(
            row[predictor_column], column=predictor_column, line=line, path=source
        )
        if predictor_column == PRIMARY_PREDICTOR and not (
            0.0 <= predictor <= MAX_LOG2_PLOIDY
        ):
            raise SchemaError(
                f"{source.name}:{line}: {PRIMARY_PREDICTOR} must be between 0 and "
                f"{MAX_LOG2_PLOIDY:g}; found {predictor!r}"
            )
        # Haldane-Anscombe correction.  The difference-of-logs form avoids
        # floating-point overflow/underflow and preserves arbitrarily large
        # exactly parsed integer counts.
        response = math.log(2 * count + 1) - math.log(2 * (denominator - count) + 1)
        observations.append(
            SpeciesObservation(species, count, denominator, predictor, response)
        )

    predictor_values = [row.predictor for row in observations]
    response_values = [row.response for row in observations]
    if max(predictor_values) == min(predictor_values):
        raise SchemaError(f"{source.name}: predictor {predictor_column!r} has no variation")
    if max(response_values) == min(response_values):
        raise SchemaError(f"{source.name}: continuity-corrected logit response has no variation")
    return observations


def build_brownian_covariance(
    tree_path: str | Path,
    species: Sequence[str],
    *,
    ultrametric_tolerance: float = 1e-6,
) -> TreeCovariance:
    """Return the Brownian shared-path covariance in the supplied species order.

    The top-level node in the Newick representation is treated as the
    biological root.  The tree must be ultrametric and every non-root branch
    must have an explicit, finite, non-negative length.
    """

    if ultrametric_tolerance < 0 or not math.isfinite(ultrametric_tolerance):
        raise SchemaError("ultrametric_tolerance must be finite and non-negative")
    np, Phylo, _, _, _ = _scientific_dependencies()
    source = Path(tree_path)
    with source.open(encoding="utf-8") as handle:
        trees = list(Phylo.parse(handle, "newick"))
    if len(trees) != 1:
        raise SchemaError(f"{source.name}: expected exactly one Newick tree, found {len(trees)}")
    tree = trees[0]
    internal_clades = [clade for clade in tree.find_clades() if not clade.is_terminal()]
    nonbinary = [clade for clade in internal_clades if len(clade.clades) != 2]
    if nonbinary:
        raise SchemaError(
            f"{source.name}: accepted PGLS time tree must be strictly bifurcating; "
            f"found {len(nonbinary)} non-binary internal node(s)"
        )
    if tree.root.branch_length not in (None, 0, 0.0):
        raise SchemaError(
            f"{source.name}: root/stem branch length is ambiguous; omit it or set it to zero"
        )

    terminals = tree.get_terminals()
    names = [str(tip.name or "").strip() for tip in terminals]
    if any(not name for name in names):
        raise SchemaError(f"{source.name}: every tree tip must have a non-empty name")
    if len(set(names)) != len(names):
        duplicates = sorted({name for name in names if names.count(name) > 1})
        raise SchemaError(f"{source.name}: duplicate tree tip name(s): {', '.join(duplicates)}")
    expected, observed = set(species), set(names)
    if expected != observed or len(species) != len(names):
        missing = sorted(expected.difference(observed))
        extra = sorted(observed.difference(expected))
        raise SchemaError(
            f"{source.name}: exact tree/data tip reconciliation failed; "
            f"missing_from_tree={','.join(missing) or 'none'}; "
            f"absent_from_data={','.join(extra) or 'none'}"
        )

    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        length = clade.branch_length
        if length is None or not math.isfinite(float(length)) or float(length) < 0:
            raise SchemaError(
                f"{source.name}: every non-root branch needs a finite, non-negative length"
            )

    tip_by_name = {str(tip.name).strip(): tip for tip in terminals}
    heights = [float(tree.distance(tree.root, tip_by_name[name])) for name in species]
    if min(heights) <= 0:
        raise SchemaError(f"{source.name}: all root-to-tip distances must be positive")
    height_difference = max(heights) - min(heights)
    tolerance_bp = ultrametric_tolerance * max(heights)
    if height_difference > tolerance_bp:
        raise SchemaError(
            f"{source.name}: time tree is not ultrametric; maximum root-to-tip difference "
            f"{height_difference:.12g} exceeds tolerance {tolerance_bp:.12g}"
        )

    count = len(species)
    covariance = np.zeros((count, count), dtype=float)
    for i, first in enumerate(species):
        covariance[i, i] = heights[i]
        for j in range(i):
            ancestor = tree.common_ancestor(tip_by_name[first], tip_by_name[species[j]])
            shared = float(tree.distance(tree.root, ancestor))
            covariance[i, j] = shared
            covariance[j, i] = shared
    if not np.all(np.isfinite(covariance)) or not np.allclose(
        covariance, covariance.T, rtol=0.0, atol=1e-12
    ):
        raise SchemaError(f"{source.name}: Brownian covariance is non-finite or asymmetric")
    off_diagonal = covariance[~np.eye(count, dtype=bool)]
    if not np.any(off_diagonal > np.finfo(float).eps * max(heights)):
        raise SchemaError(
            f"{source.name}: Pagel lambda is unidentifiable because the accepted tree has "
            "no positive shared-path covariance"
        )
    try:
        np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise SchemaError(
            f"{source.name}: Brownian covariance is singular or not positive definite"
        ) from exc
    condition_number = float(np.linalg.cond(covariance))
    if not math.isfinite(condition_number) or condition_number > MAX_COVARIANCE_CONDITION_NUMBER:
        raise SchemaError(
            f"{source.name}: Brownian covariance condition number {condition_number:.12g} "
            f"exceeds the fail-closed limit {MAX_COVARIANCE_CONDITION_NUMBER:.12g}"
        )
    return TreeCovariance(
        tuple(species),
        covariance,
        tuple(heights),
        float(sum(heights) / len(heights)),
        float(height_difference),
        condition_number,
    )


def _pagel_covariance(base_covariance: Any, lambda_value: float, np: Any) -> Any:
    transformed = np.array(base_covariance, dtype=float, copy=True)
    off_diagonal = ~np.eye(transformed.shape[0], dtype=bool)
    transformed[off_diagonal] *= lambda_value
    return transformed


def _profile_fit(y: Any, x: Any, base_covariance: Any, lambda_value: float, np: Any, linalg: Any):
    covariance = _pagel_covariance(base_covariance, lambda_value, np)
    try:
        factor, lower = linalg.cho_factor(covariance, lower=True, check_finite=True)
        inverse_y = linalg.cho_solve((factor, lower), y, check_finite=True)
        inverse_x = linalg.cho_solve((factor, lower), x, check_finite=True)
    except linalg.LinAlgError as exc:
        raise SchemaError("Pagel-lambda covariance is singular or not positive definite") from exc
    information = x.T @ inverse_x
    if np.linalg.matrix_rank(information) != x.shape[1]:
        raise SchemaError("PGLS design/information matrix is singular")
    try:
        beta = np.linalg.solve(information, x.T @ inverse_y)
    except np.linalg.LinAlgError as exc:
        raise SchemaError("PGLS coefficient system is singular") from exc
    residual = y - x @ beta
    inverse_residual = linalg.cho_solve((factor, lower), residual, check_finite=True)
    rss = float(residual.T @ inverse_residual)
    if not math.isfinite(rss) or rss <= np.finfo(float).eps * max(1.0, float(y.T @ y)):
        raise SchemaError("PGLS residual variance is zero or non-finite")
    log_determinant = 2.0 * float(np.log(np.diag(factor)).sum())
    n = len(y)
    sigma2_ml = rss / n
    log_likelihood = -0.5 * (
        n * (math.log(2.0 * math.pi) + 1.0 + math.log(sigma2_ml)) + log_determinant
    )
    return log_likelihood, beta, residual, information, covariance, rss


def fit_pgls(
    observations: Sequence[SpeciesObservation],
    base_covariance: Any,
    *,
    model_id: str,
) -> PGLSFit:
    """Fit intercept + one predictor with ML Pagel lambda in [0, 1]."""

    np, _, linalg, optimize, stats = _scientific_dependencies()
    n = len(observations)
    if n < MIN_FIT_SPECIES:
        raise SchemaError(
            f"{model_id}: PGLS requires at least {MIN_FIT_SPECIES} biological species; found {n}"
        )
    y = np.asarray([row.response for row in observations], dtype=float)
    predictor = np.asarray([row.predictor for row in observations], dtype=float)
    x = np.column_stack((np.ones(n, dtype=float), predictor))
    if np.linalg.matrix_rank(x) != 2:
        raise SchemaError(f"{model_id}: predictor is constant or design matrix is singular")
    if np.ptp(y) == 0:
        raise SchemaError(f"{model_id}: response has no variation")
    covariance = np.asarray(base_covariance, dtype=float)
    if covariance.shape != (n, n):
        raise SchemaError(
            f"{model_id}: covariance has shape {covariance.shape}; expected {(n, n)}"
        )
    if not np.all(np.isfinite(covariance)) or not np.allclose(
        covariance, covariance.T, rtol=0.0, atol=1e-12
    ):
        raise SchemaError(f"{model_id}: Brownian covariance is non-finite or asymmetric")
    try:
        np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise SchemaError(f"{model_id}: Brownian covariance is singular or non-positive") from exc
    condition_number = float(np.linalg.cond(covariance))
    if not math.isfinite(condition_number) or condition_number > MAX_COVARIANCE_CONDITION_NUMBER:
        raise SchemaError(
            f"{model_id}: Brownian covariance is numerically ill-conditioned "
            f"({condition_number:.12g})"
        )

    def objective(value: float) -> float:
        try:
            return -float(_profile_fit(y, x, covariance, value, np, linalg)[0])
        except SchemaError:
            return math.inf

    identifiability_values: list[float] = []
    for candidate in (0.0, 0.25, 0.5, 0.75, 1.0):
        try:
            identifiability_values.append(
                float(_profile_fit(y, x, covariance, candidate, np, linalg)[0])
            )
        except SchemaError:
            continue
    if len(identifiability_values) < 2:
        raise SchemaError(f"{model_id}: Pagel lambda profile could not be identified")
    likelihood_span = max(identifiability_values) - min(identifiability_values)
    likelihood_scale = max(1.0, max(abs(value) for value in identifiability_values))
    if likelihood_span <= 1e-10 * likelihood_scale:
        raise SchemaError(
            f"{model_id}: Pagel lambda is unidentifiable because its profile likelihood is flat"
        )

    optimized = optimize.minimize_scalar(
        objective,
        bounds=(0.0, 1.0),
        method="bounded",
        options={"xatol": 1e-9, "maxiter": 1000},
    )
    if not optimized.success or not math.isfinite(float(optimized.fun)):
        raise SchemaError(f"{model_id}: Pagel lambda optimization failed: {optimized.message}")
    candidates = [0.0, float(optimized.x), 1.0]
    evaluated: list[tuple[float, tuple[Any, ...]]] = []
    for candidate in candidates:
        try:
            evaluated.append(
                (candidate, _profile_fit(y, x, covariance, candidate, np, linalg))
            )
        except SchemaError:
            continue
    if not evaluated:
        raise SchemaError(f"{model_id}: no positive-definite Pagel-lambda covariance was fit")
    lambda_ml, profile = max(evaluated, key=lambda item: item[1][0])
    log_likelihood, beta, residual, information, covariance_lambda, rss = profile
    fitted = x @ beta
    degrees_freedom = n - x.shape[1]
    if degrees_freedom <= 0:
        raise SchemaError(f"{model_id}: non-positive residual degrees of freedom")
    sigma2_ml = rss / n
    sigma2_unbiased = rss / degrees_freedom
    try:
        coefficient_covariance = sigma2_unbiased * np.linalg.inv(information)
    except np.linalg.LinAlgError as exc:
        raise SchemaError(f"{model_id}: coefficient covariance is singular") from exc
    diagonal = np.diag(coefficient_covariance)
    if np.any(diagonal <= 0) or not np.all(np.isfinite(diagonal)):
        raise SchemaError(f"{model_id}: coefficient covariance is non-positive or non-finite")
    standard_errors = np.sqrt(diagonal)
    statistics = beta / standard_errors
    p_values = 2.0 * stats.t.sf(np.abs(statistics), degrees_freedom)
    critical = float(stats.t.ppf(0.975, degrees_freedom))
    lower = beta - critical * standard_errors
    upper = beta + critical * standard_errors
    marginal_scale = np.sqrt(sigma2_unbiased * np.diag(covariance_lambda))
    if np.any(marginal_scale <= 0):
        raise SchemaError(f"{model_id}: residual standardization scale is non-positive")
    return PGLSFit(
        model_id=model_id,
        species=tuple(row.biological_species for row in observations),
        beta=beta,
        standard_errors=standard_errors,
        confidence_lower=lower,
        confidence_upper=upper,
        statistics=statistics,
        p_values=p_values,
        fitted=fitted,
        residuals=residual,
        marginal_standardized_residuals=residual / marginal_scale,
        lambda_ml=float(lambda_ml),
        log_likelihood_ml=float(log_likelihood),
        residual_variance_ml=float(sigma2_ml),
        residual_variance_unbiased=float(sigma2_unbiased),
        residual_df=degrees_freedom,
        aic=float(2 * 4 - 2 * log_likelihood),
        covariance_lambda=covariance_lambda,
    )


def parse_named_sensitivities(values: Iterable[str]) -> dict[str, tuple[str, ...]]:
    """Parse repeatable ``NAME=Species[,Species]`` exclusion definitions."""

    parsed: dict[str, tuple[str, ...]] = {}
    for value in values:
        if "=" not in value:
            raise SchemaError(
                f"invalid sensitivity {value!r}; expected NAME=Species[,Species]"
            )
        label, species_text = value.split("=", 1)
        label = label.strip()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", label):
            raise SchemaError(
                f"invalid sensitivity name {label!r}; use letters, digits, dot, underscore, or dash"
            )
        if label in parsed or label == "primary":
            raise SchemaError(f"duplicate or reserved sensitivity name {label!r}")
        excluded = tuple(item.strip() for item in species_text.split(",") if item.strip())
        if not excluded:
            raise SchemaError(f"sensitivity {label!r} has no excluded species")
        if len(set(excluded)) != len(excluded):
            raise SchemaError(f"sensitivity {label!r} repeats an excluded species")
        parsed[label] = excluded
    return parsed


def _format(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def _rename_directory_no_replace(source: Path, destination: Path) -> None:
    """Atomically publish a directory without replacing a racing destination."""

    library = ctypes.CDLL(None, use_errno=True)
    source_bytes = os.fsencode(source)
    destination_bytes = os.fsencode(destination)
    result: int | None = None
    if hasattr(library, "renameat2"):
        renameat2 = library.renameat2
        renameat2.argtypes = (
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        )
        renameat2.restype = ctypes.c_int
        result = renameat2(-100, source_bytes, -100, destination_bytes, 1)
    elif hasattr(library, "renamex_np"):
        renamex = library.renamex_np
        renamex.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint)
        renamex.restype = ctypes.c_int
        result = renamex(source_bytes, destination_bytes, 0x00000004)
    if result is None:
        raise SchemaError(
            "platform lacks an atomic no-replace directory rename; refusing publication"
        )
    if result == 0:
        return
    error_number = ctypes.get_errno()
    if error_number in {errno.EEXIST, errno.ENOTEMPTY}:
        raise SchemaError(
            f"output appeared during publication; refusing overwrite: {destination.name}"
        )
    raise SchemaError(
        f"atomic no-replace publication failed for {destination.name}: "
        f"{os.strerror(error_number)}"
    )


def _model_summary_row(
    fit: PGLSFit, predictor_column: str, excluded_species: Sequence[str]
) -> dict[str, Any]:
    return {
        "model_id": fit.model_id,
        "excluded_species": ";".join(excluded_species),
        "analysis_level": ANALYSIS_LEVEL,
        "loss_scope": LOSS_SCOPE,
        "response": "log((positive_loss_count+0.5)/(callable_denominator-positive_loss_count+0.5))",
        "predictor": predictor_column,
        "n_species": len(fit.species),
        "lambda_ml": _format(fit.lambda_ml),
        "log_likelihood_ml": _format(fit.log_likelihood_ml),
        "residual_variance_ml": _format(fit.residual_variance_ml),
        "residual_variance_unbiased_for_inference": _format(fit.residual_variance_unbiased),
        "residual_degrees_of_freedom": fit.residual_df,
        "aic_ml": _format(fit.aic),
        "inference_status": INFERENCE_STATUS,
    }


def _coefficient_rows(fit: PGLSFit, predictor_column: str) -> list[dict[str, Any]]:
    rows = []
    for index, term in enumerate(("intercept", predictor_column)):
        rows.append(
            {
                "model_id": fit.model_id,
                "term": term,
                "estimate": _format(float(fit.beta[index])),
                "standard_error": _format(float(fit.standard_errors[index])),
                "confidence_lower_95": _format(float(fit.confidence_lower[index])),
                "confidence_upper_95": _format(float(fit.confidence_upper[index])),
                "t_statistic": _format(float(fit.statistics[index])),
                "degrees_of_freedom": fit.residual_df,
                "p_value_two_sided": _format(float(fit.p_values[index])),
                "inference_status": INFERENCE_STATUS,
            }
        )
    return rows


def _leave_one_out_rows(
    observations: Sequence[SpeciesObservation], base_covariance: Any, predictor_column: str
) -> list[dict[str, Any]]:
    np, _, linalg, _, _ = _scientific_dependencies()
    rows: list[dict[str, Any]] = []
    for omitted_index, omitted in enumerate(observations):
        retained_indices = [index for index in range(len(observations)) if index != omitted_index]
        retained = [observations[index] for index in retained_indices]
        retained_covariance = base_covariance[np.ix_(retained_indices, retained_indices)]
        fit = fit_pgls(
            retained,
            retained_covariance,
            model_id=f"leave_one_out:{omitted.biological_species}",
        )
        cross = base_covariance[omitted_index, retained_indices].astype(float, copy=True)
        cross *= fit.lambda_ml
        try:
            factor, lower = linalg.cho_factor(
                fit.covariance_lambda, lower=True, check_finite=True
            )
            conditional = float(
                cross @ linalg.cho_solve((factor, lower), fit.residuals, check_finite=True)
            )
            conditional_variance_factor = float(
                base_covariance[omitted_index, omitted_index]
                - cross
                @ linalg.cho_solve((factor, lower), cross, check_finite=True)
            )
        except linalg.LinAlgError as exc:
            raise SchemaError(
                f"leave_one_out:{omitted.biological_species}: predictive covariance is singular"
            ) from exc
        if conditional_variance_factor <= 0 or not math.isfinite(conditional_variance_factor):
            raise SchemaError(
                f"leave_one_out:{omitted.biological_species}: predictive variance is non-positive"
            )
        marginal_prediction = float(fit.beta[0] + fit.beta[1] * omitted.predictor)
        prediction = marginal_prediction + conditional
        prediction_se = math.sqrt(fit.residual_variance_unbiased * conditional_variance_factor)
        rows.append(
            {
                "omitted_biological_species": omitted.biological_species,
                "n_species_fit": len(retained),
                "predictor": predictor_column,
                "predictor_value": _format(omitted.predictor),
                "observed_response": _format(omitted.response),
                "predicted_response_phylogenetic": _format(prediction),
                "prediction_standard_error_residual_component": _format(prediction_se),
                "prediction_error_observed_minus_predicted": _format(omitted.response - prediction),
                "lambda_ml": _format(fit.lambda_ml),
                "log_likelihood_ml": _format(fit.log_likelihood_ml),
                "predictor_estimate": _format(float(fit.beta[1])),
                "predictor_standard_error": _format(float(fit.standard_errors[1])),
                "predictor_confidence_lower_95": _format(float(fit.confidence_lower[1])),
                "predictor_confidence_upper_95": _format(float(fit.confidence_upper[1])),
                "predictor_p_value_two_sided": _format(float(fit.p_values[1])),
                "inference_status": INFERENCE_STATUS,
            }
        )
    return rows


def run_species_pgls(
    *,
    data_path: str | Path,
    tree_path: str | Path,
    input_pass_report_path: str | Path,
    species_loss_manifest_path: str | Path,
    ploidy_ledger_pass_report_path: str | Path,
    tree_pass_report_path: str | Path,
    output_dir: str | Path,
    predictor_column: str,
    sensitivities: Mapping[str, Sequence[str]] | None = None,
    ultrametric_tolerance: float = 1e-6,
    species_column: str = "biological_species",
    count_column: str = "lineage_specific_nonshared_positive_loss_count",
    denominator_column: str = "callable_denominator",
    scope_column: str = "loss_scope",
    level_column: str = "analysis_level",
) -> dict[str, Path]:
    """Validate, fit, and atomically publish a complete species-PGLS bundle."""

    np, _, _, _, _ = _scientific_dependencies()
    output = Path(output_dir)
    if os.path.lexists(output):
        raise SchemaError(f"output directory already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    snapshots = {
        "data": _capture_snapshot(data_path),
        "tree": _capture_snapshot(tree_path),
        "input_pass_report": _capture_snapshot(input_pass_report_path),
        "species_loss_manifest": _capture_snapshot(species_loss_manifest_path),
        "ploidy_ledger_pass_report": _capture_snapshot(ploidy_ledger_pass_report_path),
        "tree_pass_report": _capture_snapshot(tree_pass_report_path),
    }
    with tempfile.TemporaryDirectory(
        prefix=f".{output.name}.validated-inputs.", dir=output.parent
    ) as snapshot_directory:
        snapshot_root = Path(snapshot_directory)
        frozen_data = _write_snapshot(snapshots["data"], snapshot_root / "analysis_data.tsv")
        frozen_tree = _write_snapshot(snapshots["tree"], snapshot_root / "time_tree.nwk")
        observations = read_species_data(
            frozen_data,
            predictor_column=predictor_column,
            species_column=species_column,
            count_column=count_column,
            denominator_column=denominator_column,
            scope_column=scope_column,
            level_column=level_column,
        )
        species = [row.biological_species for row in observations]
        _validate_species_loss_manifest(
            snapshots["species_loss_manifest"], expected_species=species
        )
        _validate_ploidy_pass_report(
            snapshots["ploidy_ledger_pass_report"], expected_species=species
        )
        _validate_input_pass_report(
            snapshots["input_pass_report"],
            data_snapshot=snapshots["data"],
            species_loss_snapshot=snapshots["species_loss_manifest"],
            ploidy_pass_snapshot=snapshots["ploidy_ledger_pass_report"],
            expected_species=species,
            predictor_column=predictor_column,
        )
        _validate_tree_pass_report(
            snapshots["tree_pass_report"],
            tree_snapshot=snapshots["tree"],
            expected_species=species,
        )
        tree = build_brownian_covariance(
            frozen_tree, species, ultrametric_tolerance=ultrametric_tolerance
        )
    base_covariance = np.asarray(tree.covariance, dtype=float)
    primary = fit_pgls(observations, base_covariance, model_id="primary")

    declared = dict(sensitivities or {})
    known = set(species)
    model_fits: list[tuple[PGLSFit, tuple[str, ...]]] = [(primary, ())]
    for label, excluded_values in declared.items():
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", label) or label == "primary":
            raise SchemaError(f"invalid or reserved sensitivity name {label!r}")
        excluded = tuple(excluded_values)
        if not excluded or len(set(excluded)) != len(excluded):
            raise SchemaError(f"sensitivity {label!r} must list unique excluded species")
        unknown = sorted(set(excluded).difference(known))
        if unknown:
            raise SchemaError(
                f"sensitivity {label!r} names species absent from data: {', '.join(unknown)}"
            )
        retained_indices = [index for index, name in enumerate(species) if name not in excluded]
        retained = [observations[index] for index in retained_indices]
        sensitivity_covariance = base_covariance[np.ix_(retained_indices, retained_indices)]
        model_fits.append(
            (
                fit_pgls(
                    retained,
                    sensitivity_covariance,
                    model_id=f"sensitivity:{label}",
                ),
                excluded,
            )
        )

    loo_rows = _leave_one_out_rows(observations, base_covariance, predictor_column)
    data_rows = [
        {
            "biological_species": row.biological_species,
            "analysis_level": ANALYSIS_LEVEL,
            "loss_scope": LOSS_SCOPE,
            "lineage_specific_nonshared_positive_loss_count": row.positive_loss_count,
            "callable_denominator": row.callable_denominator,
            "observed_loss_rate": _format(row.observed_loss_rate),
            "predictor": predictor_column,
            "predictor_value": _format(row.predictor),
            "continuity_corrected_logit_response": _format(row.response),
            "root_to_tip_distance": _format(tree.root_to_tip[index]),
        }
        for index, row in enumerate(observations)
    ]
    summary_rows = [
        _model_summary_row(fit, predictor_column, excluded)
        for fit, excluded in model_fits
    ]
    coefficient_rows = [
        row for fit, _ in model_fits for row in _coefficient_rows(fit, predictor_column)
    ]
    fitted_rows = [
        {
            "model_id": "primary",
            "biological_species": observation.biological_species,
            "observed_response": _format(observation.response),
            "fitted_response": _format(float(primary.fitted[index])),
            "raw_residual_observed_minus_fitted": _format(float(primary.residuals[index])),
            "marginal_standardized_residual": _format(
                float(primary.marginal_standardized_residuals[index])
            ),
            "positive_loss_count": observation.positive_loss_count,
            "callable_denominator": observation.callable_denominator,
            "observed_loss_rate": _format(observation.observed_loss_rate),
            "predictor": predictor_column,
            "predictor_value": _format(observation.predictor),
        }
        for index, observation in enumerate(observations)
    ]
    sensitivity_rows = [
        {
            "sensitivity_name": fit.model_id.removeprefix("sensitivity:"),
            "excluded_species": ";".join(excluded),
            "n_species": len(fit.species),
            "lambda_ml": _format(fit.lambda_ml),
            "log_likelihood_ml": _format(fit.log_likelihood_ml),
            "predictor": predictor_column,
            "predictor_estimate": _format(float(fit.beta[1])),
            "predictor_standard_error": _format(float(fit.standard_errors[1])),
            "predictor_confidence_lower_95": _format(float(fit.confidence_lower[1])),
            "predictor_confidence_upper_95": _format(float(fit.confidence_upper[1])),
            "predictor_p_value_two_sided": _format(float(fit.p_values[1])),
            "inference_status": INFERENCE_STATUS,
        }
        for fit, excluded in model_fits[1:]
    ]

    temporary = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        files: dict[str, Path] = {}
        files["analysis_data"] = write_tsv(
            temporary / "analysis_data.tsv",
            data_rows,
            [
                "biological_species",
                "analysis_level",
                "loss_scope",
                "lineage_specific_nonshared_positive_loss_count",
                "callable_denominator",
                "observed_loss_rate",
                "predictor",
                "predictor_value",
                "continuity_corrected_logit_response",
                "root_to_tip_distance",
            ],
        )
        files["model_summary"] = write_tsv(
            temporary / "model_summary.tsv",
            summary_rows,
            [
                "model_id",
                "excluded_species",
                "analysis_level",
                "loss_scope",
                "response",
                "predictor",
                "n_species",
                "lambda_ml",
                "log_likelihood_ml",
                "residual_variance_ml",
                "residual_variance_unbiased_for_inference",
                "residual_degrees_of_freedom",
                "aic_ml",
                "inference_status",
            ],
        )
        files["coefficients"] = write_tsv(
            temporary / "model_coefficients.tsv",
            coefficient_rows,
            [
                "model_id",
                "term",
                "estimate",
                "standard_error",
                "confidence_lower_95",
                "confidence_upper_95",
                "t_statistic",
                "degrees_of_freedom",
                "p_value_two_sided",
                "inference_status",
            ],
        )
        files["fitted_residuals"] = write_tsv(
            temporary / "fitted_residuals.tsv",
            fitted_rows,
            [
                "model_id",
                "biological_species",
                "observed_response",
                "fitted_response",
                "raw_residual_observed_minus_fitted",
                "marginal_standardized_residual",
                "positive_loss_count",
                "callable_denominator",
                "observed_loss_rate",
                "predictor",
                "predictor_value",
            ],
        )
        files["leave_one_species_out"] = write_tsv(
            temporary / "leave_one_species_out.tsv",
            loo_rows,
            [
                "omitted_biological_species",
                "n_species_fit",
                "predictor",
                "predictor_value",
                "observed_response",
                "predicted_response_phylogenetic",
                "prediction_standard_error_residual_component",
                "prediction_error_observed_minus_predicted",
                "lambda_ml",
                "log_likelihood_ml",
                "predictor_estimate",
                "predictor_standard_error",
                "predictor_confidence_lower_95",
                "predictor_confidence_upper_95",
                "predictor_p_value_two_sided",
                "inference_status",
            ],
        )
        files["named_exclusion_sensitivities"] = write_tsv(
            temporary / "named_exclusion_sensitivities.tsv",
            sensitivity_rows,
            [
                "sensitivity_name",
                "excluded_species",
                "n_species",
                "lambda_ml",
                "log_likelihood_ml",
                "predictor",
                "predictor_estimate",
                "predictor_standard_error",
                "predictor_confidence_lower_95",
                "predictor_confidence_upper_95",
                "predictor_p_value_two_sided",
                "inference_status",
            ],
        )
        files["publication_gate"] = write_tsv(
            temporary / "publication_gate.tsv",
            [
                {
                    "publication_status": "BLOCKED_DENOMINATOR_AWARE_MODEL_REQUIRED",
                    "ordinary_pgls_role": "exploratory_phylogenetic_trait_sensitivity_only",
                    "blocking_reason": (
                        "ordinary Gaussian PGLS does not model the unequal count precision "
                        "implied by species-specific callable denominators"
                    ),
                    "required_resolution": (
                        "predeclare and validate a denominator-aware phylogenetic count model; "
                        "do not treat an ad hoc weighted PGLS as equivalent"
                    ),
                }
            ],
            [
                "publication_status",
                "ordinary_pgls_role",
                "blocking_reason",
                "required_resolution",
            ],
        )
        try:
            import Bio
            import scipy
        except ImportError as exc:  # pragma: no cover - dependencies already checked
            raise RuntimeError("scientific dependency disappeared during PGLS run") from exc
        manifest = {
            "schema_version": PGLS_SCHEMA,
            "status": "COMPLETE_EXPLORATORY_BLOCKED_FOR_PUBLICATION",
            "publication_gate": "BLOCKED_DENOMINATOR_AWARE_MODEL_REQUIRED",
            "analysis_level": ANALYSIS_LEVEL,
            "loss_scope": LOSS_SCOPE,
            "input_data": snapshots["data"].public_binding(),
            "input_time_tree": snapshots["tree"].public_binding(),
            "input_pass_report": snapshots["input_pass_report"].public_binding(),
            "species_loss_manifest": snapshots["species_loss_manifest"].public_binding(),
            "ploidy_ledger_pass_report": snapshots[
                "ploidy_ledger_pass_report"
            ].public_binding(),
            "time_tree_pass_report": snapshots["tree_pass_report"].public_binding(),
            "species_count": len(observations),
            "predictor": predictor_column,
            "source_columns": {
                "species": species_column,
                "positive_loss_count": count_column,
                "callable_denominator": denominator_column,
                "scope": scope_column,
                "analysis_level": level_column,
            },
            "aggregation_policy": AGGREGATION_POLICY,
            "response_formula": "log((L+0.5)/(D-L+0.5))",
            "brownian_covariance": "shared root-to-MRCA branch length",
            "pagel_lambda": "ML, constrained to [0,1], off-diagonal covariance scaled",
            "coefficient_inference": (
                "exploratory conditional two-sided t inference using RSS/(n-2); "
                "lambda-estimation uncertainty is not propagated"
            ),
            "denominator_handling": (
                "callable denominator enters the transformed response but not an observation-"
                "specific variance model; therefore ordinary PGLS is not publication-ready"
            ),
            "tree_height": tree.tree_height,
            "maximum_tip_height_difference": tree.maximum_tip_height_difference,
            "brownian_covariance_condition_number": tree.condition_number,
            "ultrametric_relative_tolerance": ultrametric_tolerance,
            "named_exclusion_sensitivities": {
                label: list(excluded) for label, excluded in declared.items()
            },
            "runtime": {
                "python": platform.python_version(),
                "implementation": platform.python_implementation(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "biopython": Bio.__version__,
                "operating_system": platform.system(),
                "machine": platform.machine(),
            },
        }
        manifest_path = temporary / "analysis_manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        files["manifest"] = manifest_path
        checksum_rows = [
            {"relative_path": path.name, "sha256": _sha256(path)}
            for path in sorted(files.values(), key=lambda item: item.name)
        ]
        files["checksums"] = write_tsv(
            temporary / "checksums.sha256.tsv",
            checksum_rows,
            ["relative_path", "sha256"],
        )
        for snapshot in snapshots.values():
            _require_unchanged(snapshot)
        _rename_directory_no_replace(temporary, output)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return {label: output / path.name for label, path in files.items()}
