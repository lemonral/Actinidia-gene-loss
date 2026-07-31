"""Fail-closed chromosome-homology assignment from frozen score matrices.

This module does not run nucleotide alignments or JCVI.  It consumes four
complete, precomputed 29 x 29 matrices (nucleotide and JCVI evidence against
HY4A and HY4P), validates their arithmetic and identifier closure, solves each
matrix independently with a deterministic Hungarian assignment, and publishes
final chromosome labels only when every configured gate passes.
"""

from __future__ import annotations

import csv
import ctypes
import errno
import hashlib
import io
import json
import math
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Mapping, Sequence

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 only
    import tomli as tomllib  # type: ignore[no-redef]

from .io_utils import natural_key
from .chromosome_provenance import (
    ChromosomeProvenanceError,
    FileSnapshot,
    MatrixProvenance,
    PROVENANCE_SCHEMA_VERSION,
    capture_snapshot,
    parse_exact_int,
    read_reference_registries,
    read_target_registry,
    reject_duplicate_snapshots,
    validate_matrix_provenance,
)


WORKFLOW_VERSION = "1.0.0"
MATRIX_SCHEMA_VERSION = "1.0.0"
EXPECTED_CHROMOSOMES = 29
COORDINATE_REFERENCE_ID = "act_chinensis_hongyang_v4_hy4a"
CONFIRMATION_REFERENCE_ID = "act_chinensis_hongyang_v4_hy4p"
MINIMAP2_VERSION = "2.28-r1209"
MINIMAP2_COMMAND_TEMPLATE = (
    "minimap2",
    "-x",
    "asm5",
    "--secondary=no",
    "-c",
    "--cs=long",
    "{reference_fasta}",
    "{query_fasta}",
)
REQUIRED_PAF_TAGS = ("tp:A:P", "de:f", "cg:Z", "cs:Z")
PUBLIC_FILE_MODE = 0o644
PUBLIC_DIRECTORY_MODE = 0o755
SAFE_UNIT_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
SAFE_CHROMOSOME_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:+-]{0,127}$")
GIT_COMMIT_RE = re.compile(r"^(?:[0-9a-f]{40}|[0-9a-f]{64})$")
DECIMAL_RE = re.compile(
    r"^-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?$"
)

NUCLEOTIDE_COLUMNS = (
    "query_chromosome",
    "reference_chromosome",
    "canonical_chromosome",
    "score",
    "query_covered_bp",
    "query_length_bp",
    "query_coverage",
    "reference_covered_bp",
    "reference_length_bp",
    "reference_coverage",
    "reciprocal_coverage",
    "matching_bases",
    "weighted_divergence",
    "orientation",
)

JCVI_COLUMNS = (
    "query_chromosome",
    "reference_chromosome",
    "canonical_chromosome",
    "score",
    "query_anchored_genes",
    "query_eligible_genes",
    "query_gene_coverage",
    "reference_anchored_genes",
    "reference_eligible_genes",
    "reference_gene_coverage",
    "unique_anchor_pairs",
)

EVIDENCE_COLUMNS = (
    "query_chromosome",
    "assigned_reference_chromosome",
    "assigned_canonical_chromosome",
    "assigned_score",
    "assigned_reciprocal_coverage",
    "assigned_reciprocal_gene_coverage",
    "assigned_unique_anchor_pairs",
    "assigned_orientation",
    "row_best_reference_chromosome",
    "row_best_score",
    "row_second_score",
    "row_top_second_ratio",
    "row_normalized_margin",
    "row_unique_best",
    "column_best_query_chromosome",
    "column_best_score",
    "column_second_score",
    "column_top_second_ratio",
    "column_normalized_margin",
    "column_unique_best",
    "row_and_column_reciprocal_best",
    "separation_gate",
    "coverage_gate",
    "absolute_support_gate",
    "matrix_gate",
    "failure_reasons",
)

COMBINED_COLUMNS = (
    "query_chromosome",
    "nucleotide_hy4a_reference",
    "nucleotide_hy4a_canonical",
    "jcvi_hy4a_reference",
    "jcvi_hy4a_canonical",
    "hy4a_nucleotide_jcvi_agreement",
    "nucleotide_hy4p_reference",
    "nucleotide_hy4p_canonical",
    "jcvi_hy4p_reference",
    "jcvi_hy4p_canonical",
    "hy4p_nucleotide_jcvi_agreement",
    "hy4a_hy4p_label_agreement",
    "candidate_canonical_chromosome",
    "all_matrix_gates",
    "final_chromosome",
    "status",
    "failure_reasons",
)

FINAL_MAP_COLUMNS = (
    "query_chromosome",
    "final_chromosome",
    "coordinate_reference",
    "confirmation_reference",
    "status",
)

SUMMARY_COLUMNS = (
    "assembly_unit_id",
    "target_scope_id",
    "trusted_repository_commit",
    "status",
    "publication_gate",
    "failure_states",
    "expected_query_chromosomes",
    "observed_query_chromosomes",
    "expected_reference_chromosomes_per_reference",
    "nucleotide_jcvi_agreement_count_hy4a",
    "nucleotide_jcvi_agreement_count_hy4p",
    "hy4a_hy4p_label_agreement_count",
    "rows_passing_all_matrix_gates",
    "final_map_row_count",
)

FAILURE_PRIORITY = (
    "NON_BIJECTIVE",
    "CONFLICT_NUCLEOTIDE_JCVI",
    "HY4A_HY4P_DISAGREEMENT",
    "LOW_RECIPROCAL_COVERAGE",
    "LOW_JCVI_ABSOLUTE_SUPPORT",
    "AMBIGUOUS_RECIPROCAL_BEST",
    "AMBIGUOUS_SEPARATION",
)


class ChromosomeAssignmentError(RuntimeError):
    """Raised when score matrices or assignment policy are malformed."""


@dataclass(frozen=True)
class HomologyPolicy:
    """Frozen chromosome-homology gates loaded from the project TOML file."""

    coordinate_reference: str
    confirmation_reference: str
    matrix_schema_version: str
    matrix_provenance_schema_version: str
    reference_asset_registry_sha256: str
    reference_chromosome_map_registry_sha256: str
    expected_query_chromosomes: int
    expected_reference_chromosomes: int
    assignment_method: str
    minimap2_version: str
    minimap2_command_template: tuple[str, ...]
    minimap2_primary_only: bool
    minimum_mapq: int
    minimum_alignment_block_bp: int
    maximum_de: float
    coverage_arithmetic: str
    minimum_top_second_ratio: float
    minimum_normalized_score_margin: float
    minimum_assigned_reciprocal_nucleotide_coverage: float
    minimum_unique_anchor_pairs: int
    minimum_reciprocal_gene_coverage: float
    minimum_assigned_jcvi_score: float
    require_row_and_column_reciprocal_best: bool
    require_nucleotide_and_jcvi_assignment_agreement: bool
    require_hy4a_hy4p_label_agreement: bool
    reverse_complement_allowed: bool
    nucleotide_score_formula: str
    reciprocal_nucleotide_coverage_formula: str
    jcvi_score_formula: str
    jcvi_anchor_counting: str
    jcvi_aligner: str
    jcvi_database_type: str
    jcvi_cscore: float
    jcvi_tandem_nmax: int
    jcvi_maximum_gene_distance: int
    jcvi_minimum_anchor_block_size: int
    jcvi_coverage_anchor_source: str
    arithmetic_tolerance: float
    failure_policy: str


@dataclass(frozen=True)
class MatrixCell:
    """One validated cell in a long-form score matrix."""

    query: str
    reference: str
    canonical: str
    score: float
    raw: Mapping[str, str]


@dataclass(frozen=True)
class ScoreMatrix:
    """Validated complete matrix and its reference-to-canonical map."""

    role: str
    kind: str
    source: Path
    columns: tuple[str, ...]
    queries: tuple[str, ...]
    references: tuple[str, ...]
    canonical_by_reference: Mapping[str, str]
    cells: Mapping[tuple[str, str], MatrixCell]
    normalized_rows: tuple[Mapping[str, str], ...]


@dataclass(frozen=True)
class AssignmentEvidence:
    """One matrix's globally assigned edge and local gate evidence."""

    query: str
    reference: str
    canonical: str
    score: float
    reciprocal_coverage: float | None
    reciprocal_gene_coverage: float | None
    unique_anchor_pairs: int | None
    orientation: str
    row_best_reference: str
    row_best_score: float
    row_second_score: float
    row_ratio: float | None
    row_margin: float
    row_unique_best: bool
    column_best_query: str
    column_best_score: float
    column_second_score: float
    column_ratio: float | None
    column_margin: float
    column_unique_best: bool
    reciprocal_best: bool
    separation_gate: bool
    coverage_gate: bool
    absolute_support_gate: bool
    matrix_gate: bool
    failure_reasons: tuple[str, ...]


@dataclass(frozen=True)
class ChromosomeAssignmentResult:
    """Outcome of an atomically published assignment audit."""

    output_dir: Path
    status: str
    publication_gate: str
    failure_states: tuple[str, ...]
    final_map_row_count: int


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
    except OSError as error:
        raise ChromosomeAssignmentError(f"Cannot read input {path.name}: {error}") from error
    return digest.hexdigest()


def _finite_decimal(
    value: str,
    *,
    column: str,
    path: Path,
    line: int,
    nonnegative: bool = False,
) -> Decimal:
    if len(value) > 128 or not DECIMAL_RE.fullmatch(value):
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: {column} is not an exact decimal: {value!r}"
        )
    try:
        result = Decimal(value)
    except InvalidOperation as error:
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: {column} is not numeric: {value!r}"
        ) from error
    if not result.is_finite():
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: {column} must be finite: {value!r}"
        )
    if nonnegative and value.startswith("-"):
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: {column} must use a nonnegative decimal spelling"
        )
    try:
        converted = float(result)
    except (InvalidOperation, OverflowError, ValueError) as error:
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: {column} lies outside the supported finite numeric range"
        ) from error
    if not math.isfinite(converted) or (result != 0 and converted == 0.0):
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: {column} lies outside the supported finite numeric range"
        )
    return result


def _integer(value: str, *, column: str, path: Path, line: int, positive: bool = False) -> int:
    try:
        return parse_exact_int(
            value,
            label=f"{path.name}:{line}:{column}",
            positive=positive,
        )
    except ChromosomeProvenanceError as error:
        raise ChromosomeAssignmentError(str(error)) from error


def _format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.12g}"


def _format_bool(value: bool) -> str:
    return "true" if value else "false"


def _decimal_harmonic(first: Decimal, second: Decimal) -> Decimal:
    if first <= 0 or second <= 0:
        return Decimal(0)
    return Decimal(2) * first * second / (first + second)


def _decimal_close(observed: Decimal, expected: Decimal, tolerance: Decimal) -> bool:
    difference = abs(observed - expected)
    relative = tolerance * max(abs(observed), abs(expected))
    return difference <= max(tolerance, relative)


def _close(observed: float, expected: float, tolerance: float) -> bool:
    return math.isclose(observed, expected, rel_tol=tolerance, abs_tol=tolerance)


def _canonical_labels(count: int) -> tuple[str, ...]:
    return tuple(f"Chr{index:02d}" for index in range(1, count + 1))


def _total_key(value: str) -> tuple[list[object], str]:
    """Natural order with the original ID as a deterministic collision breaker."""

    return natural_key(value), value


def _require_policy_types(
    section: Mapping[str, Any], contract: Mapping[str, type | tuple[type, ...]], label: str
) -> None:
    missing = sorted(set(contract).difference(section))
    if missing:
        raise ChromosomeAssignmentError(f"{label} lacks: {', '.join(missing)}")
    for key, expected_type in contract.items():
        value = section[key]
        if isinstance(value, bool) and expected_type != bool:
            raise ChromosomeAssignmentError(f"{label}.{key} has the wrong type")
        if not isinstance(value, expected_type):
            raise ChromosomeAssignmentError(f"{label}.{key} has the wrong type")


def load_homology_policy(path: str | Path | FileSnapshot) -> HomologyPolicy:
    """Load strong-typed production policy from one immutable TOML snapshot."""

    snapshot = path if isinstance(path, FileSnapshot) else capture_snapshot(path)
    try:
        data = tomllib.loads(
            snapshot.data.decode("utf-8", errors="strict"), parse_float=Decimal
        )
        section = data["chromosome_homology"]
        jcvi = data["assembly_qc"]["jcvi"]
    except (
        UnicodeError,
        KeyError,
        InvalidOperation,
        ValueError,
        tomllib.TOMLDecodeError,
    ) as error:
        raise ChromosomeAssignmentError(
            f"Cannot load frozen policy from {snapshot.basename}: {error}"
        ) from error
    if not isinstance(section, dict) or not isinstance(jcvi, dict):
        raise ChromosomeAssignmentError("Frozen policy sections must be TOML tables")
    homology_contract: dict[str, type | tuple[type, ...]] = {
        "coordinate_reference": str,
        "independent_confirmation_reference": str,
        "matrix_schema_version": str,
        "matrix_provenance_schema_version": str,
        "reference_asset_registry_sha256": str,
        "reference_chromosome_map_registry_sha256": str,
        "expected_query_chromosomes": int,
        "expected_reference_chromosomes": int,
        "minimap2_version": str,
        "minimap2_command_template": list,
        "minimap2_primary_only": bool,
        "minimum_mapq": int,
        "minimum_alignment_block_bp": int,
        "maximum_de": (int, float, Decimal),
        "coverage_arithmetic": str,
        "assignment_method": str,
        "nucleotide_score_formula": str,
        "reciprocal_nucleotide_coverage_formula": str,
        "jcvi_score_formula": str,
        "jcvi_anchor_counting": str,
        "arithmetic_tolerance": (int, float, Decimal),
        "minimum_top_second_ratio": (int, float, Decimal),
        "minimum_normalized_score_margin": (int, float, Decimal),
        "minimum_assigned_reciprocal_nucleotide_coverage": (int, float, Decimal),
        "minimum_unique_anchor_pairs": int,
        "minimum_reciprocal_gene_coverage": (int, float, Decimal),
        "minimum_assigned_jcvi_score": (int, float, Decimal),
        "require_row_and_column_reciprocal_best": bool,
        "require_nucleotide_and_jcvi_assignment_agreement": bool,
        "require_hy4a_hy4p_label_agreement": bool,
        "reverse_complement_allowed": bool,
        "failure_policy": str,
    }
    jcvi_contract: dict[str, type | tuple[type, ...]] = {
        "aligner": str,
        "database_type": str,
        "cscore": (int, float, Decimal),
        "tandem_nmax": int,
        "maximum_gene_distance": int,
        "minimum_anchor_block_size": int,
        "coverage_anchor_source": str,
    }
    _require_policy_types(section, homology_contract, "[chromosome_homology]")
    _require_policy_types(jcvi, jcvi_contract, "[assembly_qc.jcvi]")
    policy_decimals = {
        "maximum_de": Decimal(section["maximum_de"]),
        "minimum_top_second_ratio": Decimal(section["minimum_top_second_ratio"]),
        "minimum_normalized_score_margin": Decimal(
            section["minimum_normalized_score_margin"]
        ),
        "minimum_assigned_reciprocal_nucleotide_coverage": Decimal(
            section["minimum_assigned_reciprocal_nucleotide_coverage"]
        ),
        "minimum_reciprocal_gene_coverage": Decimal(
            section["minimum_reciprocal_gene_coverage"]
        ),
        "minimum_assigned_jcvi_score": Decimal(
            section["minimum_assigned_jcvi_score"]
        ),
        "jcvi_cscore": Decimal(jcvi["cscore"]),
        "arithmetic_tolerance": Decimal(section["arithmetic_tolerance"]),
    }
    for label, value in policy_decimals.items():
        if not value.is_finite():
            raise ChromosomeAssignmentError(f"Frozen numeric policy {label} must be finite")

    def finite_policy_float(label: str) -> float:
        try:
            value = float(policy_decimals[label])
        except (InvalidOperation, OverflowError, ValueError) as error:
            raise ChromosomeAssignmentError(
                f"Frozen numeric policy {label} cannot be represented safely"
            ) from error
        if not math.isfinite(value):
            raise ChromosomeAssignmentError(f"Frozen numeric policy {label} must be finite")
        return value

    policy = HomologyPolicy(
        coordinate_reference=section["coordinate_reference"],
        confirmation_reference=section["independent_confirmation_reference"],
        matrix_schema_version=section["matrix_schema_version"],
        matrix_provenance_schema_version=section["matrix_provenance_schema_version"],
        reference_asset_registry_sha256=section["reference_asset_registry_sha256"],
        reference_chromosome_map_registry_sha256=section[
            "reference_chromosome_map_registry_sha256"
        ],
        expected_query_chromosomes=section["expected_query_chromosomes"],
        expected_reference_chromosomes=section["expected_reference_chromosomes"],
        assignment_method=section["assignment_method"],
        minimap2_version=section["minimap2_version"],
        minimap2_command_template=tuple(section["minimap2_command_template"]),
        minimap2_primary_only=section["minimap2_primary_only"],
        minimum_mapq=section["minimum_mapq"],
        minimum_alignment_block_bp=section["minimum_alignment_block_bp"],
        maximum_de=finite_policy_float("maximum_de"),
        coverage_arithmetic=section["coverage_arithmetic"],
        minimum_top_second_ratio=finite_policy_float("minimum_top_second_ratio"),
        minimum_normalized_score_margin=finite_policy_float(
            "minimum_normalized_score_margin"
        ),
        minimum_assigned_reciprocal_nucleotide_coverage=finite_policy_float(
            "minimum_assigned_reciprocal_nucleotide_coverage"
        ),
        minimum_unique_anchor_pairs=section["minimum_unique_anchor_pairs"],
        minimum_reciprocal_gene_coverage=finite_policy_float(
            "minimum_reciprocal_gene_coverage"
        ),
        minimum_assigned_jcvi_score=finite_policy_float(
            "minimum_assigned_jcvi_score"
        ),
        require_row_and_column_reciprocal_best=section[
            "require_row_and_column_reciprocal_best"
        ],
        require_nucleotide_and_jcvi_assignment_agreement=section[
            "require_nucleotide_and_jcvi_assignment_agreement"
        ],
        require_hy4a_hy4p_label_agreement=section[
            "require_hy4a_hy4p_label_agreement"
        ],
        reverse_complement_allowed=section["reverse_complement_allowed"],
        nucleotide_score_formula=section["nucleotide_score_formula"],
        reciprocal_nucleotide_coverage_formula=section[
            "reciprocal_nucleotide_coverage_formula"
        ],
        jcvi_score_formula=section["jcvi_score_formula"],
        jcvi_anchor_counting=section["jcvi_anchor_counting"],
        jcvi_aligner=jcvi["aligner"],
        jcvi_database_type=jcvi["database_type"],
        jcvi_cscore=finite_policy_float("jcvi_cscore"),
        jcvi_tandem_nmax=jcvi["tandem_nmax"],
        jcvi_maximum_gene_distance=jcvi["maximum_gene_distance"],
        jcvi_minimum_anchor_block_size=jcvi["minimum_anchor_block_size"],
        jcvi_coverage_anchor_source=jcvi["coverage_anchor_source"],
        arithmetic_tolerance=finite_policy_float("arithmetic_tolerance"),
        failure_policy=section["failure_policy"],
    )
    if policy.coordinate_reference != COORDINATE_REFERENCE_ID or (
        policy.confirmation_reference != CONFIRMATION_REFERENCE_ID
    ):
        raise ChromosomeAssignmentError(
            "Production policy requires the exact frozen HY4A coordinate and HY4P confirmation IDs"
        )
    for label, value in (
        ("reference_asset_registry_sha256", policy.reference_asset_registry_sha256),
        (
            "reference_chromosome_map_registry_sha256",
            policy.reference_chromosome_map_registry_sha256,
        ),
    ):
        if not re.fullmatch(r"[0-9a-f]{64}", value):
            raise ChromosomeAssignmentError(f"{label} must be a lowercase SHA-256")
    if policy.expected_query_chromosomes != EXPECTED_CHROMOSOMES or (
        policy.expected_reference_chromosomes != EXPECTED_CHROMOSOMES
    ):
        raise ChromosomeAssignmentError("Production policy requires exact 29 x 29 matrices")
    if policy.matrix_schema_version != MATRIX_SCHEMA_VERSION:
        raise ChromosomeAssignmentError("Unsupported matrix_schema_version")
    if policy.matrix_provenance_schema_version != PROVENANCE_SCHEMA_VERSION:
        raise ChromosomeAssignmentError("Unsupported matrix_provenance_schema_version")
    if policy.assignment_method != "global one-to-one Hungarian assignment":
        raise ChromosomeAssignmentError("Hungarian assignment policy is required")
    if policy.minimap2_version != MINIMAP2_VERSION:
        raise ChromosomeAssignmentError(
            f"minimap2_version must remain frozen at {MINIMAP2_VERSION}"
        )
    if policy.minimap2_command_template != MINIMAP2_COMMAND_TEMPLATE:
        raise ChromosomeAssignmentError(
            "minimap2_command_template must equal the frozen path-independent argv"
        )
    if not policy.minimap2_primary_only:
        raise ChromosomeAssignmentError("minimap2_primary_only must be true")
    if policy.minimum_mapq != 20 or policy.minimum_alignment_block_bp != 10000:
        raise ChromosomeAssignmentError("Frozen nucleotide gates require MAPQ 20 and 10 kb blocks")
    if policy_decimals["maximum_de"] != Decimal("0.15"):
        raise ChromosomeAssignmentError("maximum_de must remain frozen at 0.15")
    if policy.coverage_arithmetic != (
        "interval union in both query and reference directions; never raw PAF-row sums"
    ):
        raise ChromosomeAssignmentError("Bidirectional interval-union coverage policy is required")
    if policy.failure_policy != (
        "retain provisional PubChr labels and stop; never force a final Chr label"
    ):
        raise ChromosomeAssignmentError("The fail-closed provisional-label policy is required")
    if policy_decimals["minimum_top_second_ratio"] != Decimal("1.5"):
        raise ChromosomeAssignmentError("minimum_top_second_ratio must remain frozen at 1.5")
    if policy_decimals["minimum_normalized_score_margin"] != Decimal("0.10"):
        raise ChromosomeAssignmentError("minimum_normalized_score_margin must remain 0.10")
    if policy.minimum_unique_anchor_pairs != 30:
        raise ChromosomeAssignmentError("minimum_unique_anchor_pairs must remain 30")
    if policy_decimals["minimum_reciprocal_gene_coverage"] != Decimal("0.05"):
        raise ChromosomeAssignmentError("minimum_reciprocal_gene_coverage must remain 0.05")
    if policy_decimals["minimum_assigned_jcvi_score"] != Decimal("0.05"):
        raise ChromosomeAssignmentError("minimum_assigned_jcvi_score must remain 0.05")
    if policy_decimals[
        "minimum_assigned_reciprocal_nucleotide_coverage"
    ] != Decimal("0.05"):
        raise ChromosomeAssignmentError(
            "minimum_assigned_reciprocal_nucleotide_coverage must remain 0.05"
        )
    for label in (
        "minimum_normalized_score_margin",
        "minimum_assigned_reciprocal_nucleotide_coverage",
        "minimum_reciprocal_gene_coverage",
        "minimum_assigned_jcvi_score",
    ):
        if not Decimal(0) <= policy_decimals[label] <= Decimal(1):
            raise ChromosomeAssignmentError(f"{label} must lie in [0,1]")
    if policy_decimals["arithmetic_tolerance"] != Decimal("1e-9"):
        raise ChromosomeAssignmentError("arithmetic_tolerance must remain frozen at 1e-9")
    if not (
        policy.require_row_and_column_reciprocal_best
        and policy.require_nucleotide_and_jcvi_assignment_agreement
        and policy.require_hy4a_hy4p_label_agreement
    ):
        raise ChromosomeAssignmentError("All reciprocal-best and agreement gates must be enabled")
    if policy.reverse_complement_allowed:
        raise ChromosomeAssignmentError("This label-only workflow cannot reverse-complement")
    if policy.nucleotide_score_formula != (
        "harmonic_mean(query_coverage,reference_coverage)*(1-weighted_divergence)"
    ):
        raise ChromosomeAssignmentError("Unsupported nucleotide_score_formula")
    if policy.reciprocal_nucleotide_coverage_formula != (
        "min(query_coverage,reference_coverage)"
    ):
        raise ChromosomeAssignmentError("Unsupported reciprocal coverage formula")
    if policy.jcvi_score_formula != (
        "harmonic_mean(query_gene_coverage,reference_gene_coverage)"
    ):
        raise ChromosomeAssignmentError("Unsupported jcvi_score_formula")
    if policy.jcvi_anchor_counting != (
        "unique anchored genes and unique anchor pairs after exact chromosome-ID reconciliation"
    ):
        raise ChromosomeAssignmentError("Unsupported JCVI anchor-counting policy")
    expected_jcvi = {
        "aligner": (policy.jcvi_aligner, "LAST"),
        "database_type": (policy.jcvi_database_type, "protein"),
        "cscore": (policy_decimals["jcvi_cscore"], Decimal("0.7")),
        "tandem_nmax": (policy.jcvi_tandem_nmax, 10),
        "maximum_gene_distance": (policy.jcvi_maximum_gene_distance, 20),
        "minimum_anchor_block_size": (policy.jcvi_minimum_anchor_block_size, 4),
        "coverage_anchor_source": (policy.jcvi_coverage_anchor_source, "raw JCVI anchors"),
    }
    for label, (observed, expected) in expected_jcvi.items():
        if observed != expected:
            raise ChromosomeAssignmentError(
                f"JCVI {label} must remain frozen at {expected!r}"
            )
    return policy


def _read_exact_rows(
    snapshot: FileSnapshot, columns: Sequence[str]
) -> list[tuple[int, dict[str, str]]]:
    path = snapshot.path
    try:
        text = snapshot.data.decode("utf-8", errors="strict")
        with io.StringIO(text, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            observed = tuple(reader.fieldnames or ())
            if observed != tuple(columns):
                raise ChromosomeAssignmentError(
                    f"{path.name}: exact schema required; expected {list(columns)}, "
                    f"found {list(observed)}"
                )
            rows = []
            for line, row in enumerate(reader, start=2):
                if None in row:
                    raise ChromosomeAssignmentError(
                        f"{path.name}:{line}: row has more fields than the exact schema"
                    )
                cleaned = {key: (value or "").strip() for key, value in row.items()}
                if not any(cleaned.values()):
                    continue
                if any(value == "" for value in cleaned.values()):
                    missing = [key for key, value in cleaned.items() if value == ""]
                    raise ChromosomeAssignmentError(
                        f"{path.name}:{line}: empty required values: {','.join(missing)}"
                    )
                rows.append((line, cleaned))
    except (OSError, UnicodeError, csv.Error) as error:
        raise ChromosomeAssignmentError(f"Cannot read matrix {path.name}: {error}") from error
    if not rows:
        raise ChromosomeAssignmentError(f"{path.name}: matrix contains no data rows")
    return rows


def _validate_identifier(value: str, *, column: str, path: Path, line: int) -> None:
    if not SAFE_CHROMOSOME_ID.fullmatch(value):
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: unsafe or empty {column}: {value!r}"
        )


def _validate_common_cell(
    row: Mapping[str, str], *, path: Path, line: int, expected_count: int
) -> tuple[str, str, str, float]:
    query = row["query_chromosome"]
    reference = row["reference_chromosome"]
    canonical = row["canonical_chromosome"]
    for column, value in (
        ("query_chromosome", query),
        ("reference_chromosome", reference),
        ("canonical_chromosome", canonical),
    ):
        _validate_identifier(value, column=column, path=path, line=line)
    if canonical not in _canonical_labels(expected_count):
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: canonical_chromosome must be Chr01..Chr{expected_count:02d}"
        )
    score_decimal = _finite_decimal(
        row["score"], column="score", path=path, line=line, nonnegative=True
    )
    if score_decimal > 1:
        raise ChromosomeAssignmentError(f"{path.name}:{line}: score must lie in [0,1]")
    score = float(score_decimal)
    return query, reference, canonical, score


def _validate_nucleotide_row(
    row: Mapping[str, str], *, path: Path, line: int, policy: HomologyPolicy
) -> None:
    query_covered = _integer(
        row["query_covered_bp"], column="query_covered_bp", path=path, line=line
    )
    query_length = _integer(
        row["query_length_bp"], column="query_length_bp", path=path, line=line, positive=True
    )
    reference_covered = _integer(
        row["reference_covered_bp"], column="reference_covered_bp", path=path, line=line
    )
    reference_length = _integer(
        row["reference_length_bp"],
        column="reference_length_bp",
        path=path,
        line=line,
        positive=True,
    )
    matching_bases = _integer(
        row["matching_bases"], column="matching_bases", path=path, line=line
    )
    if query_covered > query_length or reference_covered > reference_length:
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: covered bp cannot exceed chromosome length"
        )
    query_fraction = _finite_decimal(
        row["query_coverage"],
        column="query_coverage",
        path=path,
        line=line,
        nonnegative=True,
    )
    reference_fraction = _finite_decimal(
        row["reference_coverage"],
        column="reference_coverage",
        path=path,
        line=line,
        nonnegative=True,
    )
    reciprocal = _finite_decimal(
        row["reciprocal_coverage"],
        column="reciprocal_coverage",
        path=path,
        line=line,
        nonnegative=True,
    )
    divergence = _finite_decimal(
        row["weighted_divergence"],
        column="weighted_divergence",
        path=path,
        line=line,
        nonnegative=True,
    )
    for column, value in (
        ("query_coverage", query_fraction),
        ("reference_coverage", reference_fraction),
        ("reciprocal_coverage", reciprocal),
        ("weighted_divergence", divergence),
    ):
        if value > 1:
            raise ChromosomeAssignmentError(f"{path.name}:{line}: {column} must lie in [0,1]")
    if row["orientation"] not in {"+", "-", "mixed", "none"}:
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: orientation must be +, -, mixed, or none"
        )
    if (query_covered == 0) != (reference_covered == 0):
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: bidirectional covered bp must be zero or non-zero together"
        )
    if (query_covered == 0 or reference_covered == 0) and (
        matching_bases != 0 or row["orientation"] != "none"
    ):
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: zero coverage requires zero matching bases and orientation none"
        )
    tolerance = Decimal(str(policy.arithmetic_tolerance))
    if query_covered == 0 and not _decimal_close(divergence, Decimal(1), tolerance):
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: zero coverage requires weighted_divergence = 1"
        )
    if query_covered > 0 and reference_covered > 0 and matching_bases == 0:
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: non-zero reciprocal coverage requires matching_bases > 0"
        )
    if query_covered > 0 and row["orientation"] == "none":
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: non-zero reciprocal coverage requires an orientation"
        )
    if query_covered > 0 and divergence > Decimal(str(policy.maximum_de)):
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: weighted_divergence exceeds frozen maximum_de"
        )
    expected_query = Decimal(query_covered) / Decimal(query_length)
    expected_reference = Decimal(reference_covered) / Decimal(reference_length)
    expected_reciprocal = min(query_fraction, reference_fraction)
    expected_score = _decimal_harmonic(query_fraction, reference_fraction) * (
        Decimal(1) - divergence
    )
    score = _finite_decimal(
        row["score"], column="score", path=path, line=line, nonnegative=True
    )
    checks = (
        ("query_coverage", query_fraction, expected_query),
        ("reference_coverage", reference_fraction, expected_reference),
        ("reciprocal_coverage", reciprocal, expected_reciprocal),
        ("score", score, expected_score),
    )
    for column, observed, expected in checks:
        if not _decimal_close(observed, expected, tolerance):
            raise ChromosomeAssignmentError(
                f"{path.name}:{line}: {column} arithmetic mismatch; "
                f"observed={observed:.12g}, expected={expected:.12g}"
            )


def _validate_jcvi_row(
    row: Mapping[str, str], *, path: Path, line: int, policy: HomologyPolicy
) -> None:
    query_anchored = _integer(
        row["query_anchored_genes"], column="query_anchored_genes", path=path, line=line
    )
    query_eligible = _integer(
        row["query_eligible_genes"],
        column="query_eligible_genes",
        path=path,
        line=line,
        positive=True,
    )
    reference_anchored = _integer(
        row["reference_anchored_genes"],
        column="reference_anchored_genes",
        path=path,
        line=line,
    )
    reference_eligible = _integer(
        row["reference_eligible_genes"],
        column="reference_eligible_genes",
        path=path,
        line=line,
        positive=True,
    )
    anchor_pairs = _integer(
        row["unique_anchor_pairs"], column="unique_anchor_pairs", path=path, line=line
    )
    if query_anchored > query_eligible or reference_anchored > reference_eligible:
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: anchored genes cannot exceed eligible genes"
        )
    if anchor_pairs < max(query_anchored, reference_anchored):
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: unique_anchor_pairs cannot be smaller than unique anchored genes"
        )
    if anchor_pairs > query_anchored * reference_anchored:
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: unique_anchor_pairs exceeds the Cartesian gene-pair bound"
        )
    zero_states = (query_anchored == 0, reference_anchored == 0, anchor_pairs == 0)
    if len(set(zero_states)) != 1:
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: both anchored-gene counts and anchor pairs must "
            "be zero or non-zero together"
        )
    query_fraction = _finite_decimal(
        row["query_gene_coverage"],
        column="query_gene_coverage",
        path=path,
        line=line,
        nonnegative=True,
    )
    reference_fraction = _finite_decimal(
        row["reference_gene_coverage"],
        column="reference_gene_coverage",
        path=path,
        line=line,
        nonnegative=True,
    )
    if query_fraction > 1 or reference_fraction > 1:
        raise ChromosomeAssignmentError(
            f"{path.name}:{line}: JCVI gene-coverage fractions must lie in [0,1]"
        )
    tolerance = Decimal(str(policy.arithmetic_tolerance))
    score = _finite_decimal(
        row["score"], column="score", path=path, line=line, nonnegative=True
    )
    checks = (
        (
            "query_gene_coverage",
            query_fraction,
            Decimal(query_anchored) / Decimal(query_eligible),
        ),
        (
            "reference_gene_coverage",
            reference_fraction,
            Decimal(reference_anchored) / Decimal(reference_eligible),
        ),
        ("score", score, _decimal_harmonic(query_fraction, reference_fraction)),
    )
    for column, observed, expected in checks:
        if not _decimal_close(observed, expected, tolerance):
            raise ChromosomeAssignmentError(
                f"{path.name}:{line}: {column} arithmetic mismatch; "
                f"observed={observed:.12g}, expected={expected:.12g}"
            )


def read_score_matrix(
    path: str | Path | FileSnapshot, *, role: str, kind: str, policy: HomologyPolicy
) -> ScoreMatrix:
    """Read one exact long-form 29 x 29 score matrix."""

    snapshot = path if isinstance(path, FileSnapshot) else capture_snapshot(path)
    source = snapshot.path
    if kind not in {"nucleotide", "jcvi"}:
        raise ChromosomeAssignmentError(f"Unsupported matrix kind: {kind}")
    columns = NUCLEOTIDE_COLUMNS if kind == "nucleotide" else JCVI_COLUMNS
    raw_rows = _read_exact_rows(snapshot, columns)
    expected_rows = policy.expected_query_chromosomes * policy.expected_reference_chromosomes
    if len(raw_rows) != expected_rows:
        raise ChromosomeAssignmentError(
            f"{source.name}: exact {policy.expected_query_chromosomes} x "
            f"{policy.expected_reference_chromosomes} matrix requires {expected_rows} rows; "
            f"found {len(raw_rows)}"
        )

    cells: dict[tuple[str, str], MatrixCell] = {}
    canonical_sets: dict[str, set[str]] = {}
    normalized_rows: list[Mapping[str, str]] = []
    query_lengths: dict[str, int] = {}
    reference_lengths: dict[str, int] = {}
    query_eligible: dict[str, int] = {}
    reference_eligible: dict[str, int] = {}
    for line, row in raw_rows:
        query, reference, canonical, score = _validate_common_cell(
            row, path=source, line=line, expected_count=policy.expected_reference_chromosomes
        )
        key = (query, reference)
        if key in cells:
            raise ChromosomeAssignmentError(
                f"{source.name}:{line}: duplicate matrix cell {query!r}, {reference!r}"
            )
        if kind == "nucleotide":
            _validate_nucleotide_row(row, path=source, line=line, policy=policy)
            q_length = _integer(
                row["query_length_bp"],
                column="query_length_bp",
                path=source,
                line=line,
                positive=True,
            )
            r_length = _integer(
                row["reference_length_bp"],
                column="reference_length_bp",
                path=source,
                line=line,
                positive=True,
            )
            if query in query_lengths and query_lengths[query] != q_length:
                raise ChromosomeAssignmentError(
                    f"{source.name}:{line}: query_length_bp changes within {query!r}"
                )
            if reference in reference_lengths and reference_lengths[reference] != r_length:
                raise ChromosomeAssignmentError(
                    f"{source.name}:{line}: reference_length_bp changes within {reference!r}"
                )
            query_lengths[query] = q_length
            reference_lengths[reference] = r_length
        else:
            _validate_jcvi_row(row, path=source, line=line, policy=policy)
            q_eligible = _integer(
                row["query_eligible_genes"],
                column="query_eligible_genes",
                path=source,
                line=line,
                positive=True,
            )
            r_eligible = _integer(
                row["reference_eligible_genes"],
                column="reference_eligible_genes",
                path=source,
                line=line,
                positive=True,
            )
            if query in query_eligible and query_eligible[query] != q_eligible:
                raise ChromosomeAssignmentError(
                    f"{source.name}:{line}: query_eligible_genes changes within {query!r}"
                )
            if reference in reference_eligible and reference_eligible[reference] != r_eligible:
                raise ChromosomeAssignmentError(
                    f"{source.name}:{line}: reference_eligible_genes changes within {reference!r}"
                )
            query_eligible[query] = q_eligible
            reference_eligible[reference] = r_eligible
        canonical_sets.setdefault(reference, set()).add(canonical)
        cells[key] = MatrixCell(query, reference, canonical, score, dict(row))
        normalized_rows.append(dict(row))

    queries = tuple(sorted({key[0] for key in cells}, key=_total_key))
    references = tuple(sorted({key[1] for key in cells}, key=_total_key))
    if len(queries) != policy.expected_query_chromosomes:
        raise ChromosomeAssignmentError(
            f"{source.name}: expected exactly "
            f"{policy.expected_query_chromosomes} unique query IDs; "
            f"found {len(queries)}"
        )
    if len(references) != policy.expected_reference_chromosomes:
        raise ChromosomeAssignmentError(
            f"{source.name}: expected exactly {policy.expected_reference_chromosomes} unique "
            f"reference IDs; found {len(references)}"
        )
    expected_pairs = {(query, reference) for query in queries for reference in references}
    if set(cells) != expected_pairs:
        missing = sorted(
            expected_pairs.difference(cells),
            key=lambda pair: (_total_key(pair[0]), _total_key(pair[1])),
        )
        raise ChromosomeAssignmentError(
            f"{source.name}: matrix is not a complete Cartesian product; "
            f"missing={','.join(f'{q}:{r}' for q, r in missing[:5])}"
        )
    inconsistent = [
        reference for reference, labels in canonical_sets.items() if len(labels) != 1
    ]
    if inconsistent:
        raise ChromosomeAssignmentError(
            f"{source.name}: reference-to-canonical labels change across rows: "
            + ",".join(sorted(inconsistent, key=_total_key))
        )
    canonical_map = {
        reference: next(iter(canonical_sets[reference])) for reference in references
    }
    if set(canonical_map.values()) != set(_canonical_labels(policy.expected_reference_chromosomes)):
        raise ChromosomeAssignmentError(
            f"{source.name}: reference-to-canonical map must be a bijection onto Chr01..Chr29"
        )
    normalized_rows.sort(
        key=lambda row: (
            _total_key(row["query_chromosome"]),
            _total_key(row["reference_chromosome"]),
        )
    )
    return ScoreMatrix(
        role=role,
        kind=kind,
        source=source,
        columns=tuple(columns),
        queries=queries,
        references=references,
        canonical_by_reference=canonical_map,
        cells=cells,
        normalized_rows=tuple(normalized_rows),
    )


def _hungarian_maximize(scores: Sequence[Sequence[float]]) -> tuple[int, ...]:
    """Return a deterministic maximum-weight square assignment.

    The result is indexed by row and contains the assigned column index.  The
    implementation is the O(n^3) Hungarian primal-dual algorithm; iteration in
    input order gives a stable answer when equal costs remain.  Biological
    acceptance still requires unique row and column best edges, so a stable
    tie-break can never turn ambiguous evidence into a passing map.
    """

    count = len(scores)
    if count == 0 or any(len(row) != count for row in scores):
        raise ChromosomeAssignmentError("Hungarian assignment requires a non-empty square matrix")
    if any(not math.isfinite(value) for row in scores for value in row):
        raise ChromosomeAssignmentError("Hungarian assignment scores must be finite")
    maximum = max(value for row in scores for value in row)
    costs = [[maximum - value for value in row] for row in scores]
    u = [0.0] * (count + 1)
    v = [0.0] * (count + 1)
    p = [0] * (count + 1)
    way = [0] * (count + 1)
    for row_index in range(1, count + 1):
        p[0] = row_index
        column0 = 0
        minimum = [math.inf] * (count + 1)
        used = [False] * (count + 1)
        while True:
            used[column0] = True
            active_row = p[column0]
            delta = math.inf
            column1 = 0
            for column in range(1, count + 1):
                if used[column]:
                    continue
                reduced = costs[active_row - 1][column - 1] - u[active_row] - v[column]
                if reduced < minimum[column]:
                    minimum[column] = reduced
                    way[column] = column0
                if minimum[column] < delta:
                    delta = minimum[column]
                    column1 = column
            if not math.isfinite(delta):
                raise ChromosomeAssignmentError(
                    "Hungarian assignment could not find an augmenting path"
                )
            for column in range(count + 1):
                if used[column]:
                    u[p[column]] += delta
                    v[column] -= delta
                else:
                    minimum[column] -= delta
            column0 = column1
            if p[column0] == 0:
                break
        while True:
            column1 = way[column0]
            p[column0] = p[column1]
            column0 = column1
            if column0 == 0:
                break
    assignment = [-1] * count
    for column in range(1, count + 1):
        if p[column] <= 0:
            raise ChromosomeAssignmentError("Hungarian assignment is incomplete")
        assignment[p[column] - 1] = column - 1
    if sorted(assignment) != list(range(count)):
        raise ChromosomeAssignmentError("Hungarian assignment is not one-to-one")
    return tuple(assignment)


def _top_metrics(
    identifiers: Sequence[str], values: Sequence[float], tolerance: float
) -> tuple[str, float, float, float | None, float, bool]:
    ordered = sorted(
        zip(identifiers, values), key=lambda item: (-item[1], _total_key(item[0]))
    )
    best_identifier, best = ordered[0]
    second = ordered[1][1]
    unique = not _close(best, second, tolerance)
    ratio = (
        None
        if second == 0.0 and best > 0.0
        else (best / second if second > 0.0 else 0.0)
    )
    margin = (best - second) / best if best > 0.0 else 0.0
    return best_identifier, best, second, ratio, margin, unique


def solve_score_matrix(
    matrix: ScoreMatrix, policy: HomologyPolicy
) -> dict[str, AssignmentEvidence]:
    """Solve one matrix globally and calculate every local acceptance gate."""

    scores = [
        [matrix.cells[(query, reference)].score for reference in matrix.references]
        for query in matrix.queries
    ]
    assigned_columns = _hungarian_maximize(scores)
    evidence: dict[str, AssignmentEvidence] = {}
    for row_index, query in enumerate(matrix.queries):
        reference = matrix.references[assigned_columns[row_index]]
        cell = matrix.cells[(query, reference)]
        row_values = [matrix.cells[(query, item)].score for item in matrix.references]
        (
            row_best_reference,
            row_best,
            row_second,
            row_ratio,
            row_margin,
            row_unique,
        ) = _top_metrics(matrix.references, row_values, policy.arithmetic_tolerance)
        column_values = [matrix.cells[(item, reference)].score for item in matrix.queries]
        (
            column_best_query,
            column_best,
            column_second,
            column_ratio,
            column_margin,
            column_unique,
        ) = _top_metrics(matrix.queries, column_values, policy.arithmetic_tolerance)
        reciprocal_best = (
            row_unique
            and column_unique
            and row_best_reference == reference
            and column_best_query == query
        )
        row_ratio_pass = (
            row_ratio is None
            or row_ratio + policy.arithmetic_tolerance
            >= policy.minimum_top_second_ratio
        )
        column_ratio_pass = (
            column_ratio is None
            or column_ratio + policy.arithmetic_tolerance
            >= policy.minimum_top_second_ratio
        )
        separation = (
            row_ratio_pass
            and column_ratio_pass
            and row_margin + policy.arithmetic_tolerance
            >= policy.minimum_normalized_score_margin
            and column_margin + policy.arithmetic_tolerance
            >= policy.minimum_normalized_score_margin
        )
        reciprocal_coverage: float | None = None
        reciprocal_gene_coverage: float | None = None
        unique_anchor_pairs: int | None = None
        orientation = ""
        coverage_gate = True
        absolute_support_gate = True
        if matrix.kind == "nucleotide":
            reciprocal_coverage = float(cell.raw["reciprocal_coverage"])
            orientation = cell.raw["orientation"]
            coverage_gate = (
                reciprocal_coverage + policy.arithmetic_tolerance
                >= policy.minimum_assigned_reciprocal_nucleotide_coverage
            )
        else:
            reciprocal_gene_coverage = min(
                float(cell.raw["query_gene_coverage"]),
                float(cell.raw["reference_gene_coverage"]),
            )
            unique_anchor_pairs = _integer(
                cell.raw["unique_anchor_pairs"],
                column="unique_anchor_pairs",
                path=matrix.source,
                line=0,
            )
            coverage_gate = (
                reciprocal_gene_coverage + policy.arithmetic_tolerance
                >= policy.minimum_reciprocal_gene_coverage
            )
            absolute_support_gate = (
                coverage_gate
                and unique_anchor_pairs >= policy.minimum_unique_anchor_pairs
                and cell.score + policy.arithmetic_tolerance
                >= policy.minimum_assigned_jcvi_score
            )
        failures = []
        if not reciprocal_best:
            failures.append("AMBIGUOUS_RECIPROCAL_BEST")
        if not separation:
            failures.append("AMBIGUOUS_SEPARATION")
        if not coverage_gate:
            failures.append(
                "LOW_RECIPROCAL_COVERAGE"
                if matrix.kind == "nucleotide"
                else "LOW_JCVI_ABSOLUTE_SUPPORT"
            )
        elif not absolute_support_gate:
            failures.append("LOW_JCVI_ABSOLUTE_SUPPORT")
        evidence[query] = AssignmentEvidence(
            query=query,
            reference=reference,
            canonical=cell.canonical,
            score=cell.score,
            reciprocal_coverage=reciprocal_coverage,
            reciprocal_gene_coverage=reciprocal_gene_coverage,
            unique_anchor_pairs=unique_anchor_pairs,
            orientation=orientation,
            row_best_reference=row_best_reference,
            row_best_score=row_best,
            row_second_score=row_second,
            row_ratio=row_ratio,
            row_margin=row_margin,
            row_unique_best=row_unique,
            column_best_query=column_best_query,
            column_best_score=column_best,
            column_second_score=column_second,
            column_ratio=column_ratio,
            column_margin=column_margin,
            column_unique_best=column_unique,
            reciprocal_best=reciprocal_best,
            separation_gate=separation,
            coverage_gate=coverage_gate,
            absolute_support_gate=absolute_support_gate,
            matrix_gate=(
                reciprocal_best
                and separation
                and coverage_gate
                and absolute_support_gate
            ),
            failure_reasons=tuple(failures),
        )
    if len({item.reference for item in evidence.values()}) != len(matrix.references):
        raise ChromosomeAssignmentError(f"{matrix.role}: global assignment is not bijective")
    return evidence


def _write_tsv(path: Path, columns: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _evidence_rows(evidence: Mapping[str, AssignmentEvidence]) -> list[dict[str, str]]:
    rows = []
    for query in sorted(evidence, key=_total_key):
        item = evidence[query]
        rows.append(
            {
                "query_chromosome": query,
                "assigned_reference_chromosome": item.reference,
                "assigned_canonical_chromosome": item.canonical,
                "assigned_score": _format_float(item.score),
                "assigned_reciprocal_coverage": (
                    ""
                    if item.reciprocal_coverage is None
                    else _format_float(item.reciprocal_coverage)
                ),
                "assigned_reciprocal_gene_coverage": (
                    ""
                    if item.reciprocal_gene_coverage is None
                    else _format_float(item.reciprocal_gene_coverage)
                ),
                "assigned_unique_anchor_pairs": (
                    "" if item.unique_anchor_pairs is None else str(item.unique_anchor_pairs)
                ),
                "assigned_orientation": item.orientation,
                "row_best_reference_chromosome": item.row_best_reference,
                "row_best_score": _format_float(item.row_best_score),
                "row_second_score": _format_float(item.row_second_score),
                "row_top_second_ratio": _format_float(item.row_ratio),
                "row_normalized_margin": _format_float(item.row_margin),
                "row_unique_best": _format_bool(item.row_unique_best),
                "column_best_query_chromosome": item.column_best_query,
                "column_best_score": _format_float(item.column_best_score),
                "column_second_score": _format_float(item.column_second_score),
                "column_top_second_ratio": _format_float(item.column_ratio),
                "column_normalized_margin": _format_float(item.column_margin),
                "column_unique_best": _format_bool(item.column_unique_best),
                "row_and_column_reciprocal_best": _format_bool(item.reciprocal_best),
                "separation_gate": _format_bool(item.separation_gate),
                "coverage_gate": _format_bool(item.coverage_gate),
                "absolute_support_gate": _format_bool(item.absolute_support_gate),
                "matrix_gate": _format_bool(item.matrix_gate),
                "failure_reasons": ";".join(item.failure_reasons),
            }
        )
    return rows


def _matrix_pair_compatible(first: ScoreMatrix, second: ScoreMatrix, reference: str) -> None:
    if first.queries != second.queries:
        raise ChromosomeAssignmentError(
            f"{reference}: nucleotide and JCVI query ID sets/order do not agree exactly"
        )
    if first.references != second.references:
        raise ChromosomeAssignmentError(
            f"{reference}: nucleotide and JCVI reference ID sets/order do not agree exactly"
        )
    if dict(first.canonical_by_reference) != dict(second.canonical_by_reference):
        raise ChromosomeAssignmentError(
            f"{reference}: nucleotide and JCVI reference-to-canonical maps differ"
        )


def _matrix_query_denominators(matrix: ScoreMatrix, column: str) -> dict[str, int]:
    return {
        query: _integer(
            matrix.cells[(query, matrix.references[0])].raw[column],
            column=column,
            path=matrix.source,
            line=0,
            positive=True,
        )
        for query in matrix.queries
    }


def _ordered_failure_states(states: set[str]) -> tuple[str, ...]:
    ordered = [state for state in FAILURE_PRIORITY if state in states]
    ordered.extend(sorted(states.difference(FAILURE_PRIORITY)))
    return tuple(ordered)


def _generation_parameters(policy: HomologyPolicy, kind: str) -> dict[str, Any]:
    if kind == "nucleotide":
        return {
            "minimap2_version": policy.minimap2_version,
            "minimap2_command_template": list(policy.minimap2_command_template),
            "required_paf_tags": list(REQUIRED_PAF_TAGS),
            "minimap2_primary_only": policy.minimap2_primary_only,
            "minimum_mapq": policy.minimum_mapq,
            "minimum_alignment_block_bp": policy.minimum_alignment_block_bp,
            "maximum_de": policy.maximum_de,
            "coverage_arithmetic": policy.coverage_arithmetic,
            "nucleotide_score_formula": policy.nucleotide_score_formula,
            "reciprocal_nucleotide_coverage_formula": (
                policy.reciprocal_nucleotide_coverage_formula
            ),
            "arithmetic_tolerance": policy.arithmetic_tolerance,
        }
    if kind == "jcvi":
        return {
            "aligner": policy.jcvi_aligner,
            "database_type": policy.jcvi_database_type,
            "cscore": policy.jcvi_cscore,
            "tandem_nmax": policy.jcvi_tandem_nmax,
            "maximum_gene_distance": policy.jcvi_maximum_gene_distance,
            "minimum_anchor_block_size": policy.jcvi_minimum_anchor_block_size,
            "coverage_anchor_source": policy.jcvi_coverage_anchor_source,
            "jcvi_anchor_counting": policy.jcvi_anchor_counting,
            "jcvi_score_formula": policy.jcvi_score_formula,
            "arithmetic_tolerance": policy.arithmetic_tolerance,
        }
    raise ChromosomeAssignmentError(f"Unsupported matrix kind: {kind}")


def _acquire_output_lock(output: Path) -> tuple[int, Path, tuple[int, int]]:
    output.parent.mkdir(parents=True, exist_ok=True)
    lock_path = output.parent / f".{output.name}.chromosome_assignment.lock"
    try:
        descriptor = os.open(
            lock_path,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0),
            0o600,
        )
    except FileExistsError as error:
        raise ChromosomeAssignmentError(
            f"Another assignment owns the output lock: {lock_path.name}"
        ) from error
    state = os.fstat(descriptor)
    os.write(descriptor, f"pid={os.getpid()}\n".encode("ascii"))
    return descriptor, lock_path, (state.st_dev, state.st_ino)


def _release_output_lock(
    descriptor: int, lock_path: Path, expected_inode: tuple[int, int]
) -> None:
    os.close(descriptor)
    try:
        state = lock_path.stat(follow_symlinks=False)
    except FileNotFoundError:
        return
    if (state.st_dev, state.st_ino) == expected_inode:
        lock_path.unlink()


def _rename_no_replace(source: Path, destination: Path) -> None:
    """Atomically rename a directory without replacing an existing destination."""

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
    if result is not None:
        if result == 0:
            return
        error_number = ctypes.get_errno()
        if error_number in {errno.EEXIST, errno.ENOTEMPTY}:
            raise ChromosomeAssignmentError(
                f"Output appeared during publication; refusing overwrite: {destination.name}"
            )
        raise ChromosomeAssignmentError(
            f"Atomic no-replace rename failed for {destination.name}: "
            f"{os.strerror(error_number)}"
        )
    raise ChromosomeAssignmentError(
        "This platform lacks an atomic no-replace directory rename primitive; "
        "refusing non-atomic publication"
    )


def assign_chromosome_homology(
    *,
    nucleotide_hy4a: str | Path,
    jcvi_hy4a: str | Path,
    nucleotide_hy4p: str | Path,
    jcvi_hy4p: str | Path,
    nucleotide_hy4a_provenance: str | Path,
    jcvi_hy4a_provenance: str | Path,
    nucleotide_hy4p_provenance: str | Path,
    jcvi_hy4p_provenance: str | Path,
    parameters: str | Path,
    target_asset_registry: str | Path,
    reference_asset_registry: str | Path,
    reference_chromosome_map_registry: str | Path,
    assembly_unit_id: str,
    target_scope_id: str,
    trusted_repository_commit: str,
    output_dir: str | Path,
) -> ChromosomeAssignmentResult:
    """Validate, solve, gate, and atomically publish one assignment bundle."""

    if not SAFE_UNIT_ID.fullmatch(assembly_unit_id):
        raise ChromosomeAssignmentError(
            "assembly_unit_id must be path-safe and contain only letters, numbers, '.', '_', or '-'"
        )
    if not SAFE_UNIT_ID.fullmatch(target_scope_id):
        raise ChromosomeAssignmentError(
            "target_scope_id must be path-safe and contain only letters, numbers, '.', '_', or '-'"
        )
    if not GIT_COMMIT_RE.fullmatch(trusted_repository_commit):
        raise ChromosomeAssignmentError(
            "trusted_repository_commit must be one full lowercase 40- or 64-hex Git commit ID"
        )
    output = Path(output_dir).expanduser()
    if os.path.lexists(output):
        raise ChromosomeAssignmentError(
            f"Output directory already exists; refusing overwrite: {output.name}"
        )
    matrix_paths = {
        "nucleotide_hy4a": nucleotide_hy4a,
        "jcvi_hy4a": jcvi_hy4a,
        "nucleotide_hy4p": nucleotide_hy4p,
        "jcvi_hy4p": jcvi_hy4p,
    }
    sidecar_paths = {
        "nucleotide_hy4a": nucleotide_hy4a_provenance,
        "jcvi_hy4a": jcvi_hy4a_provenance,
        "nucleotide_hy4p": nucleotide_hy4p_provenance,
        "jcvi_hy4p": jcvi_hy4p_provenance,
    }
    matrix_snapshots = {
        role: capture_snapshot(path) for role, path in matrix_paths.items()
    }
    sidecar_snapshots = {
        role: capture_snapshot(path) for role, path in sidecar_paths.items()
    }
    reject_duplicate_snapshots(matrix_snapshots, label="score matrices")
    reject_duplicate_snapshots(sidecar_snapshots, label="matrix provenance sidecars")

    parameter_snapshot = capture_snapshot(parameters)
    target_registry_snapshot = capture_snapshot(target_asset_registry)
    reference_registry_snapshot = capture_snapshot(reference_asset_registry)
    reference_map_snapshot = capture_snapshot(reference_chromosome_map_registry)
    policy = load_homology_policy(parameter_snapshot)
    if reference_registry_snapshot.sha256 != policy.reference_asset_registry_sha256:
        raise ChromosomeAssignmentError(
            "Reference asset registry checksum conflicts with the frozen policy"
        )
    if reference_map_snapshot.sha256 != policy.reference_chromosome_map_registry_sha256:
        raise ChromosomeAssignmentError(
            "Reference chromosome-map registry checksum conflicts with the frozen policy"
        )
    target_registry = read_target_registry(target_registry_snapshot)
    reference_registries = read_reference_registries(
        reference_registry_snapshot,
        reference_map_snapshot,
        coordinate_reference=policy.coordinate_reference,
        confirmation_reference=policy.confirmation_reference,
        expected_chromosomes=policy.expected_reference_chromosomes,
    )

    role_specs = {
        "nucleotide_hy4a": (
            "nucleotide",
            "HY4A",
            policy.coordinate_reference,
            "primary_coordinate_reference",
        ),
        "jcvi_hy4a": (
            "jcvi",
            "HY4A",
            policy.coordinate_reference,
            "primary_coordinate_reference",
        ),
        "nucleotide_hy4p": (
            "nucleotide",
            "HY4P",
            policy.confirmation_reference,
            "independent_confirmation_reference",
        ),
        "jcvi_hy4p": (
            "jcvi",
            "HY4P",
            policy.confirmation_reference,
            "independent_confirmation_reference",
        ),
    }

    matrices: dict[str, ScoreMatrix] = {}
    matrix_provenance: dict[str, MatrixProvenance] = {}
    for role in sorted(role_specs, key=_total_key):
        kind, slot, reference_id, reference_role = role_specs[role]
        matrix = read_score_matrix(
            matrix_snapshots[role], role=role, kind=kind, policy=policy
        )
        frozen_map = reference_registries.chromosome_maps[reference_id]
        if tuple(sorted(frozen_map, key=_total_key)) != matrix.references:
            raise ChromosomeAssignmentError(
                f"{role}: matrix reference IDs do not equal the frozen {slot} map"
            )
        if dict(matrix.canonical_by_reference) != dict(frozen_map):
            raise ChromosomeAssignmentError(
                f"{role}: matrix canonical labels conflict with the frozen {slot} map"
            )
        matrices[role] = matrix
        matrix_provenance[role] = validate_matrix_provenance(
            sidecar=sidecar_snapshots[role],
            matrix=matrix_snapshots[role],
            assembly_unit_id=assembly_unit_id,
            target_scope_id=target_scope_id,
            matrix_role=role,
            matrix_kind=kind,
            reference_slot=slot,
            reference_id=reference_id,
            reference_role=reference_role,
            reference_map_id=reference_registries.map_ids[reference_id],
            expected_generation_parameters=_generation_parameters(policy, kind),
            target_registry=target_registry,
            reference_registries=reference_registries,
        )

    report_snapshots = {
        role: item.upstream_report_snapshot
        for role, item in matrix_provenance.items()
    }
    role_bound_snapshots: dict[str, FileSnapshot] = {}
    for role in sorted(role_specs, key=_total_key):
        role_bound_snapshots[f"matrix:{role}"] = matrix_snapshots[role]
        role_bound_snapshots[f"sidecar:{role}"] = sidecar_snapshots[role]
        role_bound_snapshots[f"upstream_report:{role}"] = report_snapshots[role]
    reject_duplicate_snapshots(
        role_bound_snapshots,
        label="matrix, sidecar, and upstream-report role bindings",
    )
    _matrix_pair_compatible(matrices["nucleotide_hy4a"], matrices["jcvi_hy4a"], "HY4A")
    _matrix_pair_compatible(matrices["nucleotide_hy4p"], matrices["jcvi_hy4p"], "HY4P")
    query_sets = {matrix.queries for matrix in matrices.values()}
    if len(query_sets) != 1:
        raise ChromosomeAssignmentError(
            "All four matrices must contain the exact same query IDs"
        )
    if _matrix_query_denominators(matrices["nucleotide_hy4a"], "query_length_bp") != (
        _matrix_query_denominators(matrices["nucleotide_hy4p"], "query_length_bp")
    ):
        raise ChromosomeAssignmentError(
            "HY4A/HY4P nucleotide matrices disagree on query lengths"
        )
    if _matrix_query_denominators(matrices["jcvi_hy4a"], "query_eligible_genes") != (
        _matrix_query_denominators(matrices["jcvi_hy4p"], "query_eligible_genes")
    ):
        raise ChromosomeAssignmentError(
            "HY4A/HY4P JCVI matrices disagree on query gene denominators"
        )

    assignments = {
        role: solve_score_matrix(matrix, policy) for role, matrix in matrices.items()
    }
    queries = matrices["nucleotide_hy4a"].queries
    combined_rows: list[dict[str, str]] = []
    failure_states: set[str] = set()
    agreement_a_count = 0
    agreement_p_count = 0
    cross_reference_count = 0
    all_matrix_count = 0

    for query in queries:
        na = assignments["nucleotide_hy4a"][query]
        ja = assignments["jcvi_hy4a"][query]
        np = assignments["nucleotide_hy4p"][query]
        jp = assignments["jcvi_hy4p"][query]
        agreement_a = na.reference == ja.reference and na.canonical == ja.canonical
        agreement_p = np.reference == jp.reference and np.canonical == jp.canonical
        cross_reference = agreement_a and agreement_p and na.canonical == np.canonical
        all_matrix = all(item.matrix_gate for item in (na, ja, np, jp))
        if agreement_a:
            agreement_a_count += 1
        else:
            failure_states.add("CONFLICT_NUCLEOTIDE_JCVI")
        if agreement_p:
            agreement_p_count += 1
        else:
            failure_states.add("CONFLICT_NUCLEOTIDE_JCVI")
        if cross_reference:
            cross_reference_count += 1
        elif agreement_a and agreement_p:
            failure_states.add("HY4A_HY4P_DISAGREEMENT")
        if all_matrix:
            all_matrix_count += 1
        for item in (na, ja, np, jp):
            failure_states.update(item.failure_reasons)
        row_failures = set()
        if not agreement_a or not agreement_p:
            row_failures.add("CONFLICT_NUCLEOTIDE_JCVI")
        if agreement_a and agreement_p and not cross_reference:
            row_failures.add("HY4A_HY4P_DISAGREEMENT")
        for item in (na, ja, np, jp):
            row_failures.update(item.failure_reasons)
        candidate = na.canonical if cross_reference else ""
        combined_rows.append(
            {
                "query_chromosome": query,
                "nucleotide_hy4a_reference": na.reference,
                "nucleotide_hy4a_canonical": na.canonical,
                "jcvi_hy4a_reference": ja.reference,
                "jcvi_hy4a_canonical": ja.canonical,
                "hy4a_nucleotide_jcvi_agreement": _format_bool(agreement_a),
                "nucleotide_hy4p_reference": np.reference,
                "nucleotide_hy4p_canonical": np.canonical,
                "jcvi_hy4p_reference": jp.reference,
                "jcvi_hy4p_canonical": jp.canonical,
                "hy4p_nucleotide_jcvi_agreement": _format_bool(agreement_p),
                "hy4a_hy4p_label_agreement": _format_bool(cross_reference),
                "candidate_canonical_chromosome": candidate,
                "all_matrix_gates": _format_bool(all_matrix),
                "final_chromosome": "",  # filled only after the whole-unit gate passes
                "status": (
                    "PASS_AUTO"
                    if not row_failures
                    else _ordered_failure_states(row_failures)[0]
                ),
                "failure_reasons": ";".join(_ordered_failure_states(row_failures)),
            }
        )

    canonical_candidates = [row["candidate_canonical_chromosome"] for row in combined_rows]
    if all(canonical_candidates) and len(set(canonical_candidates)) != len(queries):
        failure_states.add("NON_BIJECTIVE")
    ordered_failures = _ordered_failure_states(failure_states)
    status = "PASS_AUTO" if not ordered_failures else ordered_failures[0]
    publication_gate = "PASS" if status == "PASS_AUTO" else "FAIL"
    if status == "PASS_AUTO":
        for row in combined_rows:
            row["final_chromosome"] = row["candidate_canonical_chromosome"]
            row["status"] = "PASS_AUTO"
            row["failure_reasons"] = ""
    else:
        for row in combined_rows:
            if row["status"] == "PASS_AUTO":
                row["status"] = "NOT_PUBLISHED_UNIT_FAILURE"
                row["failure_reasons"] = f"UNIT_GATE_FAILED:{status}"
    final_rows = [
        {
            "query_chromosome": row["query_chromosome"],
            "final_chromosome": row["final_chromosome"],
            "coordinate_reference": policy.coordinate_reference,
            "confirmation_reference": policy.confirmation_reference,
            "status": "PASS_AUTO",
        }
        for row in combined_rows
        if status == "PASS_AUTO"
    ]

    input_provenance = {
        "parameters": parameter_snapshot.public_binding(),
        "target_asset_registry": target_registry_snapshot.public_binding(),
        "reference_asset_registry": reference_registry_snapshot.public_binding(),
        "reference_chromosome_map_registry": reference_map_snapshot.public_binding(),
        "matrix_roles": {
            role: {
                "matrix": matrix_snapshots[role].public_binding(),
                "provenance_sidecar": sidecar_snapshots[role].public_binding(),
                "upstream_validation_report": report_snapshots[role].public_binding(),
                "matrix_kind": matrix_provenance[role].matrix_kind,
                "reference_slot": matrix_provenance[role].reference_slot,
                "reference_id": matrix_provenance[role].reference_id,
                "reference_role": matrix_provenance[role].reference_role,
                "reference_map_id": matrix_provenance[role].reference_map_id,
            }
            for role in sorted(role_specs, key=_total_key)
        },
    }

    lock_descriptor, lock_path, lock_inode = _acquire_output_lock(output)
    try:
        if os.path.lexists(output):
            raise ChromosomeAssignmentError(
                f"Output appeared before staging; refusing overwrite: {output.name}"
            )
        staging = Path(
            tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent)
        )
    except Exception:
        _release_output_lock(lock_descriptor, lock_path, lock_inode)
        raise
    try:
        output_files: list[str] = []
        for role, matrix in matrices.items():
            matrix_name = f"{role}.normalized_score_matrix.tsv"
            evidence_name = f"{role}.assignment_evidence.tsv"
            _write_tsv(staging / matrix_name, matrix.columns, matrix.normalized_rows)
            _write_tsv(staging / evidence_name, EVIDENCE_COLUMNS, _evidence_rows(assignments[role]))
            output_files.extend((matrix_name, evidence_name))
        assignment_name = f"{assembly_unit_id}.chromosome_assignment.tsv"
        _write_tsv(staging / assignment_name, COMBINED_COLUMNS, combined_rows)
        output_files.append(assignment_name)
        if status == "PASS_AUTO":
            final_name = f"{assembly_unit_id}.final_chromosome_map.tsv"
            _write_tsv(staging / final_name, FINAL_MAP_COLUMNS, final_rows)
            output_files.append(final_name)

        summary_name = f"{assembly_unit_id}.chromosome_assignment.summary.tsv"
        summary_row = {
            "assembly_unit_id": assembly_unit_id,
            "target_scope_id": target_scope_id,
            "trusted_repository_commit": trusted_repository_commit,
            "status": status,
            "publication_gate": publication_gate,
            "failure_states": ";".join(ordered_failures),
            "expected_query_chromosomes": str(policy.expected_query_chromosomes),
            "observed_query_chromosomes": str(len(queries)),
            "expected_reference_chromosomes_per_reference": str(
                policy.expected_reference_chromosomes
            ),
            "nucleotide_jcvi_agreement_count_hy4a": str(agreement_a_count),
            "nucleotide_jcvi_agreement_count_hy4p": str(agreement_p_count),
            "hy4a_hy4p_label_agreement_count": str(cross_reference_count),
            "rows_passing_all_matrix_gates": str(all_matrix_count),
            "final_map_row_count": str(len(final_rows)),
        }
        _write_tsv(staging / summary_name, SUMMARY_COLUMNS, [summary_row])
        output_files.append(summary_name)

        validation_name = "validation.json"
        validation = {
            "workflow": "chromosome_homology_assignment",
            "workflow_version": WORKFLOW_VERSION,
            "matrix_schema_version": MATRIX_SCHEMA_VERSION,
            "matrix_provenance_schema_version": PROVENANCE_SCHEMA_VERSION,
            "assembly_unit_id": assembly_unit_id,
            "target_scope_id": target_scope_id,
            "trusted_repository_commit": trusted_repository_commit,
            "status": status,
            "publication_gate": publication_gate,
            "failure_states": list(ordered_failures),
            "checks": {
                "exact_query_scope": len(queries) == policy.expected_query_chromosomes,
                "exact_reference_scope_each_matrix": all(
                    len(matrix.references) == policy.expected_reference_chromosomes
                    for matrix in matrices.values()
                ),
                "four_role_provenance_bindings_passed": len(matrix_provenance) == 4,
                "frozen_reference_maps_reconciled": all(
                    dict(matrices[role].canonical_by_reference)
                    == dict(
                        reference_registries.chromosome_maps[
                            role_specs[role][2]
                        ]
                    )
                    for role in role_specs
                ),
                "global_assignments_bijective": all(
                    len({item.reference for item in evidence.values()}) == len(queries)
                    for evidence in assignments.values()
                ),
                "nucleotide_jcvi_agreement_hy4a": agreement_a_count == len(queries),
                "nucleotide_jcvi_agreement_hy4p": agreement_p_count == len(queries),
                "hy4a_hy4p_label_agreement": cross_reference_count == len(queries),
                "all_local_matrix_gates": all_matrix_count == len(queries),
                "final_map_is_complete_bijection": (
                    len(final_rows) == len(queries)
                    and len({row["final_chromosome"] for row in final_rows}) == len(queries)
                ) if status == "PASS_AUTO" else False,
            },
            "counts": {
                "query_chromosomes": len(queries),
                "matrix_rows_each": len(queries) * policy.expected_reference_chromosomes,
                "nucleotide_jcvi_agreement_hy4a": agreement_a_count,
                "nucleotide_jcvi_agreement_hy4p": agreement_p_count,
                "hy4a_hy4p_label_agreement": cross_reference_count,
                "all_matrix_gates": all_matrix_count,
                "final_map_rows": len(final_rows),
            },
        }
        (staging / validation_name).write_text(
            json.dumps(
                validation,
                indent=2,
                sort_keys=True,
                ensure_ascii=True,
                allow_nan=False,
            )
            + "\n",
            encoding="utf-8",
        )
        output_files.append(validation_name)

        manifest_name = "run_manifest.json"
        manifest = {
            "workflow": "chromosome_homology_assignment",
            "workflow_version": WORKFLOW_VERSION,
            "matrix_schema_version": MATRIX_SCHEMA_VERSION,
            "matrix_provenance_schema_version": PROVENANCE_SCHEMA_VERSION,
            "assembly_unit_id": assembly_unit_id,
            "target_scope_id": target_scope_id,
            "trusted_repository_commit": trusted_repository_commit,
            "status": status,
            "publication_gate": publication_gate,
            "failure_states": list(ordered_failures),
            "inputs": input_provenance,
            "policy": {
                "coordinate_reference": policy.coordinate_reference,
                "independent_confirmation_reference": policy.confirmation_reference,
                "matrix_schema_version": policy.matrix_schema_version,
                "matrix_provenance_schema_version": (
                    policy.matrix_provenance_schema_version
                ),
                "reference_asset_registry_sha256": (
                    policy.reference_asset_registry_sha256
                ),
                "reference_chromosome_map_registry_sha256": (
                    policy.reference_chromosome_map_registry_sha256
                ),
                "expected_query_chromosomes": policy.expected_query_chromosomes,
                "expected_reference_chromosomes": policy.expected_reference_chromosomes,
                "assignment_method": policy.assignment_method,
                "minimap2_version": policy.minimap2_version,
                "minimap2_command_template": list(policy.minimap2_command_template),
                "required_paf_tags": list(REQUIRED_PAF_TAGS),
                "minimap2_primary_only": policy.minimap2_primary_only,
                "minimum_mapq": policy.minimum_mapq,
                "minimum_alignment_block_bp": policy.minimum_alignment_block_bp,
                "maximum_de": policy.maximum_de,
                "coverage_arithmetic": policy.coverage_arithmetic,
                "minimum_top_second_ratio": policy.minimum_top_second_ratio,
                "minimum_normalized_score_margin": policy.minimum_normalized_score_margin,
                "minimum_assigned_reciprocal_nucleotide_coverage": (
                    policy.minimum_assigned_reciprocal_nucleotide_coverage
                ),
                "minimum_unique_anchor_pairs": policy.minimum_unique_anchor_pairs,
                "minimum_reciprocal_gene_coverage": (
                    policy.minimum_reciprocal_gene_coverage
                ),
                "minimum_assigned_jcvi_score": policy.minimum_assigned_jcvi_score,
                "require_row_and_column_reciprocal_best": (
                    policy.require_row_and_column_reciprocal_best
                ),
                "require_nucleotide_and_jcvi_assignment_agreement": (
                    policy.require_nucleotide_and_jcvi_assignment_agreement
                ),
                "require_hy4a_hy4p_label_agreement": policy.require_hy4a_hy4p_label_agreement,
                "reverse_complement_allowed": policy.reverse_complement_allowed,
                "nucleotide_score_formula": policy.nucleotide_score_formula,
                "reciprocal_nucleotide_coverage_formula": (
                    policy.reciprocal_nucleotide_coverage_formula
                ),
                "jcvi_score_formula": policy.jcvi_score_formula,
                "jcvi_anchor_counting": policy.jcvi_anchor_counting,
                "jcvi_aligner": policy.jcvi_aligner,
                "jcvi_database_type": policy.jcvi_database_type,
                "jcvi_cscore": policy.jcvi_cscore,
                "jcvi_tandem_nmax": policy.jcvi_tandem_nmax,
                "jcvi_maximum_gene_distance": policy.jcvi_maximum_gene_distance,
                "jcvi_minimum_anchor_block_size": (
                    policy.jcvi_minimum_anchor_block_size
                ),
                "jcvi_coverage_anchor_source": policy.jcvi_coverage_anchor_source,
                "arithmetic_tolerance": policy.arithmetic_tolerance,
                "failure_policy": policy.failure_policy,
                "reverse_complement_action": "none; orientation is audit-only",
            },
            "outputs": sorted([*output_files, manifest_name], key=_total_key),
        }
        (staging / manifest_name).write_text(
            json.dumps(
                manifest,
                indent=2,
                sort_keys=True,
                ensure_ascii=True,
                allow_nan=False,
            )
            + "\n",
            encoding="utf-8",
        )
        output_files.append(manifest_name)

        checksum_rows = []
        for name in sorted(output_files, key=_total_key):
            path = staging / name
            checksum_rows.append(
                {"file": name, "bytes": str(path.stat().st_size), "sha256": _sha256(path)}
            )
            path.chmod(PUBLIC_FILE_MODE)
        _write_tsv(staging / "checksums.tsv", ("file", "bytes", "sha256"), checksum_rows)
        (staging / "checksums.tsv").chmod(PUBLIC_FILE_MODE)
        staging.chmod(PUBLIC_DIRECTORY_MODE)
        if os.path.lexists(output):
            raise ChromosomeAssignmentError(
                f"Output appeared during the run; refusing overwrite: {output.name}"
            )
        _rename_no_replace(staging, output)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    finally:
        _release_output_lock(lock_descriptor, lock_path, lock_inode)

    return ChromosomeAssignmentResult(
        output_dir=output,
        status=status,
        publication_gate=publication_gate,
        failure_states=ordered_failures,
        final_map_row_count=len(final_rows),
    )


__all__ = [
    "ChromosomeAssignmentError",
    "ChromosomeAssignmentResult",
    "HomologyPolicy",
    "JCVI_COLUMNS",
    "NUCLEOTIDE_COLUMNS",
    "assign_chromosome_homology",
    "load_homology_policy",
    "read_score_matrix",
    "solve_score_matrix",
]
