"""Immutable input snapshots and provenance gates for chromosome assignment."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import re
import stat
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Mapping, Sequence

from .io_utils import natural_key


PROVENANCE_SCHEMA_VERSION = "1.0.0"
UPSTREAM_REPORT_SCHEMA_VERSION = "1.0.0"
UPSTREAM_WORKFLOW_CONTRACT = {
    "nucleotide": ("chromosome_nucleotide_matrix_builder", "1.0.0"),
    "jcvi": ("chromosome_jcvi_matrix_builder", "1.0.0"),
}
MAX_SNAPSHOT_BYTES = 64 * 1024 * 1024
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
INTEGER_RE = re.compile(r"^(?:0|[1-9][0-9]*)$")
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:+-]{0,127}$")
MAX_INTEGER_DIGITS = 32

REFERENCE_ASSET_COLUMNS = (
    "reference_id",
    "biological_species",
    "individual_or_cultivar",
    "haplome",
    "role",
    "asset_role",
    "file_name",
    "bytes",
    "sha256",
    "source_url",
    "source_commit",
    "status",
    "notes",
)
REFERENCE_MAP_COLUMNS = (
    "reference_map_id",
    "reference_id",
    "reference_role",
    "reference_scope_id",
    "reference_chromosome",
    "canonical_chromosome",
    "status",
)
TARGET_ASSET_COLUMNS = (
    "assembly_unit_id",
    "target_scope_id",
    "asset_role",
    "file_name",
    "bytes",
    "sha256",
    "status",
)
TARGET_ASSET_ROLES = ("genome", "gff", "protein")
REFERENCE_ASSET_ROLES = ("genome", "gff", "protein", "cds")
ASSET_BINDING_KEYS = frozenset(("basename", "bytes", "sha256"))
PROVENANCE_KEYS = frozenset(
    (
        "schema_version",
        "status",
        "assembly_unit_id",
        "target_scope_id",
        "matrix_role",
        "matrix_kind",
        "reference_slot",
        "reference_id",
        "reference_role",
        "reference_map_id",
        "matrix",
        "target_asset_registry",
        "target_assets",
        "reference_asset_registry",
        "reference_assets",
        "reference_chromosome_map_registry",
        "generation_parameters",
        "upstream_validation_report",
    )
)
UPSTREAM_REPORT_KEYS = frozenset(
    (
        "schema_version",
        "workflow",
        "workflow_version",
        "status",
        "assembly_unit_id",
        "target_scope_id",
        "matrix_role",
        "matrix_kind",
        "reference_slot",
        "reference_id",
        "reference_role",
        "reference_map_id",
        "matrix_sha256",
        "target_asset_registry_sha256",
        "reference_asset_registry_sha256",
        "reference_chromosome_map_registry_sha256",
        "generation_parameters",
        "checks",
    )
)
COMMON_UPSTREAM_CHECKS = frozenset(
    (
        "target_asset_registry_reconciled",
        "reference_asset_registry_reconciled",
        "reference_chromosome_map_reconciled",
        "matrix_schema_reconciled",
        "complete_cartesian_matrix",
        "chromosome_denominators_reconciled",
        "matrix_arithmetic_reconciled",
        "matrix_checksum_reconciled",
    )
)
NUCLEOTIDE_UPSTREAM_CHECKS = COMMON_UPSTREAM_CHECKS | frozenset(
    (
        "bidirectional_paf_inputs_reconciled",
        "required_paf_tags_reconciled",
        "primary_alignment_filter_reconciled",
        "mapq_filter_reconciled",
        "block_length_filter_reconciled",
        "divergence_filter_reconciled",
        "query_interval_union_reconciled",
        "reference_interval_union_reconciled",
    )
)
JCVI_UPSTREAM_CHECKS = COMMON_UPSTREAM_CHECKS | frozenset(
    (
        "bidirectional_anchor_inputs_reconciled",
        "bed_protein_identity_reconciled",
        "eligible_gene_denominators_reconciled",
        "unique_anchor_counts_reconciled",
    )
)


class ChromosomeProvenanceError(RuntimeError):
    """Raised when an input snapshot or provenance binding is invalid."""


@dataclass(frozen=True)
class FileSnapshot:
    """One immutable byte snapshot used for both hashing and parsing."""

    path: Path
    data: bytes
    sha256: str
    device: int
    inode: int
    size: int
    mtime_ns: int

    @property
    def basename(self) -> str:
        return self.path.name

    @property
    def inode_key(self) -> tuple[int, int]:
        return self.device, self.inode

    def public_binding(self) -> dict[str, str | int]:
        return {"basename": self.basename, "bytes": self.size, "sha256": self.sha256}


@dataclass(frozen=True)
class AssetBinding:
    """Path-free expected identity of one target or reference asset."""

    basename: str
    bytes: int
    sha256: str

    def as_dict(self) -> dict[str, str | int]:
        return {"basename": self.basename, "bytes": self.bytes, "sha256": self.sha256}


@dataclass(frozen=True)
class ReferenceRegistries:
    """Frozen reference assets and chromosome truth maps."""

    asset_snapshot: FileSnapshot
    map_snapshot: FileSnapshot
    roles: Mapping[str, str]
    assets: Mapping[str, Mapping[str, AssetBinding]]
    chromosome_maps: Mapping[str, Mapping[str, str]]
    map_ids: Mapping[str, str]
    scope_ids: Mapping[str, str]


@dataclass(frozen=True)
class TargetRegistry:
    """Frozen target asset identities keyed by assembly unit."""

    snapshot: FileSnapshot
    assets: Mapping[tuple[str, str], Mapping[str, AssetBinding]]


@dataclass(frozen=True)
class MatrixProvenance:
    """Validated matrix sidecar and its upstream PASS report."""

    sidecar_snapshot: FileSnapshot
    upstream_report_snapshot: FileSnapshot
    assembly_unit_id: str
    target_scope_id: str
    matrix_role: str
    matrix_kind: str
    reference_slot: str
    reference_id: str
    reference_role: str
    reference_map_id: str
    generation_parameters: Mapping[str, Any]


def _total_key(value: str) -> tuple[list[object], str]:
    return natural_key(value), value


def parse_exact_int(value: str, *, label: str, positive: bool = False) -> int:
    """Parse an integer lexically without any float round trip."""

    if len(value) > MAX_INTEGER_DIGITS or not INTEGER_RE.fullmatch(value):
        raise ChromosomeProvenanceError(f"{label} must be an exact base-10 integer")
    try:
        parsed = int(value)
    except ValueError as error:
        raise ChromosomeProvenanceError(
            f"{label} must be an exact base-10 integer"
        ) from error
    if positive and parsed == 0:
        raise ChromosomeProvenanceError(f"{label} must be positive")
    return parsed


def _read_fd_all(file_descriptor: int, *, maximum_bytes: int) -> bytes:
    parts: list[bytes] = []
    observed = 0
    while True:
        block = os.read(file_descriptor, min(1 << 20, maximum_bytes + 1 - observed))
        if not block:
            return b"".join(parts)
        parts.append(block)
        observed += len(block)
        if observed > maximum_bytes:
            raise ChromosomeProvenanceError(
                f"Immutable input exceeds the {maximum_bytes}-byte snapshot limit"
            )


def capture_snapshot(path: str | Path) -> FileSnapshot:
    """Read one regular file once and bind its hash to exactly those parsed bytes."""

    source = Path(path).expanduser()
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NONBLOCK"):
        flags |= os.O_NONBLOCK
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(source, flags)
    except (OSError, ValueError) as error:
        raise ChromosomeProvenanceError(
            f"Cannot open immutable input {source.name}: {error}"
        ) from error
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_size <= 0:
            raise ChromosomeProvenanceError(
                f"Input must be one non-empty regular file: {source.name}"
            )
        if before.st_size > MAX_SNAPSHOT_BYTES:
            raise ChromosomeProvenanceError(
                f"Input exceeds the {MAX_SNAPSHOT_BYTES}-byte snapshot limit: {source.name}"
            )
        data = _read_fd_all(descriptor, maximum_bytes=MAX_SNAPSHOT_BYTES)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    try:
        path_state = source.stat(follow_symlinks=False)
    except (OSError, ValueError) as error:
        raise ChromosomeProvenanceError(
            f"Input disappeared during snapshot: {source.name}: {error}"
        ) from error
    before_signature = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    after_signature = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    path_signature = (
        path_state.st_dev,
        path_state.st_ino,
        path_state.st_size,
        path_state.st_mtime_ns,
    )
    if before_signature != after_signature or before_signature != path_signature:
        raise ChromosomeProvenanceError(
            f"Input changed while its immutable snapshot was read: {source.name}"
        )
    if len(data) != before.st_size:
        raise ChromosomeProvenanceError(
            f"Input byte count changed while reading: {source.name}"
        )
    return FileSnapshot(
        path=source,
        data=data,
        sha256=hashlib.sha256(data).hexdigest(),
        device=before.st_dev,
        inode=before.st_ino,
        size=before.st_size,
        mtime_ns=before.st_mtime_ns,
    )


def verify_snapshot(snapshot: FileSnapshot) -> None:
    """Require the input path to retain the exact captured inode and bytes."""

    current = capture_snapshot(snapshot.path)
    if (
        current.inode_key != snapshot.inode_key
        or current.size != snapshot.size
        or current.mtime_ns != snapshot.mtime_ns
        or current.sha256 != snapshot.sha256
    ):
        raise ChromosomeProvenanceError(
            f"Input changed after snapshot; refusing publication: {snapshot.basename}"
        )


def _decode_text(snapshot: FileSnapshot) -> str:
    try:
        return snapshot.data.decode("utf-8", errors="strict")
    except UnicodeError as error:
        raise ChromosomeProvenanceError(
            f"Input is not strict UTF-8 text: {snapshot.basename}: {error}"
        ) from error


def _reject_json_constant(value: str) -> None:
    raise ChromosomeProvenanceError(f"JSON non-finite numeric constant is forbidden: {value}")


def _unique_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ChromosomeProvenanceError(f"JSON object repeats key {key!r}")
        result[key] = value
    return result


def _validate_json_finite(value: Any, label: str = "JSON") -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise ChromosomeProvenanceError(f"{label} contains a non-finite number")
    if isinstance(value, Decimal) and not value.is_finite():
        raise ChromosomeProvenanceError(f"{label} contains a non-finite number")
    if isinstance(value, dict):
        for key, child in value.items():
            if not isinstance(key, str):
                raise ChromosomeProvenanceError(f"{label} contains a non-string object key")
            _validate_json_finite(child, f"{label}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _validate_json_finite(child, f"{label}[{index}]")


def parse_json_snapshot(snapshot: FileSnapshot) -> Mapping[str, Any]:
    """Parse strict JSON from the same bytes that supplied its SHA-256."""

    try:
        parsed = json.loads(
            _decode_text(snapshot),
            parse_constant=_reject_json_constant,
            parse_float=Decimal,
            parse_int=lambda value: parse_exact_int(value, label="JSON integer"),
            object_pairs_hook=_unique_json_object,
        )
    except (InvalidOperation, json.JSONDecodeError, TypeError, ValueError) as error:
        raise ChromosomeProvenanceError(
            f"Invalid JSON in {snapshot.basename}: {error}"
        ) from error
    _validate_json_finite(parsed, snapshot.basename)
    if not isinstance(parsed, dict):
        raise ChromosomeProvenanceError(f"{snapshot.basename}: top-level JSON must be an object")
    return parsed


def _read_exact_tsv(
    snapshot: FileSnapshot, expected_columns: Sequence[str]
) -> list[dict[str, str]]:
    try:
        reader = csv.DictReader(io.StringIO(_decode_text(snapshot)), delimiter="\t")
        observed = tuple(reader.fieldnames or ())
        if observed != tuple(expected_columns):
            raise ChromosomeProvenanceError(
                f"{snapshot.basename}: exact schema required; expected {list(expected_columns)}, "
                f"found {list(observed)}"
            )
        rows: list[dict[str, str]] = []
        for line, row in enumerate(reader, start=2):
            if None in row:
                raise ChromosomeProvenanceError(
                    f"{snapshot.basename}:{line}: row has extra fields"
                )
            cleaned = {key: (value or "").strip() for key, value in row.items()}
            if not any(cleaned.values()):
                continue
            missing = [key for key, value in cleaned.items() if value == ""]
            if missing:
                raise ChromosomeProvenanceError(
                    f"{snapshot.basename}:{line}: empty required values: {','.join(missing)}"
                )
            rows.append(cleaned)
    except csv.Error as error:
        raise ChromosomeProvenanceError(
            f"Cannot parse TSV {snapshot.basename}: {error}"
        ) from error
    if not rows:
        raise ChromosomeProvenanceError(f"{snapshot.basename}: table contains no rows")
    return rows


def _validate_safe_basename(value: str, label: str) -> None:
    if value in {"", ".", ".."} or Path(value).name != value or "/" in value or "\\" in value:
        raise ChromosomeProvenanceError(f"{label} must be one path-free basename")
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        raise ChromosomeProvenanceError(f"{label} contains a control character")


def _validate_sha256(value: str, label: str) -> None:
    if not isinstance(value, str) or not SHA256_RE.fullmatch(value):
        raise ChromosomeProvenanceError(f"{label} must be a lowercase SHA-256")


def _asset_from_values(
    *, basename: str, bytes_value: str | int, sha256: str, label: str
) -> AssetBinding:
    if not isinstance(basename, str):
        raise ChromosomeProvenanceError(f"{label}.basename must be a string")
    _validate_safe_basename(basename, f"{label}.basename")
    if isinstance(bytes_value, bool):
        raise ChromosomeProvenanceError(f"{label}.bytes must be a positive integer")
    if isinstance(bytes_value, int):
        size = bytes_value
    elif isinstance(bytes_value, str):
        size = parse_exact_int(bytes_value, label=f"{label}.bytes", positive=True)
    else:
        raise ChromosomeProvenanceError(f"{label}.bytes must be a positive integer")
    if size <= 0:
        raise ChromosomeProvenanceError(f"{label}.bytes must be positive")
    _validate_sha256(sha256, f"{label}.sha256")
    return AssetBinding(basename=basename, bytes=size, sha256=sha256)


def _binding_from_json(value: Any, label: str) -> AssetBinding:
    if not isinstance(value, dict) or set(value) != ASSET_BINDING_KEYS:
        raise ChromosomeProvenanceError(
            f"{label} must have exactly basename, bytes, and sha256"
        )
    return _asset_from_values(
        basename=value["basename"],
        bytes_value=value["bytes"],
        sha256=value["sha256"],
        label=label,
    )


def read_reference_registries(
    asset_registry: FileSnapshot,
    chromosome_map_registry: FileSnapshot,
    *,
    coordinate_reference: str,
    confirmation_reference: str,
    expected_chromosomes: int,
) -> ReferenceRegistries:
    """Read exact frozen reference assets and chromosome truth maps."""

    if coordinate_reference == confirmation_reference:
        raise ChromosomeProvenanceError("Coordinate and confirmation references must differ")
    asset_rows = _read_exact_tsv(asset_registry, REFERENCE_ASSET_COLUMNS)
    assets: dict[str, dict[str, AssetBinding]] = {}
    roles: dict[str, str] = {}
    for line, row in enumerate(asset_rows, start=2):
        reference_id = row["reference_id"]
        if reference_id not in {coordinate_reference, confirmation_reference}:
            continue
        if row["status"] != "verified_legacy_asset":
            raise ChromosomeProvenanceError(
                f"{asset_registry.basename}:{line}: reference asset is not verified"
            )
        role = row["role"]
        if reference_id in roles and roles[reference_id] != role:
            raise ChromosomeProvenanceError(
                f"{asset_registry.basename}:{line}: reference role changes within ID"
            )
        roles[reference_id] = role
        asset_role = row["asset_role"]
        if asset_role not in REFERENCE_ASSET_ROLES:
            raise ChromosomeProvenanceError(
                f"{asset_registry.basename}:{line}: unsupported reference asset role"
            )
        if asset_role in assets.setdefault(reference_id, {}):
            raise ChromosomeProvenanceError(
                f"{asset_registry.basename}:{line}: duplicate reference asset role"
            )
        assets[reference_id][asset_role] = _asset_from_values(
            basename=row["file_name"],
            bytes_value=row["bytes"],
            sha256=row["sha256"],
            label=f"{reference_id}.{asset_role}",
        )
    for reference_id, expected_role in (
        (coordinate_reference, "primary_coordinate_reference"),
        (confirmation_reference, "independent_confirmation_reference"),
    ):
        if roles.get(reference_id) != expected_role:
            raise ChromosomeProvenanceError(
                f"Reference {reference_id!r} does not have frozen role {expected_role!r}"
            )
        if set(assets.get(reference_id, {})) != set(REFERENCE_ASSET_ROLES):
            raise ChromosomeProvenanceError(
                f"Reference {reference_id!r} lacks the exact genome/GFF/protein/CDS asset set"
            )
    for asset_role in REFERENCE_ASSET_ROLES:
        if (
            assets[coordinate_reference][asset_role].sha256
            == assets[confirmation_reference][asset_role].sha256
        ):
            raise ChromosomeProvenanceError(
                f"Coordinate and confirmation references share {asset_role} bytes"
            )

    map_rows = _read_exact_tsv(chromosome_map_registry, REFERENCE_MAP_COLUMNS)
    chromosome_maps: dict[str, dict[str, str]] = {}
    map_ids: dict[str, str] = {}
    scope_ids: dict[str, str] = {}
    for line, row in enumerate(map_rows, start=2):
        reference_id = row["reference_id"]
        if reference_id not in {coordinate_reference, confirmation_reference}:
            continue
        if row["status"] != "verified":
            raise ChromosomeProvenanceError(
                f"{chromosome_map_registry.basename}:{line}: map row is not verified"
            )
        if row["reference_role"] != roles[reference_id]:
            raise ChromosomeProvenanceError(
                f"{chromosome_map_registry.basename}:{line}: reference role conflicts with asset registry"
            )
        for field, store in (
            ("reference_map_id", map_ids),
            ("reference_scope_id", scope_ids),
        ):
            value = row[field]
            if not SAFE_ID_RE.fullmatch(value):
                raise ChromosomeProvenanceError(
                    f"{chromosome_map_registry.basename}:{line}: unsafe {field}"
                )
            if reference_id in store and store[reference_id] != value:
                raise ChromosomeProvenanceError(
                    f"{chromosome_map_registry.basename}:{line}: {field} changes within reference"
                )
            store[reference_id] = value
        reference_chromosome = row["reference_chromosome"]
        canonical = row["canonical_chromosome"]
        if not SAFE_ID_RE.fullmatch(reference_chromosome):
            raise ChromosomeProvenanceError(
                f"{chromosome_map_registry.basename}:{line}: unsafe reference chromosome ID"
            )
        expected_canonical = {
            f"Chr{index:02d}" for index in range(1, expected_chromosomes + 1)
        }
        if canonical not in expected_canonical:
            raise ChromosomeProvenanceError(
                f"{chromosome_map_registry.basename}:{line}: invalid canonical chromosome"
            )
        target = chromosome_maps.setdefault(reference_id, {})
        if reference_chromosome in target:
            raise ChromosomeProvenanceError(
                f"{chromosome_map_registry.basename}:{line}: duplicate reference chromosome"
            )
        target[reference_chromosome] = canonical
    for reference_id in (coordinate_reference, confirmation_reference):
        mapping = chromosome_maps.get(reference_id, {})
        if len(mapping) != expected_chromosomes:
            raise ChromosomeProvenanceError(
                f"Reference {reference_id!r} must have exactly {expected_chromosomes} map rows"
            )
        if set(mapping.values()) != {
            f"Chr{index:02d}" for index in range(1, expected_chromosomes + 1)
        }:
            raise ChromosomeProvenanceError(
                f"Reference {reference_id!r} chromosome map is not a canonical bijection"
            )
    if map_ids[coordinate_reference] == map_ids[confirmation_reference]:
        raise ChromosomeProvenanceError("Coordinate and confirmation map IDs must differ")
    if scope_ids[coordinate_reference] == scope_ids[confirmation_reference]:
        raise ChromosomeProvenanceError("Coordinate and confirmation scope IDs must differ")
    return ReferenceRegistries(
        asset_snapshot=asset_registry,
        map_snapshot=chromosome_map_registry,
        roles=roles,
        assets=assets,
        chromosome_maps=chromosome_maps,
        map_ids=map_ids,
        scope_ids=scope_ids,
    )


def read_target_registry(snapshot: FileSnapshot) -> TargetRegistry:
    """Read a frozen assembly-unit-to-standardized-asset registry."""

    rows = _read_exact_tsv(snapshot, TARGET_ASSET_COLUMNS)
    assets: dict[tuple[str, str], dict[str, AssetBinding]] = {}
    for line, row in enumerate(rows, start=2):
        unit = row["assembly_unit_id"]
        if not SAFE_ID_RE.fullmatch(unit):
            raise ChromosomeProvenanceError(
                f"{snapshot.basename}:{line}: unsafe assembly_unit_id"
            )
        target_scope_id = row["target_scope_id"]
        if not SAFE_ID_RE.fullmatch(target_scope_id):
            raise ChromosomeProvenanceError(
                f"{snapshot.basename}:{line}: unsafe target_scope_id"
            )
        if row["status"] != "verified":
            raise ChromosomeProvenanceError(
                f"{snapshot.basename}:{line}: target asset is not verified"
            )
        role = row["asset_role"]
        if role not in TARGET_ASSET_ROLES:
            raise ChromosomeProvenanceError(
                f"{snapshot.basename}:{line}: unsupported target asset role"
            )
        target_key = (unit, target_scope_id)
        if role in assets.setdefault(target_key, {}):
            raise ChromosomeProvenanceError(
                f"{snapshot.basename}:{line}: duplicate target asset role"
            )
        assets[target_key][role] = _asset_from_values(
            basename=row["file_name"],
            bytes_value=row["bytes"],
            sha256=row["sha256"],
            label=f"{unit}.{role}",
        )
    for (unit, scope), unit_assets in assets.items():
        if set(unit_assets) != set(TARGET_ASSET_ROLES):
            raise ChromosomeProvenanceError(
                f"Target unit/scope {unit!r}/{scope!r} lacks the exact genome/GFF/protein asset set"
            )
        bindings = [unit_assets[role] for role in TARGET_ASSET_ROLES]
        if len({binding.basename for binding in bindings}) != len(bindings):
            raise ChromosomeProvenanceError(
                f"Target unit/scope {unit!r}/{scope!r} reuses one basename across asset roles"
            )
        if len({binding.sha256 for binding in bindings}) != len(bindings):
            raise ChromosomeProvenanceError(
                f"Target unit/scope {unit!r}/{scope!r} has copy-identical role assets"
            )
    if len(assets) != 1:
        raise ChromosomeProvenanceError(
            "Target asset registry must contain exactly one assembly-unit/target-scope pair"
        )
    return TargetRegistry(snapshot=snapshot, assets=assets)


def _require_typed_equal(observed: Any, expected: Any, label: str) -> None:
    if isinstance(expected, float) and isinstance(observed, Decimal):
        if (expected >= 0.0 and observed.is_signed()) or observed != Decimal(str(expected)):
            raise ChromosomeProvenanceError(
                f"{label} differs from the frozen generation parameter"
            )
        return
    if type(observed) is not type(expected):
        raise ChromosomeProvenanceError(f"{label} has the wrong JSON type")
    if isinstance(expected, dict):
        if set(observed) != set(expected):
            raise ChromosomeProvenanceError(f"{label} has an unexpected key set")
        for key in expected:
            _require_typed_equal(observed[key], expected[key], f"{label}.{key}")
    elif observed != expected:
        raise ChromosomeProvenanceError(
            f"{label} differs from the frozen generation parameter"
        )


def _validate_assets_json(
    value: Any,
    *,
    expected: Mapping[str, AssetBinding],
    expected_roles: Sequence[str],
    label: str,
) -> None:
    if not isinstance(value, dict) or set(value) != set(expected_roles):
        raise ChromosomeProvenanceError(f"{label} has the wrong asset-role set")
    for role in expected_roles:
        observed = _binding_from_json(value[role], f"{label}.{role}")
        if observed != expected[role]:
            raise ChromosomeProvenanceError(
                f"{label}.{role} conflicts with its frozen registry"
            )


def validate_matrix_provenance(
    *,
    sidecar: FileSnapshot,
    matrix: FileSnapshot,
    assembly_unit_id: str,
    target_scope_id: str,
    matrix_role: str,
    matrix_kind: str,
    reference_slot: str,
    reference_id: str,
    reference_role: str,
    reference_map_id: str,
    expected_generation_parameters: Mapping[str, Any],
    target_registry: TargetRegistry,
    reference_registries: ReferenceRegistries,
) -> MatrixProvenance:
    """Validate one exact sidecar and its checksum-bound upstream PASS report."""

    document = parse_json_snapshot(sidecar)
    if set(document) != PROVENANCE_KEYS:
        raise ChromosomeProvenanceError(
            f"{sidecar.basename}: provenance sidecar has an unexpected key set"
        )
    expected_scalars = {
        "schema_version": PROVENANCE_SCHEMA_VERSION,
        "status": "PASS",
        "assembly_unit_id": assembly_unit_id,
        "target_scope_id": target_scope_id,
        "matrix_role": matrix_role,
        "matrix_kind": matrix_kind,
        "reference_slot": reference_slot,
        "reference_id": reference_id,
        "reference_role": reference_role,
        "reference_map_id": reference_map_id,
    }
    for key, expected in expected_scalars.items():
        if type(document[key]) is not type(expected) or document[key] != expected:
            raise ChromosomeProvenanceError(
                f"{sidecar.basename}: {key} does not match the invoked matrix role"
            )
    matrix_binding = _binding_from_json(document["matrix"], f"{sidecar.basename}.matrix")
    if matrix_binding != AssetBinding(matrix.basename, matrix.size, matrix.sha256):
        raise ChromosomeProvenanceError(
            f"{sidecar.basename}: matrix binding does not match the parsed matrix bytes"
        )
    target_key = (assembly_unit_id, target_scope_id)
    if target_key not in target_registry.assets:
        raise ChromosomeProvenanceError(
            f"Target registry has no verified assets for {assembly_unit_id!r}"
        )
    target_registry_binding = _binding_from_json(
        document["target_asset_registry"], f"{sidecar.basename}.target_asset_registry"
    )
    if target_registry_binding != AssetBinding(
        target_registry.snapshot.basename,
        target_registry.snapshot.size,
        target_registry.snapshot.sha256,
    ):
        raise ChromosomeProvenanceError(
            f"{sidecar.basename}: target registry binding is not exact"
        )
    _validate_assets_json(
        document["target_assets"],
        expected=target_registry.assets[target_key],
        expected_roles=TARGET_ASSET_ROLES,
        label=f"{sidecar.basename}.target_assets",
    )
    reference_registry_binding = _binding_from_json(
        document["reference_asset_registry"],
        f"{sidecar.basename}.reference_asset_registry",
    )
    if reference_registry_binding != AssetBinding(
        reference_registries.asset_snapshot.basename,
        reference_registries.asset_snapshot.size,
        reference_registries.asset_snapshot.sha256,
    ):
        raise ChromosomeProvenanceError(
            f"{sidecar.basename}: reference asset registry binding is not exact"
        )
    map_registry_binding = _binding_from_json(
        document["reference_chromosome_map_registry"],
        f"{sidecar.basename}.reference_chromosome_map_registry",
    )
    if map_registry_binding != AssetBinding(
        reference_registries.map_snapshot.basename,
        reference_registries.map_snapshot.size,
        reference_registries.map_snapshot.sha256,
    ):
        raise ChromosomeProvenanceError(
            f"{sidecar.basename}: reference chromosome-map registry binding is not exact"
        )
    _validate_assets_json(
        document["reference_assets"],
        expected=reference_registries.assets[reference_id],
        expected_roles=REFERENCE_ASSET_ROLES,
        label=f"{sidecar.basename}.reference_assets",
    )
    if reference_registries.map_ids[reference_id] != reference_map_id:
        raise ChromosomeProvenanceError(
            f"{sidecar.basename}: reference_map_id conflicts with the frozen map registry"
        )
    _require_typed_equal(
        document["generation_parameters"],
        dict(expected_generation_parameters),
        f"{sidecar.basename}.generation_parameters",
    )

    report_binding = _binding_from_json(
        document["upstream_validation_report"],
        f"{sidecar.basename}.upstream_validation_report",
    )
    report_path = sidecar.path.parent / report_binding.basename
    report_snapshot = capture_snapshot(report_path)
    if report_binding != AssetBinding(
        report_snapshot.basename, report_snapshot.size, report_snapshot.sha256
    ):
        raise ChromosomeProvenanceError(
            f"{sidecar.basename}: upstream report checksum binding is not exact"
        )
    report = parse_json_snapshot(report_snapshot)
    if set(report) != UPSTREAM_REPORT_KEYS:
        raise ChromosomeProvenanceError(
            f"{report_snapshot.basename}: upstream report has an unexpected key set"
        )
    workflow, workflow_version = UPSTREAM_WORKFLOW_CONTRACT[matrix_kind]
    expected_report = {
        "schema_version": UPSTREAM_REPORT_SCHEMA_VERSION,
        "workflow": workflow,
        "workflow_version": workflow_version,
        "status": "PASS",
        "assembly_unit_id": assembly_unit_id,
        "target_scope_id": target_scope_id,
        "matrix_role": matrix_role,
        "matrix_kind": matrix_kind,
        "reference_slot": reference_slot,
        "reference_id": reference_id,
        "reference_role": reference_role,
        "reference_map_id": reference_map_id,
        "matrix_sha256": matrix.sha256,
        "target_asset_registry_sha256": target_registry.snapshot.sha256,
        "reference_asset_registry_sha256": reference_registries.asset_snapshot.sha256,
        "reference_chromosome_map_registry_sha256": (
            reference_registries.map_snapshot.sha256
        ),
    }
    for key, expected in expected_report.items():
        if type(report[key]) is not type(expected) or report[key] != expected:
            raise ChromosomeProvenanceError(
                f"{report_snapshot.basename}: {key} conflicts with the matrix sidecar"
            )
    _require_typed_equal(
        report["generation_parameters"],
        dict(expected_generation_parameters),
        f"{report_snapshot.basename}.generation_parameters",
    )
    checks = report["checks"]
    expected_checks = (
        NUCLEOTIDE_UPSTREAM_CHECKS
        if matrix_kind == "nucleotide"
        else JCVI_UPSTREAM_CHECKS
    )
    if not isinstance(checks, dict) or set(checks) != expected_checks:
        raise ChromosomeProvenanceError(
            f"{report_snapshot.basename}: upstream report has the wrong exact check set"
        )
    if any(type(value) is not bool or not value for value in checks.values()):
        raise ChromosomeProvenanceError(
            f"{report_snapshot.basename}: every upstream validation check must be true"
        )
    return MatrixProvenance(
        sidecar_snapshot=sidecar,
        upstream_report_snapshot=report_snapshot,
        assembly_unit_id=assembly_unit_id,
        target_scope_id=target_scope_id,
        matrix_role=matrix_role,
        matrix_kind=matrix_kind,
        reference_slot=reference_slot,
        reference_id=reference_id,
        reference_role=reference_role,
        reference_map_id=reference_map_id,
        generation_parameters=dict(expected_generation_parameters),
    )


def reject_duplicate_snapshots(
    snapshots: Mapping[str, FileSnapshot], *, label: str
) -> None:
    """Reject hardlinks and byte-identical role substitution as an extra defense."""

    inode_seen: dict[tuple[int, int], str] = {}
    hash_seen: dict[str, str] = {}
    for role in sorted(snapshots, key=_total_key):
        snapshot = snapshots[role]
        if snapshot.inode_key in inode_seen:
            raise ChromosomeProvenanceError(
                f"{label}: roles {inode_seen[snapshot.inode_key]!r} and {role!r} use one inode"
            )
        if snapshot.sha256 in hash_seen:
            raise ChromosomeProvenanceError(
                f"{label}: roles {hash_seen[snapshot.sha256]!r} and {role!r} have identical bytes"
            )
        inode_seen[snapshot.inode_key] = role
        hash_seen[snapshot.sha256] = role


__all__ = [
    "AssetBinding",
    "ChromosomeProvenanceError",
    "FileSnapshot",
    "MatrixProvenance",
    "PROVENANCE_SCHEMA_VERSION",
    "ReferenceRegistries",
    "TargetRegistry",
    "capture_snapshot",
    "parse_exact_int",
    "parse_json_snapshot",
    "read_reference_registries",
    "read_target_registry",
    "reject_duplicate_snapshots",
    "validate_matrix_provenance",
    "verify_snapshot",
]
