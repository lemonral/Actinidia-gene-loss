#!/usr/bin/env python3
"""Build validated nucleotide homology matrices after a capacity gate passes.

The controller waits for a declared upstream state to become ``PASS`` so that
its single matrix-building process does not exceed the frozen scientific-worker
ceiling.  It then waits for each immutable minimap bundle, and builds HY4A and
HY4P matrices sequentially.  Existing outputs, unsafe paths, duplicate units,
changed manifests, and non-PASS upstream bundles fail closed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


HEADER = (
    "unit",
    "target_scope_id",
    "target_genome",
    "target_asset_registry",
    "minimap_bundle_dir",
    "matrix_root",
)


class QueueError(RuntimeError):
    """Raised when the matrix queue cannot continue without weakening a gate."""


@dataclass(frozen=True)
class Entry:
    unit: str
    target_scope_id: str
    target_genome: Path
    target_asset_registry: Path
    minimap_bundle_dir: Path
    matrix_root: Path


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require_file(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink() or resolved.stat().st_size == 0:
        raise QueueError(f"{label} is missing, empty, or a symlink: {resolved}")
    return resolved


def safe_relative(raw: str, label: str) -> Path:
    value = raw.strip()
    path = Path(value)
    if not value or path.is_absolute() or ".." in path.parts:
        raise QueueError(f"{label} must be a safe non-empty relative path: {raw!r}")
    return path


def resolve_beneath(root: Path, relative: Path, label: str) -> Path:
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise QueueError(f"{label} escapes data root: {relative}") from error
    return candidate


def read_json(path: Path, label: str) -> dict[str, object]:
    try:
        payload = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=lambda value: (_ for _ in ()).throw(
                QueueError(f"non-finite JSON constant in {label}: {value}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise QueueError(f"cannot read {label}: {error}") from error
    if not isinstance(payload, dict):
        raise QueueError(f"{label} is not a JSON object")
    return payload


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(temporary, path)


def read_manifest(path: Path, data_root: Path) -> list[Entry]:
    manifest = require_file(path, "matrix queue manifest")
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != HEADER:
            raise QueueError(f"manifest header must be exactly {HEADER}")
        rows = list(reader)
    if not rows:
        raise QueueError("matrix queue manifest is empty")
    entries: list[Entry] = []
    units: set[str] = set()
    roots: set[Path] = set()
    for line_number, row in enumerate(rows, start=2):
        unit = row["unit"].strip()
        scope = row["target_scope_id"].strip()
        if not unit or any(character.isspace() for character in unit):
            raise QueueError(f"line {line_number}: unsafe unit {unit!r}")
        if not scope or any(character.isspace() for character in scope):
            raise QueueError(f"line {line_number}: unsafe target scope {scope!r}")
        if unit in units:
            raise QueueError(f"line {line_number}: duplicate unit {unit}")
        units.add(unit)
        paths = {
            column: resolve_beneath(
                data_root,
                safe_relative(row[column], f"line {line_number} {column}"),
                column,
            )
            for column in HEADER[2:]
        }
        matrix_root = paths["matrix_root"]
        if matrix_root in roots:
            raise QueueError(f"line {line_number}: duplicate matrix root")
        roots.add(matrix_root)
        entries.append(
            Entry(
                unit=unit,
                target_scope_id=scope,
                target_genome=paths["target_genome"],
                target_asset_registry=paths["target_asset_registry"],
                minimap_bundle_dir=paths["minimap_bundle_dir"],
                matrix_root=matrix_root,
            )
        )
    return entries


def wait_for_capacity_gate(path: Path, poll_seconds: int) -> dict[str, object]:
    while True:
        payload = read_json(path, "capacity-prerequisite state")
        status = payload.get("status")
        if status == "PASS":
            return payload
        if status not in {"running", "starting"}:
            raise QueueError(f"capacity prerequisite ended with status {status!r}")
        time.sleep(poll_seconds)


def wait_for_bundle(entry: Entry, poll_seconds: int) -> Path:
    validation = entry.minimap_bundle_dir / "bundle_validation.json"
    while not validation.is_file():
        status_path = entry.minimap_bundle_dir / "status.json"
        if status_path.is_file():
            status = read_json(status_path, f"minimap status for {entry.unit}").get("status")
            if status in {"failed", "interrupted"}:
                raise QueueError(f"minimap bundle for {entry.unit} ended with {status!r}")
        time.sleep(poll_seconds)
    payload = read_json(validation, f"bundle validation for {entry.unit}")
    if (
        payload.get("status") != "PASS"
        or payload.get("workflow") != "bidirectional_chromosome_minimap_bundle_validation"
        or payload.get("unit") != entry.unit
    ):
        raise QueueError(f"minimap bundle is not the expected PASS for {entry.unit}")
    return require_file(validation, f"bundle validation for {entry.unit}")


def run_builder(
    *,
    entry: Entry,
    role: str,
    python: Path,
    builder: Path,
    bundle: Path,
    reference_genome: Path,
    reference_assets: Path,
    reference_maps: Path,
    parameters: Path,
) -> dict[str, object]:
    suffix = "hy4a" if role == "nucleotide_hy4a" else "hy4p"
    forward = entry.minimap_bundle_dir / f"target_to_{suffix}.paf"
    reverse = entry.minimap_bundle_dir / f"{suffix}_to_target.paf"
    output = entry.matrix_root / role
    if output.exists():
        raise QueueError(f"refusing to overwrite matrix output: {output}")
    command = [
        str(python),
        str(builder),
        "--assembly-unit-id", entry.unit,
        "--target-scope-id", entry.target_scope_id,
        "--matrix-role", role,
        "--bundle-validation", str(bundle),
        "--forward-paf", str(require_file(forward, f"{entry.unit} {role} forward PAF")),
        "--reverse-paf", str(require_file(reverse, f"{entry.unit} {role} reverse PAF")),
        "--target-genome", str(require_file(entry.target_genome, f"{entry.unit} genome")),
        "--reference-genome", str(reference_genome),
        "--parameters", str(parameters),
        "--target-asset-registry", str(
            require_file(entry.target_asset_registry, f"{entry.unit} target registry")
        ),
        "--reference-asset-registry", str(reference_assets),
        "--reference-chromosome-map-registry", str(reference_maps),
        "--output-dir", str(output),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise QueueError(
            f"matrix builder failed for {entry.unit}/{role}: {completed.stderr.strip()}"
        )
    provenance = require_file(output / f"{role}.provenance.json", "matrix provenance")
    matrix = require_file(output / f"{role}.tsv", "matrix")
    payload = read_json(provenance, f"{entry.unit} {role} provenance")
    if (
        payload.get("status") != "PASS"
        or payload.get("assembly_unit_id") != entry.unit
        or payload.get("target_scope_id") != entry.target_scope_id
        or payload.get("matrix_role") != role
    ):
        raise QueueError(f"matrix provenance closure failed for {entry.unit}/{role}")
    return {
        "role": role,
        "matrix_sha256": sha256(matrix),
        "provenance_sha256": sha256(provenance),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--builder", required=True, type=Path)
    parser.add_argument("--parameters", required=True, type=Path)
    parser.add_argument("--reference-assets", required=True, type=Path)
    parser.add_argument("--reference-maps", required=True, type=Path)
    parser.add_argument("--hy4a-genome", required=True, type=Path)
    parser.add_argument("--hy4p-genome", required=True, type=Path)
    parser.add_argument("--capacity-prerequisite-state", required=True, type=Path)
    parser.add_argument("--state", required=True, type=Path)
    parser.add_argument("--lock", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=30)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    state: dict[str, object] | None = None
    try:
        if args.poll_seconds < 10 or args.poll_seconds > 60:
            raise QueueError("--poll-seconds must be between 10 and 60")
        args.data_root = args.data_root.expanduser().resolve()
        args.manifest = require_file(args.manifest, "matrix queue manifest")
        args.python = require_file(args.python, "Python executable")
        args.builder = require_file(args.builder, "nucleotide matrix builder")
        args.parameters = require_file(args.parameters, "analysis parameters")
        args.reference_assets = require_file(args.reference_assets, "reference assets")
        args.reference_maps = require_file(args.reference_maps, "reference maps")
        args.hy4a_genome = require_file(args.hy4a_genome, "HY4A genome")
        args.hy4p_genome = require_file(args.hy4p_genome, "HY4P genome")
        args.capacity_prerequisite_state = require_file(
            args.capacity_prerequisite_state, "capacity-prerequisite state"
        )
        args.state = args.state.expanduser().resolve()
        args.lock = args.lock.expanduser().resolve()
        if args.state.exists():
            raise QueueError(f"refusing existing state: {args.state}")
        args.lock.parent.mkdir(parents=True, exist_ok=True)
        descriptor = os.open(args.lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(f"pid\t{os.getpid()}\nmanifest_sha256\t{sha256(args.manifest)}\n")
        entries = read_manifest(args.manifest, args.data_root)
        state = {
            "schema_version": 1,
            "workflow": "sequential_nucleotide_homology_matrix_queue",
            "status": "waiting_for_capacity_prerequisite",
            "started_at_utc": utc_now(),
            "manifest_sha256": sha256(args.manifest),
            "expected_units": [entry.unit for entry in entries],
            "completed": [],
        }
        atomic_json(args.state, state)
        prerequisite = wait_for_capacity_gate(
            args.capacity_prerequisite_state, args.poll_seconds
        )
        state["capacity_prerequisite_sha256"] = sha256(args.capacity_prerequisite_state)
        state["capacity_prerequisite_status"] = prerequisite["status"]
        state["status"] = "running"
        atomic_json(args.state, state)
        for entry in entries:
            bundle = wait_for_bundle(entry, args.poll_seconds)
            results = []
            for role, reference in (
                ("nucleotide_hy4a", args.hy4a_genome),
                ("nucleotide_hy4p", args.hy4p_genome),
            ):
                results.append(
                    run_builder(
                        entry=entry,
                        role=role,
                        python=args.python,
                        builder=args.builder,
                        bundle=bundle,
                        reference_genome=reference,
                        reference_assets=args.reference_assets,
                        reference_maps=args.reference_maps,
                        parameters=args.parameters,
                    )
                )
            completed = state["completed"]
            assert isinstance(completed, list)
            completed.append(
                {
                    "unit": entry.unit,
                    "target_scope_id": entry.target_scope_id,
                    "bundle_validation_sha256": sha256(bundle),
                    "matrices": results,
                    "finished_at_utc": utc_now(),
                }
            )
            atomic_json(args.state, state)
        if len(state["completed"]) != len(entries):
            raise QueueError("completed units do not close against the manifest")
        state.update({"status": "PASS", "finished_at_utc": utc_now()})
        atomic_json(args.state, state)
        return 0
    except (OSError, subprocess.SubprocessError, QueueError) as error:
        if state is not None:
            state.update({"status": "failed", "finished_at_utc": utc_now(), "error": str(error)})
            atomic_json(args.state, state)
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
