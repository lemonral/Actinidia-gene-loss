#!/usr/bin/env python3
"""Sequentially materialize validated HY4A/HY4P JCVI matrices for new units."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


class QueueError(RuntimeError):
    pass


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink() or resolved.stat().st_size == 0:
        raise QueueError(f"missing, empty, or symlink file: {resolved}")
    return resolved


def under(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute() or ".." in Path(relative).parts:
        raise QueueError(f"unsafe data-root-relative path: {relative!r}")
    resolved = (root / relative).resolve()
    if not resolved.is_relative_to(root):
        raise QueueError(f"path escapes data root: {relative!r}")
    return resolved


def write_state(path: Path, document: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def load_pass_state(path: Path, label: str, expected_units: set[str]) -> dict[str, object]:
    document = json.loads(regular(path).read_text(encoding="utf-8"))
    if document.get("status") not in {"PASS", "completed"}:
        raise QueueError(f"{label} prerequisite is not complete/PASS")
    observed = {row.get("unit") for row in document.get("completed", [])}
    if not expected_units.issubset(observed):
        raise QueueError(f"{label} prerequisite lacks expected units")
    return document


def validate_role_output(path: Path, unit: str, scope: str, role: str) -> dict[str, str]:
    matrix = regular(path / f"{role}.tsv")
    provenance = regular(path / f"{role}.provenance.json")
    upstream = regular(path / f"{role}.upstream_validation.json")
    audit = regular(path / f"{role}.input_audit.json")
    document = json.loads(provenance.read_text(encoding="utf-8"))
    if (
        document.get("status") != "PASS"
        or document.get("assembly_unit_id") != unit
        or document.get("target_scope_id") != scope
        or document.get("matrix_role") != role
        or document.get("matrix", {}).get("sha256") != sha256(matrix)
    ):
        raise QueueError(f"invalid published matrix bundle: {unit}/{role}")
    return {
        "matrix_sha256": sha256(matrix),
        "provenance_sha256": sha256(provenance),
        "upstream_validation_sha256": sha256(upstream),
        "input_audit_sha256": sha256(audit),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--builder", required=True, type=Path)
    parser.add_argument("--parameters", required=True, type=Path)
    parser.add_argument("--reference-assets", required=True, type=Path)
    parser.add_argument("--reference-maps", required=True, type=Path)
    parser.add_argument("--hy4a-protein", required=True, type=Path)
    parser.add_argument("--hy4a-bed", required=True, type=Path)
    parser.add_argument("--hy4a-manifest", required=True, type=Path)
    parser.add_argument("--hy4p-protein", required=True, type=Path)
    parser.add_argument("--hy4p-bed", required=True, type=Path)
    parser.add_argument("--hy4p-manifest", required=True, type=Path)
    parser.add_argument("--jcvi-prerequisite-state", required=True, type=Path)
    parser.add_argument("--nucleotide-prerequisite-state", required=True, type=Path)
    parser.add_argument("--state", required=True, type=Path)
    parser.add_argument("--lock", required=True, type=Path)
    args = parser.parse_args()
    lock_acquired = False
    state: dict[str, object] | None = None
    try:
        root = args.data_root.expanduser().resolve()
        manifest = regular(args.manifest)
        with manifest.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        required = {
            "assembly_unit_id", "target_scope_id", "target_alias", "target_protein",
            "target_gff", "target_bed", "target_asset_registry", "jcvi_run_root", "matrix_root",
        }
        if not rows or not required.issubset(rows[0]):
            raise QueueError("manifest is empty or lacks required columns")
        units = [row["assembly_unit_id"] for row in rows]
        if len(units) != len(set(units)):
            raise QueueError("duplicate assembly units in manifest")
        expected = set(units)
        jcvi_state = load_pass_state(args.jcvi_prerequisite_state, "JCVI", expected)
        nucleotide_state = load_pass_state(args.nucleotide_prerequisite_state, "nucleotide", expected)

        args.lock.parent.mkdir(parents=True, exist_ok=True)
        lock_fd = os.open(args.lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        os.write(lock_fd, f"{os.getpid()}\n".encode())
        os.close(lock_fd)
        lock_acquired = True
        state = {
            "schema_version": 1,
            "workflow": "sequential_jcvi_homology_matrix_queue",
            "status": "running",
            "started_at_utc": now(),
            "manifest_sha256": sha256(manifest),
            "jcvi_prerequisite_sha256": sha256(args.jcvi_prerequisite_state),
            "nucleotide_prerequisite_sha256": sha256(args.nucleotide_prerequisite_state),
            "completed": [],
        }
        write_state(args.state, state)
        for row in rows:
            unit, scope, alias = (
                row["assembly_unit_id"], row["target_scope_id"], row["target_alias"]
            )
            run_root = under(root, row["jcvi_run_root"])
            matrix_root = under(root, row["matrix_root"])
            role_results = []
            for role, slot, protein, bed, canonical in (
                ("jcvi_hy4a", "hy4a", args.hy4a_protein, args.hy4a_bed, args.hy4a_manifest),
                ("jcvi_hy4p", "hy4p", args.hy4p_protein, args.hy4p_bed, args.hy4p_manifest),
            ):
                output = matrix_root / role
                if not output.exists():
                    command = [
                        str(regular(args.python)), str(regular(args.builder)),
                        "--assembly-unit-id", unit, "--target-scope-id", scope,
                        "--target-alias", alias, "--matrix-role", role,
                        "--forward-run-dir", str(run_root / f"{slot}_reference_target_query"),
                        "--reverse-run-dir", str(run_root / f"target_reference_{slot}_query"),
                        "--target-protein", str(regular(under(root, row["target_protein"]))),
                        "--target-gff", str(regular(under(root, row["target_gff"]))),
                        "--target-bed", str(regular(under(root, row["target_bed"]))),
                        "--reference-protein", str(regular(protein)),
                        "--reference-bed", str(regular(bed)),
                        "--reference-canonical-manifest", str(regular(canonical)),
                        "--parameters", str(regular(args.parameters)),
                        "--target-asset-registry", str(regular(under(root, row["target_asset_registry"]))),
                        "--reference-asset-registry", str(regular(args.reference_assets)),
                        "--reference-chromosome-map-registry", str(regular(args.reference_maps)),
                        "--output-dir", str(output),
                    ]
                    subprocess.run(command, check=True)
                hashes = validate_role_output(output, unit, scope, role)
                role_results.append({"role": role, **hashes})
            state["completed"].append({"unit": unit, "roles": role_results, "finished_at_utc": now()})
            write_state(args.state, state)
        state["status"] = "PASS"
        state["finished_at_utc"] = now()
        state["jcvi_prerequisite_status"] = jcvi_state["status"]
        state["nucleotide_prerequisite_status"] = nucleotide_state["status"]
        write_state(args.state, state)
        print(f"PASS\t{len(rows)} units")
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, QueueError, subprocess.CalledProcessError) as error:
        if args.state.parent.exists():
            failure = state or {
                "schema_version": 1,
                "workflow": "sequential_jcvi_homology_matrix_queue",
                "completed": [],
            }
            failure.update({"status": "FAIL", "error": str(error), "failed_at_utc": now()})
            write_state(args.state, failure)
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if lock_acquired:
            try:
                args.lock.unlink()
            except OSError:
                pass


if __name__ == "__main__":
    sys.exit(main())
