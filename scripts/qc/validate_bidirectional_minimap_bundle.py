#!/usr/bin/env python3
"""Validate a completed four-direction chromosome minimap2 bundle.

This is the immutable binding gate before PAF rows are converted to nucleotide
score matrices.  It checks the fixed workflow/version/roles, all current input
and output sizes and SHA-256 values, and detects files changing while hashed.
It does not interpret PAF alignments; that is the matrix builder's job.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path


EXPECTED_COMPARISONS = {
    "target_to_hy4a": ("target", "hy4a"),
    "hy4a_to_target": ("hy4a", "target"),
    "target_to_hy4p": ("target", "hy4p"),
    "hy4p_to_target": ("hy4p", "target"),
}
EXPECTED_ARGV = [
    "minimap2",
    "-x",
    "asm5",
    "--secondary=no",
    "-c",
    "--cs=long",
    "{reference_fasta}",
    "{query_fasta}",
]


class ValidationError(RuntimeError):
    pass


def stable_record(path: Path) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise ValidationError(f"not a regular non-symlink file: {path}")
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    after = path.stat()
    before_key = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    after_key = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if before_key != after_key:
        raise ValidationError(f"file changed while being hashed: {path}")
    return {"basename": path.name, "bytes": after.st_size, "sha256": digest.hexdigest()}


def require_binding(
    *, expected: object, observed: dict[str, object], label: str
) -> None:
    if not isinstance(expected, dict):
        raise ValidationError(f"status is missing {label} binding")
    for field in ("basename", "bytes", "sha256"):
        if expected.get(field) != observed.get(field):
            raise ValidationError(
                f"{label} {field} mismatch: status={expected.get(field)!r}, "
                f"current={observed.get(field)!r}"
            )


def atomic_json(path: Path, payload: dict[str, object], overwrite: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise ValidationError(f"refusing to overwrite existing output: {path}")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent, text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def validate(
    *,
    status_path: Path,
    target_genome: Path,
    hy4a_genome: Path,
    hy4p_genome: Path,
    expected_unit: str,
) -> dict[str, object]:
    try:
        status = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValidationError(f"cannot read status JSON: {error}") from error
    if status.get("schema_version") != 1:
        raise ValidationError("status schema_version must be 1")
    if status.get("workflow") != "bidirectional_chromosome_minimap":
        raise ValidationError("unexpected workflow")
    if status.get("status") != "completed" or not status.get("finished_at_utc"):
        raise ValidationError("workflow is not completed")
    if status.get("unit") != expected_unit:
        raise ValidationError("unit does not match --expected-unit")
    if status.get("minimap2_version") != "2.28-r1209":
        raise ValidationError("minimap2 version is not 2.28-r1209")
    if status.get("fixed_argv") != EXPECTED_ARGV:
        raise ValidationError("fixed minimap2 argv differs from the frozen command")

    input_paths = {
        "target": target_genome,
        "hy4a": hy4a_genome,
        "hy4p": hy4p_genome,
    }
    inputs: dict[str, dict[str, object]] = {}
    status_inputs = status.get("inputs")
    if not isinstance(status_inputs, dict) or set(status_inputs) != set(input_paths):
        raise ValidationError("status input roles are incomplete or unexpected")
    for role, path in input_paths.items():
        record = stable_record(path)
        require_binding(expected=status_inputs.get(role), observed=record, label=f"input {role}")
        inputs[role] = record

    comparisons = status.get("comparisons")
    if not isinstance(comparisons, dict) or set(comparisons) != set(EXPECTED_COMPARISONS):
        raise ValidationError("comparison roles are incomplete or unexpected")
    outputs: dict[str, dict[str, object]] = {}
    root = status_path.parent.resolve()
    for comparison, (query_role, reference_role) in EXPECTED_COMPARISONS.items():
        item = comparisons[comparison]
        if not isinstance(item, dict):
            raise ValidationError(f"{comparison}: status record is not an object")
        if item.get("exit_code") != 0 or not item.get("finished_at_utc"):
            raise ValidationError(f"{comparison}: comparison did not complete successfully")
        if item.get("query_role") != query_role or item.get("reference_role") != reference_role:
            raise ValidationError(f"{comparison}: query/reference roles are incorrect")
        comparison_outputs: dict[str, object] = {}
        for kind in ("paf", "stderr"):
            binding = item.get(kind)
            if not isinstance(binding, dict) or not isinstance(binding.get("basename"), str):
                raise ValidationError(f"{comparison}: missing {kind} binding")
            basename = str(binding["basename"])
            if Path(basename).name != basename:
                raise ValidationError(f"{comparison}: unsafe {kind} basename")
            path = root / basename
            current = stable_record(path)
            require_binding(expected=binding, observed=current, label=f"{comparison} {kind}")
            if kind == "paf" and current["bytes"] == 0:
                raise ValidationError(f"{comparison}: empty PAF")
            comparison_outputs[kind] = current
        outputs[comparison] = comparison_outputs

    return {
        "schema_version": 1,
        "status": "PASS",
        "workflow": "bidirectional_chromosome_minimap_bundle_validation",
        "unit": expected_unit,
        "minimap2_version": "2.28-r1209",
        "fixed_argv": EXPECTED_ARGV,
        "status_json": stable_record(status_path),
        "inputs": inputs,
        "comparisons": outputs,
        "checks": {
            "completed_status": True,
            "fixed_version": True,
            "fixed_argv": True,
            "exact_roles": True,
            "input_size_sha256_closure": True,
            "output_size_sha256_closure": True,
            "stable_files_during_hashing": True,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--target-genome", required=True, type=Path)
    parser.add_argument("--hy4a-genome", required=True, type=Path)
    parser.add_argument("--hy4p-genome", required=True, type=Path)
    parser.add_argument("--expected-unit", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    try:
        payload = validate(
            status_path=args.status.resolve(),
            target_genome=args.target_genome.resolve(),
            hy4a_genome=args.hy4a_genome.resolve(),
            hy4p_genome=args.hy4p_genome.resolve(),
            expected_unit=args.expected_unit,
        )
        atomic_json(args.output, payload, args.overwrite)
        print(f"PASS\t{args.expected_unit}\t{args.output}")
        return 0
    except (ValidationError, OSError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
