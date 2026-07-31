#!/usr/bin/env python3
"""Run the prepared CAFE5 base and gamma sensitivity models sequentially."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


class RunError(RuntimeError):
    pass


EXPECTED = {
    "base_poisson": (
        "Base_report.cafe", "Base_results.txt", "Base_family_likelihoods.txt", "Base_asr.tre",
        "Base_count.tab", "Base_change.tab", "Base_family_results.txt", "Base_clade_results.txt",
        "Base_branch_probabilities.tab",
    ),
    "gamma3_poisson": (
        "Gamma_report.cafe", "Gamma_results.txt", "Gamma_family_likelihoods.txt", "Gamma_asr.tre",
        "Gamma_count.tab", "Gamma_change.tab", "Gamma_family_results.txt", "Gamma_clade_results.txt",
        "Gamma_branch_probabilities.tab",
    ),
}


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path, *, allow_empty: bool = False) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink():
        raise RunError(f"missing or symlink file: {resolved}")
    if not allow_empty and resolved.stat().st_size == 0:
        raise RunError(f"empty file: {resolved}")
    return resolved


def write_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def validate_bundle(bundle: Path) -> tuple[dict[str, object], dict[str, str]]:
    if not bundle.is_dir() or bundle.is_symlink():
        raise RunError(f"invalid CAFE5 bundle directory: {bundle}")
    manifest_path = regular(bundle / "run_manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS_PREPARED_CAFE5" or manifest.get("workflow") != (
        "cafe5_timetree_secondary_run_bundle"
    ):
        raise RunError("CAFE5 bundle is not PASS_PREPARED_CAFE5")
    with regular(bundle / "checksums.tsv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["file", "sha256"]:
            raise RunError("invalid CAFE5 bundle checksum header")
        rows = list(reader)
    checksums: dict[str, str] = {}
    for row in rows:
        name, digest = row["file"], row["sha256"]
        if not name or Path(name).name != name or name in checksums or not re.fullmatch(
            r"[0-9a-f]{64}", digest
        ):
            raise RunError("invalid CAFE5 bundle checksum row")
        checksums[name] = digest
    inventory = {
        path.name for path in bundle.iterdir() if path.is_file() and path.name != "checksums.tsv"
    }
    if set(checksums) != inventory:
        raise RunError("CAFE5 bundle checksum inventory does not close")
    for name, digest in checksums.items():
        if sha256(regular(bundle / name)) != digest:
            raise RunError(f"CAFE5 bundle checksum mismatch: {name}")
    executable = regular(Path(str(manifest.get("cafe5", {}).get("path", ""))))
    if not os.access(executable, os.X_OK) or sha256(executable) != manifest["cafe5"]["sha256"]:
        raise RunError("CAFE5 executable binding changed")
    if manifest.get("cores_per_model") != 1 or manifest.get("models_run_sequentially") is not True:
        raise RunError("CAFE5 bundle violates the one-core sequential design")
    if manifest.get("large_families_in_rate_estimation") is not False:
        raise RunError("large families entered CAFE5 rate estimation")
    return manifest, checksums


def run_model(
    *, model: dict[str, object], bundle: Path, output: Path, executable: Path,
) -> dict[str, object]:
    model_id = str(model["model_id"])
    if model_id not in EXPECTED:
        raise RunError(f"unknown CAFE5 model: {model_id}")
    stage = output / "runs" / model_id
    stage.mkdir(parents=True)
    results = stage / "results"
    command = [
        str(executable), "--infile", str(bundle / "cafe5_primary_lt100.tsv"),
        "--tree", str(bundle / "dated_tree.mean_ma.tre"), "--cores", "1",
        *[str(value) for value in model["arguments"]], "--output_prefix", str(results),
    ]
    stdout = stage / "console.stdout"
    stderr = stage / "console.stderr"
    started = now()
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        completed = subprocess.run(command, cwd=stage, stdout=out, stderr=err, text=True, check=False)
    if completed.returncode != 0:
        raise RunError(f"{model_id}: CAFE5 exited {completed.returncode}")
    if not results.is_dir() or results.is_symlink():
        raise RunError(f"{model_id}: CAFE5 result directory is missing")
    for name in EXPECTED[model_id]:
        regular(results / name)
    inventory = sorted(path for path in results.rglob("*") if path.is_file())
    if not inventory or any(path.is_symlink() for path in inventory):
        raise RunError(f"{model_id}: invalid CAFE5 result inventory")
    return {
        "model_id": model_id,
        "started_at_utc": started,
        "finished_at_utc": now(),
        "returncode": completed.returncode,
        "command": command,
        "stdout": {"bytes": stdout.stat().st_size, "sha256": sha256(stdout)},
        "stderr": {"bytes": stderr.stat().st_size, "sha256": sha256(stderr)},
        "result_files": [
            {
                "relative_path": str(path.relative_to(results)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in inventory
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    output = args.output_dir.expanduser().resolve()
    state: dict[str, object] | None = None
    try:
        bundle = args.bundle.expanduser().resolve()
        manifest, checksums = validate_bundle(bundle)
        if output.exists():
            raise RunError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.mkdir()
        (output / "runs").mkdir()
        executable = regular(Path(str(manifest["cafe5"]["path"])))
        state = {
            "schema_version": 1,
            "workflow": "sequential_cafe5_timetree_secondary",
            "status": "running",
            "started_at_utc": now(),
            "calibration_claim": manifest["calibration_claim"],
            "bundle": {
                "path": str(bundle),
                "manifest_sha256": checksums["run_manifest.json"],
                "checksums_sha256": sha256(bundle / "checksums.tsv"),
            },
            "cafe5": manifest["cafe5"],
            "cores": 1,
            "large_family_count_excluded_from_rate_estimation": manifest["large_family_count"],
            "runtime_outlier_count_excluded_from_rate_estimation": manifest.get(
                "runtime_outlier_count", 0
            ),
            "completed": [],
        }
        state_path = output / "state.json"
        write_json(state_path, state)
        for model in manifest["models"]:
            state["active_model"] = model["model_id"]
            write_json(state_path, state)
            row = run_model(
                model=model, bundle=bundle, output=output, executable=executable
            )
            state["completed"].append(row)
            if validate_bundle(bundle)[1] != checksums:
                raise RunError("CAFE5 bundle changed during execution")
            if sha256(executable) != manifest["cafe5"]["sha256"]:
                raise RunError("CAFE5 executable changed during execution")
            write_json(state_path, state)
        state.pop("active_model", None)
        state["status"] = "PASS_RUN_COMPLETE"
        state["finished_at_utc"] = now()
        write_json(state_path, state)
        print(json.dumps(state, indent=2, sort_keys=True))
        return 0
    except (OSError, UnicodeError, ValueError, KeyError, json.JSONDecodeError, RunError) as error:
        if state is not None:
            state["status"] = "ERROR"
            state["error"] = str(error)
            state["finished_at_utc"] = now()
            write_json(output / "state.json", state)
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
