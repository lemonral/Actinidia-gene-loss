#!/usr/bin/env python3
"""Wait for species-loss aggregation PASS, then build and run exact matching-tree PGLS."""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


class QueueError(RuntimeError):
    pass


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def regular(path: Path, *, executable: bool = False) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink() or resolved.stat().st_size == 0:
        raise QueueError(f"missing, empty, or symlink input: {resolved}")
    if executable and not os.access(resolved, os.X_OK):
        raise QueueError(f"not executable: {resolved}")
    return resolved


def read_json(path: Path) -> dict[str, object]:
    value = json.loads(regular(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise QueueError(f"JSON object required: {path}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def run_worker(command: list[str], stdout_path: Path, stderr_path: Path) -> None:
    completed = subprocess.run(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False
    )
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise QueueError(f"worker exited {completed.returncode}: {Path(command[1]).name}")


def validate_checksums(directory: Path) -> int:
    directory = directory.expanduser().resolve()
    checksum_path = regular(directory / "checksums.sha256.tsv")
    with checksum_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames not in (["file", "sha256"], ["relative_path", "sha256"]):
            raise QueueError("invalid PGLS checksum header")
        name_column = reader.fieldnames[0]
        rows = list(reader)
    observed: set[str] = set()
    for row in rows:
        name = row[name_column]
        candidate = Path(name)
        if (
            not name or candidate.is_absolute() or ".." in candidate.parts or name in observed
            or sha256(regular(directory / candidate)) != row["sha256"]
        ):
            raise QueueError("PGLS checksum closure failed")
        observed.add(name)
    inventory = {
        str(path.relative_to(directory)) for path in directory.rglob("*")
        if path.is_file() and path.resolve() != checksum_path
    }
    if observed != inventory:
        raise QueueError("PGLS checksum inventory differs from output files")
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-state", required=True, type=Path)
    parser.add_argument("--species-loss-dir", required=True, type=Path)
    parser.add_argument("--ploidy-ledger", required=True, type=Path)
    parser.add_argument("--time-tree-dir", required=True, type=Path)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--input-builder", required=True, type=Path)
    parser.add_argument("--pgls-worker", required=True, type=Path)
    parser.add_argument("--queue-root", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()
    lock_handle = None
    try:
        if args.poll_seconds < 1:
            raise QueueError("poll-seconds must be positive")
        python = regular(args.python, executable=True)
        input_builder = regular(args.input_builder)
        pgls_worker = regular(args.pgls_worker)
        ploidy = regular(args.ploidy_ledger)
        tree = regular(args.time_tree_dir / "species_time_tree.nwk")
        tree_pass = regular(args.time_tree_dir / "species_time_tree_pass.json")
        tree_report = read_json(tree_pass)
        if (
            tree_report.get("status") != "PASS"
            or tree_report.get("workflow") != "species_time_tree_validation"
            or not isinstance(tree_report.get("biological_species"), list)
        ):
            raise QueueError("PGLS time-tree bundle is not PASS")
        if os.path.lexists(args.input_dir) or os.path.lexists(args.output_dir):
            raise QueueError("PGLS output already exists")
        queue = args.queue_root.expanduser().resolve()
        queue.mkdir(parents=True, exist_ok=True)
        lock_handle = (queue / "controller.lock").open("w")
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise QueueError("another species-PGLS controller owns this queue") from error
        state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "species_pgls_after_primary_loss_queue",
            "status": "running",
            "started_at_utc": now(),
            "time_tree_pass_sha256": sha256(tree_pass),
            "ploidy_ledger_sha256": sha256(ploidy),
        }
        state_path = queue / "state.json"
        write_json(state_path, state)
        while True:
            if args.upstream_state.is_file():
                upstream = read_json(args.upstream_state)
                if upstream.get("workflow") != "primary_complete_loss_integration_queue":
                    raise QueueError("unexpected upstream integration workflow")
                if upstream.get("status") == "PASS":
                    break
                if upstream.get("status") not in {
                    "running", "merging_primary_complete_matrices",
                    "aggregating_primary_species_and_subgenomes",
                }:
                    raise QueueError(f"upstream terminal status is {upstream.get('status')!r}")
            state["waiting_for"] = str(args.upstream_state)
            write_json(state_path, state)
            time.sleep(args.poll_seconds)
        state.pop("waiting_for", None)
        state["status"] = "building_pgls_input"
        write_json(state_path, state)
        run_worker(
            [
                str(python), str(input_builder),
                "--species-loss-dir", str(args.species_loss_dir),
                "--ploidy-ledger", str(ploidy),
                "--output-dir", str(args.input_dir),
            ],
            queue / "input_builder.stdout.log",
            queue / "input_builder.stderr.log",
        )
        input_report = read_json(args.input_dir / "pgls_input_pass.json")
        if (
            input_report.get("status") != "PASS"
            or input_report.get("workflow") != "species_pgls_input_builder"
            or input_report.get("biological_species") != tree_report["biological_species"]
        ):
            raise QueueError("PGLS input does not close to the exact time-tree species")
        validate_checksums(args.input_dir)
        state["status"] = "running_species_pgls"
        write_json(state_path, state)
        run_worker(
            [
                str(python), str(pgls_worker),
                "--data", str(args.input_dir / "pgls_input.tsv"),
                "--time-tree", str(tree),
                "--input-pass-report", str(args.input_dir / "pgls_input_pass.json"),
                "--species-loss-manifest", str(args.species_loss_dir / "species_loss_summary.json"),
                "--ploidy-ledger-pass-report", str(args.input_dir / "ploidy_ledger_pass.json"),
                "--time-tree-pass-report", str(tree_pass),
                "--predictor-column", "log2_ploidy",
                "--sensitivity", "without_rufa=Actinidia rufa",
                "--output-dir", str(args.output_dir),
            ],
            queue / "pgls.stdout.log",
            queue / "pgls.stderr.log",
        )
        manifest = read_json(args.output_dir / "analysis_manifest.json")
        checksum_count = validate_checksums(args.output_dir)
        if (
            manifest.get("status") != "COMPLETE_EXPLORATORY_BLOCKED_FOR_PUBLICATION"
            or manifest.get("publication_gate") != "BLOCKED_DENOMINATOR_AWARE_MODEL_REQUIRED"
            or manifest.get("species_count") != len(tree_report["biological_species"])
            or manifest.get("named_exclusion_sensitivities")
            != {"without_rufa": ["Actinidia rufa"]}
        ):
            raise QueueError("species PGLS output does not match the frozen exploratory design")
        state.update({
            "status": "PASS_EXPLORATORY_PGLS",
            "finished_at_utc": now(),
            "species_count": manifest["species_count"],
            "publication_gate": manifest["publication_gate"],
            "output_checksum_count": checksum_count,
            "analysis_manifest_sha256": sha256(args.output_dir / "analysis_manifest.json"),
        })
        write_json(state_path, state)
        print("PASS_EXPLORATORY_PGLS")
        return 0
    except (OSError, UnicodeError, ValueError, KeyError, TypeError, json.JSONDecodeError, QueueError) as error:
        if args.queue_root.exists():
            try:
                current = read_json(args.queue_root / "state.json")
            except Exception:
                current = {"schema_version": 1, "workflow": "species_pgls_after_primary_loss_queue"}
            current.update({"status": "ERROR", "finished_at_utc": now(), "error": str(error)})
            write_json(args.queue_root / "state.json", current)
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if lock_handle is not None:
            lock_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
