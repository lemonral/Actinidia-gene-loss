#!/usr/bin/env python3
"""Wait for new-unit matrices, then merge and aggregate the 23 primary units."""

from __future__ import annotations

import argparse
import fcntl
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


def read_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise QueueError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise QueueError(f"{path}: JSON root is not an object")
    return value


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def run(command: list[str], stdout: Path, stderr: Path) -> None:
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    stdout.write_text(completed.stdout, encoding="utf-8")
    stderr.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise QueueError(f"worker failed with exit code {completed.returncode}: {command[1]}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-state", required=True, type=Path)
    parser.add_argument("--sources", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--expected-total-units", required=True, type=int)
    parser.add_argument("--merged-output", required=True, type=Path)
    parser.add_argument("--aggregation-metadata", required=True, type=Path)
    parser.add_argument("--aggregation-output", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--merge-worker", required=True, type=Path)
    parser.add_argument("--aggregation-worker", required=True, type=Path)
    parser.add_argument("--queue-root", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()
    lock_handle = None
    try:
        if args.poll_seconds < 1 or args.expected_total_units < 1:
            raise QueueError("poll-seconds and expected-total-units must be positive")
        for path in (
            args.sources, args.aggregation_metadata, args.python,
            args.merge_worker, args.aggregation_worker,
        ):
            if not path.is_file() or path.stat().st_size == 0:
                raise QueueError(f"missing input/executable: {path}")
        if not os.access(args.python, os.X_OK):
            raise QueueError(f"python is not executable: {args.python}")
        queue = args.queue_root.resolve()
        queue.mkdir(parents=True, exist_ok=True)
        lock_handle = (queue / "controller.lock").open("w")
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise QueueError("another primary loss integration controller owns this queue") from error
        state_path = queue / "state.json"
        state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "primary_complete_loss_integration_queue",
            "status": "running",
            "started_at_utc": now(),
        }
        write_json(state_path, state)
        while True:
            if args.upstream_state.is_file():
                upstream = read_json(args.upstream_state)
                if upstream.get("workflow") != "postsearch_complete_loss_matrix_and_aggregation_queue":
                    raise QueueError("unexpected upstream matrix workflow")
                if upstream.get("status") == "PASS":
                    break
                if upstream.get("status") not in {
                    "running", "building_complete_matrix", "aggregating_species_and_subgenomes"
                }:
                    raise QueueError(f"upstream terminal status is {upstream.get('status')!r}")
            state["waiting_for"] = str(args.upstream_state)
            write_json(state_path, state)
            time.sleep(args.poll_seconds)

        state.pop("waiting_for", None)
        state["status"] = "merging_primary_complete_matrices"
        write_json(state_path, state)
        run(
            [
                str(args.python), str(args.merge_worker),
                "--sources", str(args.sources),
                "--data-root", str(args.data_root),
                "--expected-total-units", str(args.expected_total_units),
                "--output-dir", str(args.merged_output),
            ],
            queue / "merge.stdout.log",
            queue / "merge.stderr.log",
        )
        merged = read_json(args.merged_output / "run_manifest.json")
        if (
            merged.get("status") != "PASS"
            or merged.get("workflow") != "merged_primary_complete_loss_matrix"
            or merged.get("assembly_unit_count") != args.expected_total_units
            or merged.get("matrix_rows") != merged.get("expected_matrix_rows")
        ):
            raise QueueError("merged primary matrix did not pass exact closure")

        state["status"] = "aggregating_primary_species_and_subgenomes"
        write_json(state_path, state)
        run(
            [
                str(args.python), str(args.aggregation_worker),
                "--unit-call-matrix", str(args.merged_output / "complete_unit_loss_matrix.tsv"),
                "--unit-metadata", str(args.aggregation_metadata),
                "--output-dir", str(args.aggregation_output),
            ],
            queue / "aggregation.stdout.log",
            queue / "aggregation.stderr.log",
        )
        aggregation = read_json(args.aggregation_output / "species_loss_summary.json")
        if (
            aggregation.get("status") != "PASS"
            or aggregation.get("assembly_unit_count") != args.expected_total_units
            or aggregation.get("reference_gene_count") != merged.get("reference_gene_count")
        ):
            raise QueueError("primary species aggregation did not close to merged matrix")
        state.update(
            {
                "status": "PASS",
                "finished_at_utc": now(),
                "assembly_unit_count": args.expected_total_units,
                "reference_gene_count": merged.get("reference_gene_count"),
                "matrix_rows": merged.get("matrix_rows"),
                "biological_species_count": aggregation.get("biological_species_count"),
                "shared_positive_complete_gene_count": aggregation.get(
                    "shared_positive_complete_gene_count"
                ),
            }
        )
        write_json(state_path, state)
        print(f"PASS\t{args.expected_total_units} units\t{state['biological_species_count']} lineages")
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, QueueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if lock_handle is not None:
            lock_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
