#!/usr/bin/env python3
"""Wait for translated-search PASS, then build and aggregate the complete loss matrix."""

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
    parser.add_argument("--search-queue-state", required=True, type=Path)
    parser.add_argument("--matrix-manifest", required=True, type=Path)
    parser.add_argument("--reference-cds", required=True)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--matrix-output", required=True, type=Path)
    parser.add_argument("--aggregation-metadata", required=True, type=Path)
    parser.add_argument("--aggregation-output", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--matrix-worker", required=True, type=Path)
    parser.add_argument("--aggregation-worker", required=True, type=Path)
    parser.add_argument("--queue-root", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()
    lock_handle = None
    try:
        if args.poll_seconds < 1:
            raise QueueError("poll-seconds must be positive")
        for path in (
            args.matrix_manifest,
            args.aggregation_metadata,
            args.python,
            args.matrix_worker,
            args.aggregation_worker,
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
            raise QueueError("another post-search loss controller owns this queue") from error
        state_path = queue / "state.json"
        state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "postsearch_complete_loss_matrix_and_aggregation_queue",
            "status": "running",
            "started_at_utc": now(),
        }
        write_json(state_path, state)
        while True:
            if args.search_queue_state.is_file():
                search = read_json(args.search_queue_state)
                if search.get("workflow") != "callable_aware_translated_search_queue":
                    raise QueueError("unexpected translated-search workflow")
                if search.get("status") == "PASS":
                    break
                if search.get("status") not in {"running", "RUNNING"}:
                    raise QueueError(f"translated-search terminal status is {search.get('status')!r}")
            state["waiting_for"] = str(args.search_queue_state)
            write_json(state_path, state)
            time.sleep(args.poll_seconds)

        state.pop("waiting_for", None)
        state["status"] = "building_complete_matrix"
        write_json(state_path, state)
        run(
            [
                str(args.python), str(args.matrix_worker),
                "--manifest", str(args.matrix_manifest),
                "--data-root", str(args.data_root),
                "--reference-cds", args.reference_cds,
                "--output-dir", str(args.matrix_output),
            ],
            queue / "matrix.stdout.log",
            queue / "matrix.stderr.log",
        )
        matrix_report = read_json(args.matrix_output / "run_manifest.json")
        if (
            matrix_report.get("status") != "PASS"
            or matrix_report.get("workflow") != "complete_callable_aware_new_unit_loss_matrix"
            or matrix_report.get("matrix_rows") != matrix_report.get("expected_matrix_rows")
        ):
            raise QueueError("complete matrix did not pass exact closure")

        state["status"] = "aggregating_species_and_subgenomes"
        write_json(state_path, state)
        run(
            [
                str(args.python), str(args.aggregation_worker),
                "--unit-call-matrix", str(args.matrix_output / "complete_unit_loss_matrix.tsv"),
                "--unit-metadata", str(args.aggregation_metadata),
                "--output-dir", str(args.aggregation_output),
            ],
            queue / "aggregation.stdout.log",
            queue / "aggregation.stderr.log",
        )
        aggregation = read_json(args.aggregation_output / "species_loss_summary.json")
        if (
            aggregation.get("status") != "PASS"
            or aggregation.get("assembly_unit_count") != matrix_report.get("assembly_unit_count")
            or aggregation.get("reference_gene_count") != matrix_report.get("reference_gene_count")
        ):
            raise QueueError("species/subgenome aggregation did not close to the complete matrix")
        state.update(
            {
                "status": "PASS",
                "finished_at_utc": now(),
                "assembly_unit_count": matrix_report.get("assembly_unit_count"),
                "reference_gene_count": matrix_report.get("reference_gene_count"),
                "matrix_rows": matrix_report.get("matrix_rows"),
                "biological_species_count": aggregation.get("biological_species_count"),
                "shared_positive_complete_gene_count": aggregation.get("shared_positive_complete_gene_count"),
            }
        )
        write_json(state_path, state)
        print(f"PASS\t{state['assembly_unit_count']} units\t{state['reference_gene_count']} genes")
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, QueueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if lock_handle is not None:
            lock_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
