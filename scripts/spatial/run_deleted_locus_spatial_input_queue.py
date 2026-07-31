#!/usr/bin/env python3
"""Wait for exact search/relabel PASS states, then build deleted-locus spatial inputs."""

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


def utc_now() -> str:
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


def wait_for_pass(path: Path, expected_workflow: str, state: dict[str, object], state_path: Path, poll: int) -> dict[str, object]:
    while True:
        if path.is_file():
            payload = read_json(path)
            status = payload.get("status")
            workflow = payload.get("workflow")
            if workflow != expected_workflow:
                raise QueueError(f"{path}: unexpected workflow {workflow!r}")
            if status == "PASS":
                return payload
            if status not in {"running", "RUNNING"}:
                raise QueueError(f"{path}: prerequisite entered non-PASS terminal status {status!r}")
        state["waiting_for"] = str(path)
        write_json(state_path, state)
        time.sleep(poll)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search-queue-state", required=True, type=Path)
    parser.add_argument("--relabel-queue-state", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--worker", required=True, type=Path)
    parser.add_argument("--queue-root", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()
    lock_handle = None
    try:
        if args.poll_seconds < 1:
            raise QueueError("poll-seconds must be positive")
        for path in (args.manifest, args.python, args.worker):
            if not path.is_file() or path.stat().st_size == 0:
                raise QueueError(f"missing input/executable: {path}")
        if not os.access(args.python, os.X_OK):
            raise QueueError(f"python is not executable: {args.python}")
        queue_root = args.queue_root.resolve()
        queue_root.mkdir(parents=True, exist_ok=True)
        lock_handle = (queue_root / "controller.lock").open("w")
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise QueueError("another deleted-locus input controller owns this queue") from error
        state_path = queue_root / "state.json"
        state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "deleted_locus_spatial_input_queue",
            "status": "running",
            "started_at_utc": utc_now(),
        }
        write_json(state_path, state)
        search = wait_for_pass(
            args.search_queue_state,
            "callable_aware_translated_search_queue",
            state,
            state_path,
            args.poll_seconds,
        )
        relabel = wait_for_pass(
            args.relabel_queue_state,
            "sequential_hy4a_similarity_relabelling_queue",
            state,
            state_path,
            args.poll_seconds,
        )
        state.pop("waiting_for", None)
        state["status"] = "building_inputs"
        write_json(state_path, state)
        completed = subprocess.run(
            [
                str(args.python),
                str(args.worker),
                "--manifest",
                str(args.manifest),
                "--data-root",
                str(args.data_root),
                "--output-dir",
                str(args.output_dir),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        (queue_root / "worker.stdout.log").write_text(completed.stdout, encoding="utf-8")
        (queue_root / "worker.stderr.log").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode:
            raise QueueError(f"spatial input worker failed with exit code {completed.returncode}")
        report = read_json(args.output_dir / "run_manifest.json")
        if report.get("status") != "PASS" or report.get("workflow") != "callable_positive_deleted_expected_locus_spatial_inputs":
            raise QueueError("spatial input worker did not produce the exact PASS bundle")
        state.update(
            {
                "status": "PASS",
                "finished_at_utc": utc_now(),
                "translated_search_completed_units": len(search.get("completed", [])),
                "relabel_completed_units": len(relabel.get("completed", [])),
                "spatial_input_units": report.get("unit_count"),
                "positive_deleted_count": report.get("positive_deleted_count"),
            }
        )
        write_json(state_path, state)
        print(f"PASS\t{report.get('unit_count')} units\t{report.get('positive_deleted_count')} loci")
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, QueueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if lock_handle is not None:
            lock_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
