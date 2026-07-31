#!/usr/bin/env python3
"""Run a one-worker chromosome relabelling queue after a prerequisite PID exits."""

from __future__ import annotations

import argparse
import csv
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


REQUIRED = ("unit", "label_map", "genome", "gff", "cds", "protein", "output_dir")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def write_state(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def resolve(root: Path, value: str, *, input_file: bool = False) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise QueueError(f"unsafe data-root-relative path: {value!r}")
    result = (root / relative).resolve()
    if not result.is_relative_to(root):
        raise QueueError(f"path escapes data root: {value!r}")
    if input_file and (not result.is_file() or result.stat().st_size <= 0):
        raise QueueError(f"missing input: {value!r}")
    return result


def live_pid(pid: int) -> bool:
    return pid > 1 and Path(f"/proc/{pid}").exists()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--worker", required=True, type=Path)
    parser.add_argument("--queue-root", required=True, type=Path)
    parser.add_argument("--wait-for-pid-file", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=30)
    args = parser.parse_args()
    lock_handle = None
    try:
        if args.poll_seconds < 1:
            raise QueueError("poll seconds must be positive")
        root = args.data_root.resolve()
        with args.manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != REQUIRED:
                raise QueueError("manifest columns differ from exact schema")
            rows = list(reader)
        if not rows or len({row["unit"] for row in rows}) != len(rows):
            raise QueueError("manifest is empty or has duplicate units")
        queue = args.queue_root.resolve()
        queue.mkdir(parents=True, exist_ok=True)
        lock_handle = (queue / "controller.lock").open("w")
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise QueueError("another relabelling controller owns the queue") from error
        state_path = queue / "state.json"
        state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "sequential_hy4a_similarity_relabelling_queue",
            "status": "waiting_for_prerequisite",
            "started_at_utc": utc_now(),
            "pending": [row["unit"] for row in rows],
            "completed": [],
        }
        write_state(state_path, state)
        try:
            prerequisite_pid = int(args.wait_for_pid_file.read_text(encoding="utf-8").strip())
        except (OSError, ValueError) as error:
            raise QueueError("cannot read prerequisite PID") from error
        while live_pid(prerequisite_pid):
            time.sleep(args.poll_seconds)
        state["status"] = "running"
        state["prerequisite_pid_exited"] = prerequisite_pid
        write_state(state_path, state)
        for row in rows:
            output = resolve(root, row["output_dir"])
            if output.exists():
                raise QueueError(f"{row['unit']}: refusing to overwrite output")
            state["active_unit"] = row["unit"]
            write_state(state_path, state)
            command = [
                str(args.python), str(args.worker), "--unit", row["unit"],
                "--label-map", str(resolve(root, row["label_map"], input_file=True)),
                "--genome", str(resolve(root, row["genome"], input_file=True)),
                "--gff", str(resolve(root, row["gff"], input_file=True)),
                "--cds", str(resolve(root, row["cds"], input_file=True)),
                "--protein", str(resolve(root, row["protein"], input_file=True)),
                "--output-dir", str(output),
            ]
            completed = subprocess.run(command, check=False)
            if completed.returncode:
                raise QueueError(f"{row['unit']}: relabelling worker failed")
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            if manifest.get("status") != "PASS" or manifest.get("unit") != row["unit"]:
                raise QueueError(f"{row['unit']}: output manifest is not PASS")
            state["completed"].append({  # type: ignore[union-attr]
                "unit": row["unit"], "finished_at_utc": utc_now(),
            })
            state["pending"] = [item for item in state["pending"] if item != row["unit"]]  # type: ignore[union-attr]
            state.pop("active_unit", None)
            write_state(state_path, state)
        state["status"] = "PASS"
        state["finished_at_utc"] = utc_now()
        write_state(state_path, state)
        return 0
    except (OSError, UnicodeError, csv.Error, json.JSONDecodeError, QueueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if lock_handle is not None:
            lock_handle.close()


if __name__ == "__main__":
    sys.exit(main())
