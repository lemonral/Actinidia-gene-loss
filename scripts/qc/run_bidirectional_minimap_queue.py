#!/usr/bin/env python3
"""Run validated chromosome-minimap jobs in bounded server-side lanes.

The controller first waits for every ``prerequisite`` row to finish and passes
each completed four-direction bundle through the immutable validator.  It then
runs ``queued`` rows sequentially within each declared lane.  Different lanes
run concurrently, so a two-lane manifest uses at most two minimap2 worker
pools.  Existing outputs, duplicate units, unsafe paths, stale prerequisites,
and a second controller are rejected rather than overwritten.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


HEADER = ("stage", "lane", "unit", "target_genome", "output_dir", "controller_pid")
STAGES = {"prerequisite", "queued"}


class QueueError(RuntimeError):
    """Raised when the queue cannot proceed without weakening a gate."""


@dataclass(frozen=True)
class Entry:
    stage: str
    lane: str
    unit: str
    target_genome: Path
    output_dir: Path
    controller_pid: Path | None


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
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise QueueError(f"{label} is missing or empty: {resolved}")
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


def read_manifest(path: Path, data_root: Path) -> list[Entry]:
    manifest = require_file(path, "queue manifest")
    try:
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != HEADER:
                raise QueueError(
                    f"manifest header must be exactly {HEADER}; found {reader.fieldnames}"
                )
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise QueueError(f"cannot read queue manifest: {error}") from error
    if not rows:
        raise QueueError("queue manifest has no entries")

    entries: list[Entry] = []
    units: set[str] = set()
    outputs: set[Path] = set()
    for line_number, row in enumerate(rows, start=2):
        stage = row["stage"].strip()
        lane = row["lane"].strip()
        unit = row["unit"].strip()
        if stage not in STAGES:
            raise QueueError(f"line {line_number}: invalid stage {stage!r}")
        if not unit or any(character.isspace() for character in unit):
            raise QueueError(f"line {line_number}: unsafe or empty unit {unit!r}")
        if unit in units:
            raise QueueError(f"line {line_number}: duplicate unit {unit!r}")
        units.add(unit)
        if stage == "queued" and not lane:
            raise QueueError(f"line {line_number}: queued entry requires a lane")
        if stage == "prerequisite" and lane:
            raise QueueError(f"line {line_number}: prerequisite lane must be empty")

        target_rel = safe_relative(row["target_genome"], f"line {line_number} target_genome")
        output_rel = safe_relative(row["output_dir"], f"line {line_number} output_dir")
        target = resolve_beneath(data_root, target_rel, "target genome")
        output = resolve_beneath(data_root, output_rel, "output directory")
        if output in outputs:
            raise QueueError(f"line {line_number}: duplicate output directory {output_rel}")
        outputs.add(output)

        pid_value = row["controller_pid"].strip()
        if stage == "prerequisite":
            pid_rel = safe_relative(pid_value, f"line {line_number} controller_pid")
            pid_path: Path | None = resolve_beneath(data_root, pid_rel, "controller PID")
        elif pid_value:
            raise QueueError(f"line {line_number}: queued controller_pid must be empty")
        else:
            pid_path = None
        entries.append(Entry(stage, lane, unit, target, output, pid_path))
    return entries


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def read_json(path: Path, label: str) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise QueueError(f"cannot read {label} {path}: {error}") from error
    if not isinstance(payload, dict):
        raise QueueError(f"{label} is not a JSON object: {path}")
    return payload


def pid_is_alive(path: Path) -> bool:
    try:
        value = int(path.read_text(encoding="utf-8").strip())
        os.kill(value, 0)
    except (OSError, ValueError):
        return False
    return True


def run_validator(
    *, entry: Entry, python: Path, validator: Path, hy4a: Path, hy4p: Path
) -> dict[str, object]:
    status = entry.output_dir / "status.json"
    output = entry.output_dir / "bundle_validation.json"
    if output.exists():
        raise QueueError(f"refusing to overwrite validation output: {output}")
    command = [
        str(python), str(validator),
        "--status", str(status),
        "--target-genome", str(entry.target_genome),
        "--hy4a-genome", str(hy4a),
        "--hy4p-genome", str(hy4p),
        "--expected-unit", entry.unit,
        "--output", str(output),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise QueueError(
            f"validator failed for {entry.unit} with exit {completed.returncode}: "
            f"{completed.stderr.strip()}"
        )
    payload = read_json(output, f"validation output for {entry.unit}")
    if payload.get("status") != "PASS" or payload.get("unit") != entry.unit:
        raise QueueError(f"validator did not publish exact PASS for {entry.unit}")
    return {
        "unit": entry.unit,
        "validation_sha256": sha256(output),
        "validated_at_utc": utc_now(),
    }


class Controller:
    def __init__(self, args: argparse.Namespace, entries: list[Entry]) -> None:
        self.args = args
        self.entries = entries
        self.state_lock = threading.Lock()
        self.state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "bidirectional_chromosome_minimap_queue",
            "status": "starting",
            "started_at_utc": utc_now(),
            "manifest_sha256": sha256(args.manifest),
            "max_lanes": args.max_lanes,
            "validated_prerequisites": [],
            "completed": [],
        }

    def publish(self, **updates: object) -> None:
        with self.state_lock:
            self.state.update(updates)
            atomic_json(self.args.state, self.state)

    def append(self, field: str, item: dict[str, object]) -> None:
        with self.state_lock:
            values = self.state.setdefault(field, [])
            if not isinstance(values, list):
                raise QueueError(f"internal state field is not a list: {field}")
            values.append(item)
            atomic_json(self.args.state, self.state)

    def wait_and_validate_prerequisites(self) -> None:
        prerequisites = [entry for entry in self.entries if entry.stage == "prerequisite"]
        pending = {entry.unit: entry for entry in prerequisites}
        self.publish(status="waiting_for_prerequisites", pending_prerequisites=sorted(pending))
        while pending:
            for unit, entry in list(pending.items()):
                status_path = entry.output_dir / "status.json"
                status = read_json(status_path, f"status for prerequisite {unit}")
                run_status = status.get("status")
                if run_status == "completed":
                    record = run_validator(
                        entry=entry,
                        python=self.args.python,
                        validator=self.args.validator,
                        hy4a=self.args.hy4a_genome,
                        hy4p=self.args.hy4p_genome,
                    )
                    self.append("validated_prerequisites", record)
                    del pending[unit]
                    self.publish(pending_prerequisites=sorted(pending))
                elif run_status in {"failed", "interrupted"}:
                    raise QueueError(f"prerequisite {unit} ended with status {run_status!r}")
                elif run_status != "running":
                    raise QueueError(f"prerequisite {unit} has unexpected status {run_status!r}")
                elif entry.controller_pid is None or not pid_is_alive(entry.controller_pid):
                    raise QueueError(f"prerequisite {unit} is running but its controller is not alive")
            if pending:
                time.sleep(self.args.poll_seconds)

    def run_entry(self, entry: Entry) -> dict[str, object]:
        if entry.output_dir.exists():
            raise QueueError(f"refusing existing queued output: {entry.output_dir}")
        stdout_log = entry.output_dir.parent / f"{entry.unit}.controller.stdout.log"
        stderr_log = entry.output_dir.parent / f"{entry.unit}.controller.stderr.log"
        for path in (stdout_log, stderr_log):
            if path.exists():
                raise QueueError(f"refusing existing queued controller log: {path}")
        command = [
            str(self.args.python), str(self.args.runner),
            "--unit", entry.unit,
            "--minimap2", str(self.args.minimap2),
            "--target-genome", str(entry.target_genome),
            "--hy4a-genome", str(self.args.hy4a_genome),
            "--hy4p-genome", str(self.args.hy4p_genome),
            "--output-dir", str(entry.output_dir),
        ]
        started = utc_now()
        stdout_log.parent.mkdir(parents=True, exist_ok=True)
        with stdout_log.open("wb") as stdout, stderr_log.open("wb") as stderr:
            completed = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
        if completed.returncode != 0:
            raise QueueError(f"runner failed for {entry.unit} with exit {completed.returncode}")
        validation = run_validator(
            entry=entry,
            python=self.args.python,
            validator=self.args.validator,
            hy4a=self.args.hy4a_genome,
            hy4p=self.args.hy4p_genome,
        )
        return {
            **validation,
            "lane": entry.lane,
            "started_at_utc": started,
            "finished_at_utc": utc_now(),
        }

    def run_lane(self, lane: str, entries: list[Entry]) -> None:
        for entry in entries:
            self.publish(status="running_queued_lanes")
            record = self.run_entry(entry)
            self.append("completed", record)

    def run(self) -> None:
        self.wait_and_validate_prerequisites()
        queued = [entry for entry in self.entries if entry.stage == "queued"]
        lanes: dict[str, list[Entry]] = {}
        for entry in queued:
            lanes.setdefault(entry.lane, []).append(entry)
        if not lanes:
            raise QueueError("manifest has no queued entries")
        if len(lanes) > self.args.max_lanes:
            raise QueueError(
                f"manifest requests {len(lanes)} lanes, exceeding --max-lanes {self.args.max_lanes}"
            )
        self.publish(status="running_queued_lanes", lanes=sorted(lanes))
        errors: list[str] = []
        with ThreadPoolExecutor(max_workers=len(lanes)) as executor:
            futures = {
                executor.submit(self.run_lane, lane, entries): lane
                for lane, entries in sorted(lanes.items())
            }
            for future in as_completed(futures):
                lane = futures[future]
                try:
                    future.result()
                except Exception as error:  # keep the other disjoint lane auditable
                    errors.append(f"{lane}: {error}")
        if errors:
            raise QueueError("; ".join(errors))
        completed = self.state.get("completed")
        if not isinstance(completed, list) or len(completed) != len(queued):
            raise QueueError("completed row count does not close against queued manifest")
        self.publish(status="completed", finished_at_utc=utc_now())


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--runner", required=True, type=Path)
    parser.add_argument("--validator", required=True, type=Path)
    parser.add_argument("--minimap2", required=True, type=Path)
    parser.add_argument("--hy4a-genome", required=True, type=Path)
    parser.add_argument("--hy4p-genome", required=True, type=Path)
    parser.add_argument("--state", required=True, type=Path)
    parser.add_argument("--lock", required=True, type=Path)
    parser.add_argument("--max-lanes", type=int, default=2)
    parser.add_argument("--poll-seconds", type=int, default=30)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.max_lanes < 1 or args.max_lanes > 2:
            raise QueueError("--max-lanes must be 1 or 2")
        if args.poll_seconds < 10 or args.poll_seconds > 60:
            raise QueueError("--poll-seconds must be between 10 and 60")
        args.data_root = args.data_root.expanduser().resolve()
        args.manifest = require_file(args.manifest, "queue manifest")
        args.python = require_file(args.python, "Python executable")
        args.runner = require_file(args.runner, "minimap runner")
        args.validator = require_file(args.validator, "bundle validator")
        args.minimap2 = require_file(args.minimap2, "minimap2 executable")
        args.hy4a_genome = require_file(args.hy4a_genome, "HY4A genome")
        args.hy4p_genome = require_file(args.hy4p_genome, "HY4P genome")
        args.state = args.state.expanduser().resolve()
        args.lock = args.lock.expanduser().resolve()
        if args.state.exists():
            raise QueueError(f"refusing existing state file: {args.state}")
        args.lock.parent.mkdir(parents=True, exist_ok=True)
        descriptor = os.open(args.lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(f"pid\t{os.getpid()}\nmanifest_sha256\t{sha256(args.manifest)}\n")
        entries = read_manifest(args.manifest, args.data_root)
        controller = Controller(args, entries)
        try:
            controller.run()
        except BaseException as error:
            controller.publish(status="failed", finished_at_utc=utc_now(), error=str(error))
            raise
        print(f"Completed chromosome minimap queue: {args.state}")
        return 0
    except (OSError, subprocess.SubprocessError, QueueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
