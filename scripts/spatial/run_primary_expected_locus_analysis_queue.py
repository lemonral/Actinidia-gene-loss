#!/usr/bin/env python3
"""Wait for legacy/new spatial inputs, then merge and analyze end distance."""

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
    parser.add_argument("--new-spatial-state", required=True, type=Path)
    parser.add_argument("--legacy-input-dir", required=True, type=Path)
    parser.add_argument("--sources", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--expected-total-units", type=int, default=23)
    parser.add_argument("--merged-input-dir", required=True, type=Path)
    parser.add_argument("--analysis-output-dir", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--merge-worker", required=True, type=Path)
    parser.add_argument("--analysis-worker", required=True, type=Path)
    parser.add_argument("--queue-root", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()
    lock_handle = None
    try:
        if args.poll_seconds < 1 or args.expected_total_units < 1:
            raise QueueError("poll-seconds and expected-total-units must be positive")
        for path in (args.sources, args.python, args.merge_worker, args.analysis_worker):
            if not path.is_file() or path.stat().st_size == 0:
                raise QueueError(f"missing input/executable: {path}")
        queue = args.queue_root.resolve()
        queue.mkdir(parents=True, exist_ok=True)
        lock_handle = (queue / "controller.lock").open("w")
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise QueueError("another primary spatial controller owns this queue") from error
        state_path = queue / "state.json"
        state: dict[str, object] = {
            "schema_version": 1, "workflow": "primary_expected_locus_end_distance_queue",
            "status": "running", "started_at_utc": now(),
        }
        write_json(state_path, state)
        while True:
            new_ready = False
            legacy_ready = False
            if args.new_spatial_state.is_file():
                new_state = read_json(args.new_spatial_state)
                if new_state.get("workflow") != "deleted_locus_spatial_input_queue":
                    raise QueueError("unexpected new spatial queue workflow")
                new_ready = new_state.get("status") == "PASS"
                if not new_ready and new_state.get("status") not in {"running"}:
                    raise QueueError(f"new spatial terminal status is {new_state.get('status')!r}")
            legacy_report_path = args.legacy_input_dir / "run_manifest.json"
            if legacy_report_path.is_file():
                legacy = read_json(legacy_report_path)
                legacy_ready = (
                    legacy.get("status") == "PASS"
                    and legacy.get("workflow") == "legacy_conservative_deleted_expected_locus_spatial_inputs"
                )
                if not legacy_ready:
                    raise QueueError("legacy spatial input report is not exact PASS")
            if new_ready and legacy_ready:
                break
            state["waiting_for"] = [str(args.new_spatial_state), str(legacy_report_path)]
            write_json(state_path, state)
            time.sleep(args.poll_seconds)

        state.pop("waiting_for", None)
        state["status"] = "merging_spatial_inputs"
        write_json(state_path, state)
        run(
            [
                str(args.python), str(args.merge_worker), "--sources", str(args.sources),
                "--data-root", str(args.data_root), "--expected-total-units",
                str(args.expected_total_units), "--output-dir", str(args.merged_input_dir),
            ],
            queue / "merge.stdout.log", queue / "merge.stderr.log",
        )
        merged = read_json(args.merged_input_dir / "run_manifest.json")
        if (
            merged.get("status") != "PASS"
            or merged.get("workflow") != "merged_primary_expected_deleted_locus_spatial_inputs"
            or merged.get("unit_count") != args.expected_total_units
        ):
            raise QueueError("merged spatial inputs did not pass exact closure")

        state["status"] = "analyzing_expected_locus_end_distance"
        write_json(state_path, state)
        run(
            [
                str(args.python), str(args.analysis_worker),
                "--positive-calls", str(args.merged_input_dir / "positive_deleted_calls.tsv"),
                "--feature-coordinates", str(args.merged_input_dir / "expected_deleted_locus_coordinates.tsv"),
                "--assembly-manifest", str(args.merged_input_dir / "assembly_manifest.tsv"),
                "--output-dir", str(args.analysis_output_dir),
                "--analysis-label", "primary_callable_positive_deleted_expected_locus",
                "--positive-classes", "positive_deleted", "--number-of-bins", "5",
                "--legacy-reproduction",
                "--call-unit-column", "assembly_unit_id",
                "--call-gene-column", "reference_gene_id",
                "--call-classification-column", "classification",
                "--coordinate-unit-column", "assembly_unit_id",
                "--coordinate-gene-column", "reference_gene_id",
                "--coordinate-chromosome-column", "chromosome",
                "--coordinate-start-column", "expected_locus_start_1based",
                "--coordinate-end-column", "expected_locus_end_1based",
                "--coordinate-classification-column", "classification",
            ],
            queue / "analysis.stdout.log", queue / "analysis.stderr.log",
        )
        summary = read_json(args.analysis_output_dir / "run_summary.json")
        reconciliation = summary.get("reconciliation")
        if (
            not isinstance(reconciliation, dict)
            or reconciliation.get("assembly_unit_count") != args.expected_total_units
            or reconciliation.get("positive_call_count") != merged.get("positive_deleted_count")
            or reconciliation.get("emitted_position_count") != merged.get("positive_deleted_count")
        ):
            raise QueueError("spatial analysis did not close to merged inputs")
        state.update(
            {
                "status": "PASS", "finished_at_utc": now(),
                "assembly_unit_count": args.expected_total_units,
                "positive_deleted_count": merged.get("positive_deleted_count"),
                "centromere_policy": summary.get("centromere_policy"),
            }
        )
        write_json(state_path, state)
        print("PASS")
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, QueueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if lock_handle is not None:
            lock_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
