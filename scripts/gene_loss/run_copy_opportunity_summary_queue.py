#!/usr/bin/env python3
"""Wait for primary aggregation PASS, then summarize callable opportunities."""

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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-state", required=True, type=Path)
    parser.add_argument("--aggregation-dir", required=True, type=Path)
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
        for path in (args.python, args.worker):
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
            raise QueueError("another copy-opportunity controller owns this queue") from error
        state_path = queue / "state.json"
        state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "callable_copy_opportunity_summary_queue",
            "status": "running",
            "started_at_utc": now(),
        }
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
        state["status"] = "summarizing_callable_copy_opportunities"
        write_json(state_path, state)
        completed = subprocess.run(
            [
                str(args.python), str(args.worker),
                "--unit-calls", str(args.aggregation_dir / "unit_calls_long.tsv"),
                "--species-matrix", str(args.aggregation_dir / "species_gene_matrix.tsv"),
                "--shared-genes", str(args.aggregation_dir / "shared_positive_complete_genes.tsv"),
                "--output-dir", str(args.output_dir),
            ],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        (queue / "worker.stdout.log").write_text(completed.stdout, encoding="utf-8")
        (queue / "worker.stderr.log").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode:
            raise QueueError(f"summary worker exited {completed.returncode}")
        report = read_json(args.output_dir / "run_manifest.json")
        upstream = read_json(args.upstream_state)
        if (
            report.get("status") != "PASS"
            or report.get("workflow") != "callable_copy_opportunity_and_loss_mode_summary"
            or report.get("assembly_unit_count") != upstream.get("assembly_unit_count")
            or report.get("biological_species_count") != upstream.get("biological_species_count")
            or report.get("reference_gene_count") != upstream.get("reference_gene_count")
        ):
            raise QueueError("copy-opportunity summary did not close to upstream integration")
        state.update(
            {
                "status": "PASS", "finished_at_utc": now(),
                "assembly_unit_count": report.get("assembly_unit_count"),
                "biological_species_count": report.get("biological_species_count"),
                "reference_gene_count": report.get("reference_gene_count"),
                "shared_positive_complete_gene_count": report.get("shared_positive_complete_gene_count"),
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
