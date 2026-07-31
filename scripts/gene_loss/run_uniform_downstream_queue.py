#!/usr/bin/env python3
"""Wait for the uniform search and build matrix, spatial, aggregation and PGLS outputs."""

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


def strict_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise QueueError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise QueueError(f"{path}: JSON root is not an object")
    return value


def write_state(path: Path, value: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def run_stage(name: str, command: list[str], queue: Path, state: dict[str, object], cwd: Path) -> None:
    state["status"] = "running"
    state["active_stage"] = name
    write_state(queue / "state.json", state)
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    with (queue / f"{name}.stdout.log").open("wb") as stdout, (queue / f"{name}.stderr.log").open("wb") as stderr:
        completed = subprocess.run(
            command, cwd=cwd, env=environment, stdout=stdout, stderr=stderr, check=False
        )
    if completed.returncode:
        raise QueueError(f"{name} failed with exit code {completed.returncode}")
    state.setdefault("completed_stages", []).append(name)
    state.pop("active_stage", None)
    write_state(queue / "state.json", state)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--upstream-state", required=True, type=Path)
    parser.add_argument("--queue-root", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=30)
    args = parser.parse_args()
    lock_handle = None
    state: dict[str, object] = {}
    try:
        if args.poll_seconds < 1:
            raise QueueError("poll interval must be positive")
        repository = args.repository.resolve()
        data = args.data_root.resolve()
        queue = args.queue_root.resolve()
        queue.mkdir(parents=True, exist_ok=True)
        lock_handle = (queue / "controller.lock").open("w")
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise QueueError("another downstream controller owns this queue") from error
        state = {
            "schema_version": 1,
            "workflow": "uniform_loss_matrix_spatial_species_pgls_queue",
            "status": "waiting_for_uniform_search",
            "started_at_utc": utc_now(),
            "completed_stages": [],
        }
        write_state(queue / "state.json", state)
        while True:
            upstream = strict_json(args.upstream_state)
            status = upstream.get("status")
            if status == "PASS":
                break
            if status not in {"running"}:
                raise QueueError(f"upstream uniform search ended with status {status!r}")
            time.sleep(args.poll_seconds)

        matrix = data / "qc/primary_uniform_complete_loss_matrix_v1_20260722"
        species = data / "qc/primary_uniform_species_loss_aggregation_v1_20260722"
        pgls_input = data / "qc/primary_uniform_species_pgls_input_v1_20260722"
        pgls_output = data / "results/statistics/species_pgls_uniform_v2_20260722"
        spatial = data / "results/statistics/uniform_loss_positions_v1_20260722"
        protected = [matrix, species, pgls_input, pgls_output, spatial]
        existing = [str(path) for path in protected if path.exists() or path.is_symlink()]
        if existing:
            raise QueueError("refusing to reuse/overwrite downstream outputs: " + ", ".join(existing))

        python = sys.executable
        run_stage(
            "complete_matrix",
            [python, "scripts/gene_loss/build_uniform_complete_loss_matrix.py", "--config", "config/primary_uniform_loss_matrix_v1.tsv", "--data-root", str(data), "--reference-protein", str(data / "legacy_linked/clem_scandens_reference_legacy/protein.faa"), "--output-dir", str(matrix)],
            queue, state, repository,
        )
        run_stage(
            "species_aggregation",
            [python, "scripts/gene_loss/aggregate_species_loss.py", "--unit-call-matrix", str(matrix / "complete_unit_loss_matrix.tsv"), "--unit-metadata", "config/primary_species_loss_aggregation.tsv", "--output-dir", str(species)],
            queue, state, repository,
        )
        run_stage(
            "pgls_input",
            [python, "scripts/downstream/build_species_pgls_input.py", "--species-loss-dir", str(species), "--ploidy-ledger", "config/species_ploidy.tsv", "--output-dir", str(pgls_input)],
            queue, state, repository,
        )
        tree = data / "phylogeny/species_pgls_time_tree_v1_20260720"
        run_stage(
            "pgls",
            [python, "scripts/downstream/species_pgls.py", "--data", str(pgls_input / "pgls_input.tsv"), "--time-tree", str(tree / "species_time_tree.nwk"), "--input-pass-report", str(pgls_input / "pgls_input_pass.json"), "--species-loss-manifest", str(species / "species_loss_summary.json"), "--ploidy-ledger-pass-report", str(pgls_input / "ploidy_ledger_pass.json"), "--time-tree-pass-report", str(tree / "species_time_tree_pass.json"), "--predictor-column", "log2_ploidy", "--sensitivity", "without_rufa=Actinidia rufa", "--output-dir", str(pgls_output)],
            queue, state, repository,
        )
        run_stage(
            "spatial_model",
            [python, "scripts/spatial/analyze_uniform_loss_positions.py", "--matrix-dir", str(matrix), "--uniform-config", "config/primary_uniform_miniprot_loss_v1.tsv", "--data-root", str(data), "--output-dir", str(spatial)],
            queue, state, repository,
        )
        state.update(
            {
                "status": "PASS",
                "finished_at_utc": utc_now(),
                "outputs": {
                    "complete_matrix": str(matrix.relative_to(data)),
                    "spatial_model": str(spatial.relative_to(data)),
                    "species_aggregation": str(species.relative_to(data)),
                    "pgls": str(pgls_output.relative_to(data)),
                },
            }
        )
        write_state(queue / "state.json", state)
        print("PASS")
        return 0
    except (OSError, QueueError) as error:
        state.update({"status": "ERROR", "finished_at_utc": utc_now(), "error": str(error)})
        try:
            write_state(args.queue_root.resolve() / "state.json", state)
        except OSError:
            pass
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if lock_handle is not None:
            lock_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
