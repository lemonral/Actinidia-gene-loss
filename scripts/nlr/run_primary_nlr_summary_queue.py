#!/usr/bin/env python3
"""Wait for the primary NLR batch, then validate, summarize, and render it.

This controller owns one file lock and never starts NLR-Annotator.  It remains
idle until the upstream runner atomically publishes its final directory.  The
three downstream workers then run sequentially with one Python process each.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping


SHA256 = re.compile(r"^[0-9a-f]{64}$")


class QueueError(RuntimeError):
    pass


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def regular(path: Path, *, executable: bool = False, allow_symlink: bool = False) -> Path:
    expanded = path.expanduser()
    if expanded.is_symlink() and not allow_symlink:
        raise QueueError(f"symlink input is not allowed: {expanded}")
    resolved = expanded.resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise QueueError(f"missing, empty, or symlink input: {resolved}")
    if executable and not os.access(resolved, os.X_OK):
        raise QueueError(f"not executable: {resolved}")
    return resolved


def read_json(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(regular(path).read_text(encoding="utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise QueueError(f"invalid JSON at {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise QueueError(f"JSON object required: {path}")
    return payload


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def pid_alive(path: Path) -> int:
    try:
        text = regular(path).read_text(encoding="utf-8").strip()
        pid = int(text)
    except (OSError, UnicodeError, ValueError, QueueError) as exc:
        raise QueueError(f"cannot read upstream controller PID: {path}") from exc
    if pid < 2:
        raise QueueError(f"unsafe upstream controller PID: {pid}")
    try:
        os.kill(pid, 0)
    except ProcessLookupError as exc:
        raise QueueError(f"upstream controller PID {pid} exited before atomic publication") from exc
    except PermissionError as exc:
        raise QueueError(f"cannot inspect upstream controller PID {pid}") from exc
    return pid


def run_worker(
    command: list[str], stdout_path: Path, stderr_path: Path, *, environment: dict[str, str] | None = None
) -> None:
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
        env=environment,
    )
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise QueueError(
            f"worker {Path(command[1]).name} exited {completed.returncode}; inspect queue logs"
        )


def validate_figure_bundle(directory: Path, basename: str, expected_units: int) -> tuple[int, str]:
    validation_path = regular(directory / f"{basename}.validation.json")
    validation = read_json(validation_path)
    if validation.get("status") != "pass" or validation.get("assembly_unit_count") != expected_units:
        raise QueueError("NLR figure validation does not match the exact primary cohort")
    manifest_path = regular(directory / f"{basename}.manifest.json")
    manifest = read_json(manifest_path)
    if manifest.get("bundle_basename") != basename:
        raise QueueError("NLR figure manifest basename mismatch")
    outputs = manifest.get("outputs")
    if not isinstance(outputs, list) or len(outputs) != 5:
        raise QueueError("NLR figure manifest must contain exactly five checksummed outputs")
    observed: set[str] = set()
    for item in outputs:
        if not isinstance(item, dict):
            raise QueueError("Invalid NLR figure manifest output row")
        name, digest = item.get("basename"), item.get("sha256")
        if (
            not isinstance(name, str) or Path(name).name != name or name in observed
            or not isinstance(digest, str) or SHA256.fullmatch(digest) is None
        ):
            raise QueueError("Invalid or duplicate NLR figure manifest output")
        path = regular(directory / name)
        if sha256_file(path) != digest or item.get("bytes") != path.stat().st_size:
            raise QueueError(f"NLR figure checksum/size mismatch for {name}")
        observed.add(name)
    inventory = {path.name for path in directory.iterdir() if path.is_file()}
    if inventory != observed | {manifest_path.name}:
        raise QueueError("NLR figure directory inventory is not exact")
    return len(outputs), sha256_file(manifest_path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-controller-pid-file", type=Path, required=True)
    parser.add_argument("--nlr-root", type=Path, required=True)
    parser.add_argument("--input-bundle", type=Path, required=True)
    parser.add_argument("--unit-metadata", type=Path, required=True)
    parser.add_argument("--loss-matrix", type=Path, required=True)
    parser.add_argument("--shared-positive-genes", type=Path, required=True)
    parser.add_argument("--input-output-dir", type=Path, required=True)
    parser.add_argument("--summary-output-dir", type=Path, required=True)
    parser.add_argument("--figure-output-dir", type=Path, required=True)
    parser.add_argument("--queue-root", type=Path, required=True)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--adapter", type=Path, required=True)
    parser.add_argument("--summarizer", type=Path, required=True)
    parser.add_argument("--figure-worker", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--figure-basename", default="primary_nlr_nonshared")
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    parser.add_argument("--expected-shared-positive", type=int, default=68)
    parser.add_argument("--expected-worker-threads", type=int, default=8)
    parser.add_argument("--expected-jar-sha256", required=True)
    parser.add_argument("--expected-motifs-sha256", required=True)
    parser.add_argument("--expected-store-sha256", required=True)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args(argv)
    lock_handle = None
    try:
        if args.poll_seconds < 1 or min(
            args.expected_units, args.expected_reference_genes, args.expected_worker_threads
        ) < 1 or args.expected_shared_positive < 0:
            raise QueueError("expected counts, workers, and poll interval are invalid")
        expected_hashes = [
            args.expected_jar_sha256, args.expected_motifs_sha256, args.expected_store_sha256
        ]
        if any(SHA256.fullmatch(value.lower()) is None for value in expected_hashes):
            raise QueueError("expected tool checksums must be lowercase SHA-256 values")
        python = regular(args.python, executable=True, allow_symlink=True)
        adapter = regular(args.adapter)
        summarizer = regular(args.summarizer)
        figure_worker = regular(args.figure_worker)
        for path in (args.unit_metadata, args.loss_matrix, args.shared_positive_genes):
            regular(path)
        if args.input_bundle.is_symlink() or not args.input_bundle.is_dir():
            raise QueueError("plain NLR input bundle is missing or a symlink")
        for output in (args.input_output_dir, args.summary_output_dir, args.figure_output_dir):
            if os.path.lexists(output):
                raise QueueError(f"refusing existing downstream NLR output: {output}")
        if args.source_root.is_symlink() or not args.source_root.is_dir():
            raise QueueError("source root is missing or a symlink")

        queue = args.queue_root.expanduser().resolve()
        queue.mkdir(parents=True, exist_ok=True)
        lock_handle = (queue / "controller.lock").open("w")
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise QueueError("another primary NLR summary controller owns this queue") from exc
        state_path = queue / "state.json"
        state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "primary_nlr_summary_queue",
            "status": "waiting_for_atomic_nlr_batch",
            "started_at_utc": now(),
            "expected_units": args.expected_units,
            "maximum_downstream_worker_processes": 1,
        }
        write_json(state_path, state)
        while True:
            if args.nlr_root.is_symlink():
                raise QueueError("upstream NLR output root is a symlink")
            if args.nlr_root.is_dir():
                break
            state["upstream_controller_pid"] = pid_alive(args.upstream_controller_pid_file)
            write_json(state_path, state)
            time.sleep(args.poll_seconds)

        state["status"] = "preparing_primary_nonshared_nlr_inputs"
        state.pop("upstream_controller_pid", None)
        write_json(state_path, state)
        run_worker(
            [
                str(python), str(adapter), "--nlr-root", str(args.nlr_root),
                "--input-bundle", str(args.input_bundle), "--unit-metadata", str(args.unit_metadata),
                "--loss-matrix", str(args.loss_matrix), "--shared-positive-genes",
                str(args.shared_positive_genes), "--output-dir", str(args.input_output_dir),
                "--expected-units", str(args.expected_units), "--expected-reference-genes",
                str(args.expected_reference_genes), "--expected-shared-positive",
                str(args.expected_shared_positive), "--expected-worker-threads",
                str(args.expected_worker_threads), "--expected-jar-sha256",
                args.expected_jar_sha256.lower(), "--expected-motifs-sha256",
                args.expected_motifs_sha256.lower(), "--expected-store-sha256",
                args.expected_store_sha256.lower(),
            ],
            queue / "adapter.stdout.log", queue / "adapter.stderr.log",
        )
        adapter_report = read_json(args.input_output_dir / "validation.json")
        if (
            adapter_report.get("status") != "PASS_PRIMARY_NLR_SUMMARY_INPUTS"
            or adapter_report.get("assembly_unit_count") != args.expected_units
            or adapter_report.get("reference_gene_count") != args.expected_reference_genes
            or adapter_report.get("shared_positive_complete_gene_count")
            != args.expected_shared_positive
        ):
            raise QueueError("primary NLR input adapter did not close to expected counts")

        state["status"] = "summarizing_primary_nlr_repertoire_and_loss"
        write_json(state_path, state)
        run_worker(
            [
                str(python), str(summarizer), "--metadata",
                str(args.input_output_dir / "assembly_units.tsv"), "--repertoire-counts",
                str(args.input_output_dir / "repertoire_counts.tsv"), "--positive-loss-calls",
                str(args.input_output_dir / "positive_reference_nlr_loss_calls.tsv"),
                "--callable-denominators",
                str(args.input_output_dir / "callable_reference_nlr_denominators.tsv"),
                "--analysis-cohort", "primary_23_units_nonshared_v1", "--cohort-role",
                "primary", "--output-dir", str(args.summary_output_dir),
            ],
            queue / "summarizer.stdout.log", queue / "summarizer.stderr.log",
        )
        summary_report = read_json(args.summary_output_dir / "validation.json")
        if (
            summary_report.get("status") != "pass"
            or summary_report.get("analysis_cohort") != "primary_23_units_nonshared_v1"
            or summary_report.get("cohort_role") != "primary"
            or summary_report.get("assembly_unit_count") != args.expected_units
            or summary_report.get("denominator_input_mode") != "catalog"
        ):
            raise QueueError("primary NLR summary validation is not exact")

        state["status"] = "rendering_primary_nlr_figure"
        write_json(state_path, state)
        environment = os.environ.copy()
        environment["PYTHONPATH"] = (
            str(args.source_root)
            if not environment.get("PYTHONPATH")
            else str(args.source_root) + os.pathsep + environment["PYTHONPATH"]
        )
        run_worker(
            [
                str(python), str(figure_worker), "--unit-summary",
                str(args.summary_output_dir / "nlr_unit_summary.tsv"), "--output-dir",
                str(args.figure_output_dir), "--basename", args.figure_basename,
                "--dpi", "300", "--abbreviate-genus",
            ],
            queue / "figure.stdout.log", queue / "figure.stderr.log", environment=environment,
        )
        figure_checksums, figure_manifest_sha = validate_figure_bundle(
            args.figure_output_dir, args.figure_basename, args.expected_units
        )
        state.update(
            {
                "status": "PASS",
                "finished_at_utc": now(),
                "reference_nlr_gene_count": adapter_report.get("reference_nlr_gene_count"),
                "nonshared_reference_nlr_gene_count": adapter_report.get(
                    "nonshared_reference_nlr_gene_count"
                ),
                "positive_reference_nlr_loss_call_count": summary_report.get(
                    "positive_reference_nlr_loss_call_count"
                ),
                "callable_reference_nlr_denominator_sum": summary_report.get(
                    "callable_reference_nlr_denominator_sum"
                ),
                "figure_output_checksum_count": figure_checksums,
                "figure_manifest_sha256": figure_manifest_sha,
            }
        )
        write_json(state_path, state)
        print("PASS_PRIMARY_NLR_SUMMARY_QUEUE")
        return 0
    except (OSError, UnicodeError, ValueError, TypeError, KeyError, QueueError) as exc:
        try:
            queue = args.queue_root.expanduser().resolve()
            queue.mkdir(parents=True, exist_ok=True)
            state_path = queue / "state.json"
            try:
                current = read_json(state_path)
            except Exception:
                current = {"schema_version": 1, "workflow": "primary_nlr_summary_queue"}
            current.update({"status": "ERROR", "finished_at_utc": now(), "error": str(exc)})
            write_json(state_path, current)
        except Exception:
            pass
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    finally:
        if lock_handle is not None:
            lock_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
