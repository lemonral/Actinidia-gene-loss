#!/usr/bin/env python3
"""Run or validate a fail-closed sequential SynOrths queue.

The controller is intended for ``nohup`` execution on a compute server.  An
exclusive non-blocking lock prevents duplicate queues.  Existing completed
runs are never relaunched: their status, provenance, hashes, parameters,
current inputs, pair columns, coordinates and duplicate counts are validated
before the controller advances to the next row.
"""

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


REQUIRED = ("unit", "query_protein", "query_coords", "output_dir", "output_name")
ZERO_METRICS = (
    "duplicate_pair_rows",
    "query_fasta_duplicate_id_records",
    "reference_fasta_duplicate_id_records",
    "query_coordinate_duplicate_id_rows",
    "reference_coordinate_duplicate_id_rows",
    "query_anchor_ids_absent_from_fasta_count",
    "reference_anchor_ids_absent_from_fasta_count",
    "query_anchor_ids_absent_from_coordinates_count",
    "reference_anchor_ids_absent_from_coordinates_count",
    "query_coordinate_ids_absent_from_fasta_count",
)


class QueueError(RuntimeError):
    pass


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise QueueError(f"{path}: missing header")
        missing = [field for field in REQUIRED if field not in reader.fieldnames]
        if missing:
            raise QueueError(f"{path}: missing columns: {', '.join(missing)}")
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
            if any((value or "").strip() for value in row.values())
        ]
    if not rows:
        raise QueueError(f"{path}: queue is empty")
    units = [row["unit"] for row in rows]
    if len(units) != len(set(units)):
        raise QueueError(f"{path}: duplicate unit")
    return rows


def resolve(data_root: Path, value: str, *, require_file: bool = False) -> Path:
    candidate = Path(value)
    path = (candidate if candidate.is_absolute() else data_root / candidate).resolve()
    if require_file and (not path.is_file() or path.stat().st_size == 0):
        raise QueueError(f"missing or empty input: {path}")
    return path


def write_state(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def read_json(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise QueueError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(payload, dict):
        raise QueueError(f"{path}: JSON root is not an object")
    return payload


def live_clean_runner(output_path: Path) -> bool:
    needle = str(output_path).encode()
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        try:
            command = (proc / "cmdline").read_bytes()
        except OSError:
            continue
        if b"run_clean_synorth.py" in command and needle in command:
            return True
    return False


def wait_if_running(status_path: Path, output_path: Path, poll_seconds: int) -> None:
    while status_path.exists():
        status = read_json(status_path).get("status")
        if status != "running":
            return
        if not live_clean_runner(output_path):
            raise QueueError(f"{status_path}: status is running but no matching runner is alive")
        time.sleep(poll_seconds)


def validate_completed(
    *,
    row: dict[str, str],
    data_root: Path,
    python: Path,
    summarize_script: Path,
    reference_protein: Path,
    reference_coords: Path,
    allowed_reference_coordinate_only_ids: set[str],
    blast_bin: Path,
) -> dict[str, object]:
    unit = row["unit"]
    query_protein = resolve(data_root, row["query_protein"], require_file=True)
    query_coords = resolve(data_root, row["query_coords"], require_file=True)
    output_dir = resolve(data_root, row["output_dir"])
    output_path = output_dir / row["output_name"]
    status_path = output_dir / f"{row['output_name']}.status.json"
    provenance_path = output_dir / f"{row['output_name']}.provenance.json"
    status = read_json(status_path)
    provenance = read_json(provenance_path)
    if status.get("status") != "completed" or status.get("exit_code") != 0:
        raise QueueError(f"{unit}: run status is not completed/zero")
    if provenance.get("status") != "completed" or provenance.get("exit_code") != 0:
        raise QueueError(f"{unit}: provenance is not completed/zero")
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise QueueError(f"{unit}: output missing or empty")
    output_record = status.get("output")
    if not isinstance(output_record, dict) or output_record.get("sha256") != sha256(output_path):
        raise QueueError(f"{unit}: output checksum does not match completed status")
    parameters = provenance.get("parameters")
    expected_parameters = {"m": 20, "n": 100, "r": 0.2, "blast_threads": 2}
    if parameters != expected_parameters:
        raise QueueError(f"{unit}: unexpected SynOrths parameters: {parameters}")
    external_tools = provenance.get("external_blast_tools_final")
    if external_tools is not None:
        if provenance.get("external_blast_tools_unchanged") is not True:
            raise QueueError(f"{unit}: BLAST tools were not stable during the run")
        expected_blast = {
            "blastp": blast_bin.resolve() / "blastp",
            "makeblastdb": blast_bin.resolve() / "makeblastdb",
        }
        if not isinstance(external_tools, list):
            raise QueueError(f"{unit}: malformed external BLAST tool provenance")
        observed_blast = {
            str(item.get("role")): item
            for item in external_tools
            if isinstance(item, dict)
        }
        for role, path in expected_blast.items():
            record = observed_blast.get(role)
            if record is None or Path(str(record.get("path", ""))).resolve() != path:
                raise QueueError(f"{unit}: {role} is not bound to --blast-bin")
            if record.get("sha256") != sha256(path):
                raise QueueError(f"{unit}: current {role} checksum changed")

    expected_inputs = {
        "query_protein": query_protein,
        "query_coords": query_coords,
        "reference_protein": reference_protein,
        "reference_coords": reference_coords,
    }
    input_records = provenance.get("inputs")
    if not isinstance(input_records, list):
        raise QueueError(f"{unit}: provenance inputs are missing")
    observed: dict[str, dict[str, object]] = {}
    for item in input_records:
        if isinstance(item, dict) and isinstance(item.get("role"), str):
            observed[item["role"]] = item
    for role, path in expected_inputs.items():
        record = observed.get(role)
        if record is None or Path(str(record.get("path", ""))).resolve() != path:
            raise QueueError(f"{unit}: {role} path is not bound to the queue manifest")
        if record.get("sha256") != sha256(path):
            raise QueueError(f"{unit}: current {role} checksum changed")

    summary_json = output_dir / f"{unit}.synorth_summary.json"
    summary_tsv = output_dir / f"{unit}.synorth_summary.tsv"
    command = [
        str(python), str(summarize_script), "--sample", unit,
        "--pairs", str(output_path), "--query-fasta", str(query_protein),
        "--reference-fasta", str(reference_protein), "--query-coords", str(query_coords),
        "--reference-coords", str(reference_coords), "--query-column", "1",
        "--reference-column", "5", "--output-json", str(summary_json),
        "--output-tsv", str(summary_tsv),
    ]
    completed = subprocess.run(command, check=False)
    if completed.returncode:
        raise QueueError(f"{unit}: summarize_synorth failed with {completed.returncode}")
    summary = read_json(summary_json)
    metrics = summary.get("metrics")
    if not isinstance(metrics, dict):
        raise QueueError(f"{unit}: summary metrics missing")
    failures = {name: metrics.get(name) for name in ZERO_METRICS if metrics.get(name) != 0}
    if failures:
        raise QueueError(f"{unit}: reconciliation failed: {failures}")
    absent_ids = summary.get("absent_ids")
    if not isinstance(absent_ids, dict):
        raise QueueError(f"{unit}: summary absent-ID section missing")
    observed_reference_coordinate_only = set(
        absent_ids.get("reference_coordinate_ids_absent_from_fasta") or []
    )
    if observed_reference_coordinate_only != allowed_reference_coordinate_only_ids:
        raise QueueError(
            f"{unit}: reference coordinate-only IDs differ from the frozen allowed set; "
            f"observed={len(observed_reference_coordinate_only)}, "
            f"allowed={len(allowed_reference_coordinate_only_ids)}"
        )
    if metrics.get("query_anchors_in_coordinates") != metrics.get("unique_query_anchors"):
        raise QueueError(f"{unit}: query anchor/coordinate reconciliation failed")
    if metrics.get("reference_anchors_in_coordinates") != metrics.get("unique_reference_anchors"):
        raise QueueError(f"{unit}: reference anchor/coordinate reconciliation failed")
    return {
        "unit": unit,
        "output_sha256": sha256(output_path),
        "unique_pairs": metrics.get("unique_pairs"),
        "query_coverage_percent": metrics.get("query_coverage_percent"),
        "reference_coverage_percent": metrics.get("reference_coverage_percent"),
        "validated_at_utc": utc_now(),
    }


def run_one(
    *,
    row: dict[str, str],
    data_root: Path,
    python: Path,
    runner: Path,
    synorth_dir: Path,
    reference_protein: Path,
    reference_coords: Path,
    blast_bin: Path,
) -> int:
    query_protein = resolve(data_root, row["query_protein"], require_file=True)
    query_coords = resolve(data_root, row["query_coords"], require_file=True)
    output_dir = resolve(data_root, row["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [
        str(python), str(runner), "--synorth-dir", str(synorth_dir),
        "--blast-bin", str(blast_bin),
        "--query-protein", str(query_protein), "--query-coords", str(query_coords),
        "--reference-protein", str(reference_protein), "--reference-coords", str(reference_coords),
        "--output-dir", str(output_dir), "--output-name", row["output_name"],
        "--m", "20", "--n", "100", "--r", "0.2", "--blast-threads", "2",
    ]
    return subprocess.run(command, check=False).returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--runner", required=True, type=Path)
    parser.add_argument("--summarize-script", required=True, type=Path)
    parser.add_argument("--synorth-dir", required=True, type=Path)
    parser.add_argument("--blast-bin", required=True, type=Path)
    parser.add_argument("--reference-protein", required=True, type=Path)
    parser.add_argument("--reference-coords", required=True, type=Path)
    parser.add_argument(
        "--allowed-reference-coordinate-only-ids",
        type=Path,
        help="optional one-ID-per-line allowlist for reference coordinates without proteins",
    )
    parser.add_argument("--state", required=True, type=Path)
    parser.add_argument("--lock", required=True, type=Path)
    parser.add_argument("--poll-seconds", type=int, default=60)
    args = parser.parse_args()
    if args.poll_seconds < 5:
        raise SystemExit("--poll-seconds must be at least 5")

    args.lock.parent.mkdir(parents=True, exist_ok=True)
    lock_handle = args.lock.open("w", encoding="utf-8")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print(f"ERROR: another queue owns {args.lock}", file=sys.stderr)
        return 3
    lock_handle.write(f"pid\t{os.getpid()}\nstarted_at_utc\t{utc_now()}\n")
    lock_handle.flush()

    rows = load_rows(args.manifest)
    allowed_reference_coordinate_only_ids: set[str] = set()
    if args.allowed_reference_coordinate_only_ids:
        if not args.allowed_reference_coordinate_only_ids.is_file():
            raise SystemExit(
                f"missing --allowed-reference-coordinate-only-ids file: "
                f"{args.allowed_reference_coordinate_only_ids}"
            )
        allowed_reference_coordinate_only_ids = {
            line.strip()
            for line in args.allowed_reference_coordinate_only_ids.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
    validated: list[dict[str, object]] = []
    base_state: dict[str, object] = {
        "schema_version": 1,
        "controller_pid": os.getpid(),
        "manifest_sha256": sha256(args.manifest),
        "queue_size": len(rows),
        "started_at_utc": utc_now(),
    }
    try:
        for index, row in enumerate(rows, start=1):
            unit = row["unit"]
            output_dir = resolve(args.data_root, row["output_dir"])
            output_path = output_dir / row["output_name"]
            status_path = output_dir / f"{row['output_name']}.status.json"
            write_state(args.state, {**base_state, "status": "checking", "current_unit": unit, "queue_index": index, "validated": validated})
            wait_if_running(status_path, output_path, args.poll_seconds)
            should_run = not status_path.exists()
            if status_path.exists():
                prior_status = str(read_json(status_path).get("status", ""))
                if prior_status == "completed":
                    should_run = False
                elif prior_status == "failed_no_output" and (
                    not output_path.exists() or output_path.stat().st_size == 0
                ):
                    # Safe recovery for a wrapper that returned no biological
                    # result (for example, a missing executable on PATH).
                    should_run = True
                else:
                    raise QueueError(
                        f"{unit}: existing non-completed status {prior_status!r} "
                        "requires manual review"
                    )
            if should_run:
                write_state(args.state, {**base_state, "status": "running", "current_unit": unit, "queue_index": index, "validated": validated})
                if run_one(
                    row=row, data_root=args.data_root, python=args.python, runner=args.runner,
                    synorth_dir=args.synorth_dir, reference_protein=args.reference_protein,
                    reference_coords=args.reference_coords, blast_bin=args.blast_bin,
                ):
                    raise QueueError(f"{unit}: SynOrths runner failed")
            validated.append(
                validate_completed(
                    row=row, data_root=args.data_root, python=args.python,
                    summarize_script=args.summarize_script,
                    reference_protein=args.reference_protein,
                    reference_coords=args.reference_coords,
                    allowed_reference_coordinate_only_ids=allowed_reference_coordinate_only_ids,
                    blast_bin=args.blast_bin,
                )
            )
        write_state(args.state, {**base_state, "status": "completed", "finished_at_utc": utc_now(), "validated": validated})
        print(f"Completed and validated {len(validated)} SynOrths entries")
        return 0
    except QueueError as error:
        write_state(args.state, {**base_state, "status": "failed", "finished_at_utc": utc_now(), "error": str(error), "validated": validated})
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
