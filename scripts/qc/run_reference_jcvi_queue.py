#!/usr/bin/env python3
"""Run one-direction JCVI comparisons from a fixed reference to many units.

Each lane is sequential and lanes run in parallel.  The underlying
``run_jcvi_pair.py`` workflow writes every comparison atomically and records
exact input/tool checksums.  This controller refuses pre-existing outputs so a
completed comparison cannot be silently duplicated or overwritten.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import subprocess


HEADER = ("lane", "unit", "alias", "display_name", "accession", "protein", "bed")


class QueueError(RuntimeError):
    """Raised when the queue or one of its completed outputs is invalid."""


@dataclass(frozen=True)
class Entry:
    lane: str
    unit: str
    alias: str
    display_name: str
    accession: str
    protein: Path
    bed: Path


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise QueueError(f"{label} is missing or empty: {resolved}")
    return resolved


def safe_relative(root: Path, raw: str, label: str) -> Path:
    relative = Path(raw.strip())
    if not raw.strip() or relative.is_absolute() or ".." in relative.parts:
        raise QueueError(f"{label} must be a safe relative path")
    resolved = (root / relative).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise QueueError(f"{label} escapes the data root") from error
    return regular(resolved, label)


def read_manifest(path: Path, root: Path) -> list[Entry]:
    with regular(path, "queue manifest").open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != HEADER:
            raise QueueError(f"manifest header must be exactly {HEADER}")
        rows = list(reader)
    if not rows:
        raise QueueError("queue manifest is empty")
    entries: list[Entry] = []
    units: set[str] = set()
    aliases: set[str] = set()
    for line_number, row in enumerate(rows, 2):
        lane = row["lane"].strip()
        unit = row["unit"].strip()
        alias = row["alias"].strip()
        if not lane or not unit or not alias:
            raise QueueError(f"line {line_number}: lane, unit, and alias are required")
        if any(character.isspace() for character in lane + unit + alias):
            raise QueueError(f"line {line_number}: unsafe whitespace in lane/unit/alias")
        if not alias[0].isalpha() or not alias.replace("_", "").isalnum():
            raise QueueError(f"line {line_number}: alias is not JCVI-safe")
        if unit in units or alias in aliases:
            raise QueueError(f"line {line_number}: duplicate unit or alias")
        units.add(unit)
        aliases.add(alias)
        display_name = row["display_name"].strip()
        accession = row["accession"].strip()
        if not display_name or not accession:
            raise QueueError(f"line {line_number}: display name and accession are required")
        entries.append(
            Entry(
                lane=lane,
                unit=unit,
                alias=alias,
                display_name=display_name,
                accession=accession,
                protein=safe_relative(root, row["protein"], "query protein"),
                bed=safe_relative(root, row["bed"], "query BED"),
            )
        )
    return entries


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(temporary, path)


def validate_output(
    output: Path, reference_id: str, entry: Entry, threads: int
) -> dict[str, object]:
    manifest = regular(output / "run_manifest.json", "JCVI run manifest")
    coverage = regular(output / "jcvi_bidirectional_coverage.json", "JCVI coverage JSON")
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    if (
        payload.get("reference_id") != reference_id
        or payload.get("query_id") != entry.alias
        or payload.get("threads") != threads
    ):
        raise QueueError(f"run identity closure failed: {output}")
    anchors = regular(
        output / "work" / f"{reference_id}.{entry.alias}.anchors", "raw anchors"
    )
    return {
        "unit": entry.unit,
        "lane": entry.lane,
        "output": output.name,
        "manifest_sha256": sha256(manifest),
        "coverage_sha256": sha256(coverage),
        "anchors_sha256": sha256(anchors),
        "finished_at_utc": now(),
    }


def parser() -> argparse.ArgumentParser:
    argument_parser = argparse.ArgumentParser(description=__doc__)
    argument_parser.add_argument("--manifest", required=True, type=Path)
    argument_parser.add_argument("--data-root", required=True, type=Path)
    argument_parser.add_argument("--output-root", required=True, type=Path)
    argument_parser.add_argument("--state", required=True, type=Path)
    argument_parser.add_argument("--lock", required=True, type=Path)
    argument_parser.add_argument("--python", required=True, type=Path)
    argument_parser.add_argument("--runner", required=True, type=Path)
    argument_parser.add_argument("--reference-id", required=True)
    argument_parser.add_argument("--reference-display-name", required=True)
    argument_parser.add_argument("--reference-protein", required=True, type=Path)
    argument_parser.add_argument("--reference-bed", required=True, type=Path)
    argument_parser.add_argument("--allowed-reference-bed-only-ids", type=Path)
    argument_parser.add_argument("--threads", type=int, default=4)
    argument_parser.add_argument("--max-lanes", type=int, default=2)
    return argument_parser


def main() -> int:
    args = parser().parse_args()
    try:
        if args.threads < 1 or args.threads > 8:
            raise QueueError("--threads must be between 1 and 8")
        if args.max_lanes < 1 or args.max_lanes > 3:
            raise QueueError("--max-lanes must be between 1 and 3")
        root = args.data_root.expanduser().resolve()
        entries = read_manifest(args.manifest, root)
        lanes: dict[str, list[Entry]] = {}
        for entry in entries:
            lanes.setdefault(entry.lane, []).append(entry)
        if len(lanes) > args.max_lanes:
            raise QueueError("manifest exceeds --max-lanes")

        python = regular(args.python, "Python executable")
        runner = regular(args.runner, "JCVI pair runner")
        reference_protein = regular(args.reference_protein, "reference protein")
        reference_bed = regular(args.reference_bed, "reference BED")
        allowlist = (
            regular(args.allowed_reference_bed_only_ids, "reference BED-only allow-list")
            if args.allowed_reference_bed_only_ids
            else None
        )
        output_root = args.output_root.expanduser().resolve()
        state_path = args.state.expanduser().resolve()
        lock_path = args.lock.expanduser().resolve()
        if output_root.exists() or state_path.exists() or lock_path.exists():
            raise QueueError("output root, state, and lock must all be absent")
        output_root.mkdir(parents=True)
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(f"pid\t{os.getpid()}\n")
            handle.write(f"manifest_sha256\t{sha256(regular(args.manifest, 'manifest'))}\n")

        state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "fixed_reference_jcvi_queue",
            "status": "running",
            "started_at_utc": now(),
            "reference_id": args.reference_id,
            "reference_protein_sha256": sha256(reference_protein),
            "reference_bed_sha256": sha256(reference_bed),
            "manifest_sha256": sha256(regular(args.manifest, "manifest")),
            "threads_per_pair": args.threads,
            "max_lanes": args.max_lanes,
            "completed": [],
        }
        write_json(state_path, state)

        def run_entry(entry: Entry) -> dict[str, object]:
            output = output_root / entry.unit
            command = [
                str(python),
                str(runner),
                "--reference-id",
                args.reference_id,
                "--query-id",
                entry.alias,
                "--reference-display-name",
                args.reference_display_name,
                "--query-display-name",
                entry.display_name,
                "--query-accession",
                entry.accession,
                "--reference-protein",
                str(reference_protein),
                "--reference-bed",
                str(reference_bed),
                "--query-protein",
                str(entry.protein),
                "--query-bed",
                str(entry.bed),
                "--output-dir",
                str(output),
                "--threads",
                str(args.threads),
                "--python-bin",
                str(python),
                "--minimum-block-size",
                "4",
            ]
            if allowlist:
                command.extend(["--allowed-reference-bed-only-ids", str(allowlist)])
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            if completed.returncode != 0:
                raise QueueError(
                    f"{entry.unit} failed with exit {completed.returncode}: "
                    f"{completed.stderr.strip()}"
                )
            return validate_output(output, args.reference_id, entry, args.threads)

        def run_lane(lane_entries: list[Entry]) -> list[dict[str, object]]:
            return [run_entry(entry) for entry in lane_entries]

        errors: list[str] = []
        records: list[dict[str, object]] = []
        with ThreadPoolExecutor(max_workers=len(lanes)) as executor:
            futures = {
                executor.submit(run_lane, lane_entries): lane
                for lane, lane_entries in sorted(lanes.items())
            }
            for future in as_completed(futures):
                try:
                    lane_records = future.result()
                    records.extend(lane_records)
                    state["completed"] = sorted(records, key=lambda row: str(row["unit"]))
                    write_json(state_path, state)
                except Exception as error:  # noqa: BLE001 - recorded in atomic state
                    errors.append(f"{futures[future]}: {error}")
        if errors:
            state.update({"status": "failed", "finished_at_utc": now(), "errors": errors})
            write_json(state_path, state)
            raise QueueError("; ".join(errors))
        if len(records) != len(entries):
            raise QueueError("completed comparison count does not close against manifest")
        state.update({"status": "PASS", "finished_at_utc": now()})
        write_json(state_path, state)
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, QueueError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
