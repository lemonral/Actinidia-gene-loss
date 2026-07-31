#!/usr/bin/env python3
"""Run four-direction JCVI evidence for new chromosome units in two lanes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


HEADER = ("lane", "unit", "alias", "display_name", "accession", "protein", "bed", "output_root")


class QueueError(RuntimeError):
    pass


@dataclass(frozen=True)
class Entry:
    lane: str
    unit: str
    alias: str
    display_name: str
    accession: str
    protein: Path
    bed: Path
    output_root: Path


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
    if not resolved.is_file() or resolved.is_symlink() or resolved.stat().st_size == 0:
        raise QueueError(f"{label} is missing, empty, or a symlink: {resolved}")
    return resolved


def relative(raw: str, label: str) -> Path:
    value = raw.strip()
    path = Path(value)
    if not value or path.is_absolute() or ".." in path.parts:
        raise QueueError(f"{label} must be a safe relative path")
    return path


def beneath(root: Path, path: Path, label: str) -> Path:
    resolved = (root / path).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise QueueError(f"{label} escapes data root") from error
    return resolved


def read_manifest(path: Path, root: Path) -> list[Entry]:
    manifest = regular(path, "queue manifest")
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != HEADER:
            raise QueueError(f"manifest header must be exactly {HEADER}")
        rows = list(reader)
    if not rows:
        raise QueueError("queue manifest is empty")
    entries: list[Entry] = []
    units: set[str] = set()
    aliases: set[str] = set()
    outputs: set[Path] = set()
    for line, row in enumerate(rows, 2):
        lane, unit, alias = (row[key].strip() for key in ("lane", "unit", "alias"))
        if not lane or not unit or not alias or not alias[0].isalpha() or not alias.replace("_", "").isalnum():
            raise QueueError(f"line {line}: unsafe lane, unit, or alias")
        if any(character.isspace() for character in unit + alias):
            raise QueueError(f"line {line}: unit/alias contains whitespace")
        if unit in units or alias in aliases:
            raise QueueError(f"line {line}: duplicate unit or alias")
        units.add(unit)
        aliases.add(alias)
        display = row["display_name"].strip()
        accession = row["accession"].strip()
        if not display or not accession:
            raise QueueError(f"line {line}: empty display name or accession")
        protein = beneath(root, relative(row["protein"], "protein"), "protein")
        bed = beneath(root, relative(row["bed"], "BED"), "BED")
        output = beneath(root, relative(row["output_root"], "output root"), "output root")
        if output in outputs:
            raise QueueError(f"line {line}: duplicate output root")
        outputs.add(output)
        entries.append(Entry(lane, unit, alias, display, accession, protein, bed, output))
    return entries


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(temporary, path)


def validate_output(path: Path, reference_id: str, query_id: str, threads: int) -> dict[str, object]:
    manifest_path = regular(path / "run_manifest.json", "JCVI run manifest")
    coverage = regular(path / "jcvi_bidirectional_coverage.json", "JCVI coverage JSON")
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        payload.get("reference_id") != reference_id
        or payload.get("query_id") != query_id
        or payload.get("threads") != threads
    ):
        raise QueueError(f"JCVI manifest identity closure failed: {path}")
    work = path / "work"
    anchors = regular(work / f"{reference_id}.{query_id}.anchors", "raw anchors")
    return {
        "directory": path.name,
        "reference_id": reference_id,
        "query_id": query_id,
        "manifest_sha256": sha256(manifest_path),
        "coverage_sha256": sha256(coverage),
        "raw_anchors_sha256": sha256(anchors),
    }


def run_pair(
    *, entry: Entry, reference: dict[str, str | Path], target_is_reference: bool,
    runner: Path, python: Path, threads: int, output: Path,
) -> dict[str, object]:
    if target_is_reference:
        reference_id, reference_display, reference_protein, reference_bed = (
            entry.alias, entry.display_name, entry.protein, entry.bed
        )
        query_id = str(reference["id"])
        query_display = str(reference["display"])
        query_accession = str(reference["accession"])
        query_protein = Path(reference["protein"])
        query_bed = Path(reference["bed"])
    else:
        reference_id = str(reference["id"])
        reference_display = str(reference["display"])
        reference_protein = Path(reference["protein"])
        reference_bed = Path(reference["bed"])
        query_id, query_display, query_accession, query_protein, query_bed = (
            entry.alias, entry.display_name, entry.accession, entry.protein, entry.bed
        )
    command = [
        str(python), str(runner),
        "--reference-id", reference_id,
        "--query-id", query_id,
        "--reference-display-name", reference_display,
        "--query-display-name", query_display,
        "--query-accession", query_accession,
        "--reference-protein", str(reference_protein),
        "--reference-bed", str(reference_bed),
        "--query-protein", str(query_protein),
        "--query-bed", str(query_bed),
        "--output-dir", str(output),
        "--threads", str(threads),
        "--python-bin", str(python),
        "--minimum-block-size", "4",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise QueueError(
            f"JCVI pair failed for {entry.unit}/{output.name}: {completed.stderr.strip()}"
        )
    return validate_output(output, reference_id, query_id, threads)


class Controller:
    def __init__(self, args: argparse.Namespace, entries: list[Entry]) -> None:
        self.args = args
        self.entries = entries
        self.lock = threading.Lock()
        self.state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "new_unit_four_direction_jcvi_queue",
            "status": "running",
            "started_at_utc": now(),
            "manifest_sha256": sha256(args.manifest),
            "threads_per_pair": args.threads,
            "max_lanes": args.max_lanes,
            "completed": [],
        }
        write_json(args.state, self.state)

    def append(self, record: dict[str, object]) -> None:
        with self.lock:
            completed = self.state["completed"]
            assert isinstance(completed, list)
            completed.append(record)
            write_json(self.args.state, self.state)

    def run_entry(self, entry: Entry) -> dict[str, object]:
        if entry.output_root.exists():
            raise QueueError(f"refusing existing unit output root: {entry.output_root}")
        entry.output_root.mkdir(parents=True)
        references = (
            {
                "id": "HY4A", "display": "Actinidia chinensis HY4A",
                "accession": "Hongyang_v4_HY4A",
                "protein": self.args.hy4a_protein, "bed": self.args.hy4a_bed,
            },
            {
                "id": "HY4P", "display": "Actinidia chinensis HY4P",
                "accession": "Hongyang_v4_HY4P",
                "protein": self.args.hy4p_protein, "bed": self.args.hy4p_bed,
            },
        )
        results = []
        for reference in references:
            suffix = str(reference["id"]).lower()
            for target_is_reference, name in (
                (False, f"{suffix}_reference_target_query"),
                (True, f"target_reference_{suffix}_query"),
            ):
                results.append(
                    run_pair(
                        entry=entry,
                        reference=reference,
                        target_is_reference=target_is_reference,
                        runner=self.args.runner,
                        python=self.args.python,
                        threads=self.args.threads,
                        output=entry.output_root / name,
                    )
                )
        return {"unit": entry.unit, "lane": entry.lane, "runs": results, "finished_at_utc": now()}

    def run_lane(self, entries: list[Entry]) -> None:
        for entry in entries:
            self.append(self.run_entry(entry))

    def run(self) -> None:
        lanes: dict[str, list[Entry]] = {}
        for entry in self.entries:
            lanes.setdefault(entry.lane, []).append(entry)
        if len(lanes) > self.args.max_lanes:
            raise QueueError("manifest exceeds maximum JCVI lanes")
        errors = []
        with ThreadPoolExecutor(max_workers=len(lanes)) as executor:
            futures = {executor.submit(self.run_lane, rows): lane for lane, rows in sorted(lanes.items())}
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as error:
                    errors.append(f"{futures[future]}: {error}")
        if errors:
            raise QueueError("; ".join(errors))
        completed = self.state["completed"]
        if not isinstance(completed, list) or len(completed) != len(self.entries):
            raise QueueError("completed unit count does not close against manifest")
        self.state.update({"status": "PASS", "finished_at_utc": now()})
        write_json(self.args.state, self.state)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--data-root", required=True, type=Path)
    p.add_argument("--python", required=True, type=Path)
    p.add_argument("--runner", required=True, type=Path)
    p.add_argument("--hy4a-protein", required=True, type=Path)
    p.add_argument("--hy4a-bed", required=True, type=Path)
    p.add_argument("--hy4p-protein", required=True, type=Path)
    p.add_argument("--hy4p-bed", required=True, type=Path)
    p.add_argument("--state", required=True, type=Path)
    p.add_argument("--lock", required=True, type=Path)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--max-lanes", type=int, default=2)
    return p


def main() -> int:
    args = parser().parse_args()
    controller: Controller | None = None
    try:
        if args.threads != 4 or args.max_lanes < 1 or args.max_lanes > 2:
            raise QueueError("frozen queue requires four threads and one or two lanes")
        args.data_root = args.data_root.expanduser().resolve()
        args.manifest = regular(args.manifest, "queue manifest")
        args.python = regular(args.python, "Python")
        args.runner = regular(args.runner, "JCVI pair runner")
        for name in ("hy4a_protein", "hy4a_bed", "hy4p_protein", "hy4p_bed"):
            setattr(args, name, regular(getattr(args, name), name))
        args.state = args.state.expanduser().resolve()
        args.lock = args.lock.expanduser().resolve()
        if args.state.exists():
            raise QueueError(f"refusing existing state: {args.state}")
        args.lock.parent.mkdir(parents=True, exist_ok=True)
        descriptor = os.open(args.lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(f"pid\t{os.getpid()}\nmanifest_sha256\t{sha256(args.manifest)}\n")
        entries = read_manifest(args.manifest, args.data_root)
        controller = Controller(args, entries)
        controller.run()
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, subprocess.SubprocessError, QueueError) as error:
        if controller is not None:
            controller.state.update({"status": "failed", "finished_at_utc": now(), "error": str(error)})
            write_json(args.state, controller.state)
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
