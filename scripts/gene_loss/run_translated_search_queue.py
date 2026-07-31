#!/usr/bin/env python3
"""Run a detached, fail-closed translated-search queue for new Actinidia units."""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


class SearchError(RuntimeError):
    pass


REQUIRED_MANIFEST = ("unit", "target_genome", "candidate_dir", "output_dir")
HIT_COLUMNS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qframe",
    "sframe",
)
STATE_COLUMNS = (
    "unit",
    "reference_gene",
    "callable",
    "callability_reason",
    "target_chromosome",
    "target_interval_start_1based",
    "target_interval_end_1based",
    "qualifying_genome_hit_count",
    "qualifying_local_hit_count",
    "best_hit_subject",
    "best_hit_percent_identity",
    "best_hit_alignment_length",
    "best_hit_evalue",
    "best_hit_bitscore",
    "primary_state",
    "positive_loss",
    "historical_reproduction_state",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    source = path.resolve()
    if not source.is_file() or (source.stat().st_size <= 0 and not allow_empty):
        raise SearchError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def strict_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise SearchError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise SearchError(f"{path}: JSON root is not an object")
    return value


def resolve(root: Path, value: str, *, file: bool = False) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise SearchError(f"unsafe data-root-relative path: {value!r}")
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise SearchError(f"path escapes data root: {value!r}")
    if file and (not path.is_file() or path.stat().st_size <= 0):
        raise SearchError(f"missing input file: {value!r}")
    return path


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != REQUIRED_MANIFEST:
            raise SearchError("queue manifest columns differ from exact schema")
        rows = list(reader)
    if not rows or len({row["unit"] for row in rows}) != len(rows):
        raise SearchError("queue manifest is empty or has duplicate units")
    return rows


def read_fasta_ids(path: Path) -> set[str]:
    identifiers: set[str] = set()
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            identifier = raw[1:].strip().split()[0]
            if not identifier or identifier in identifiers:
                raise SearchError(f"{path.name}:{line_number}: duplicate/empty FASTA ID")
            identifiers.add(identifier)
    if not identifiers:
        raise SearchError(f"{path.name}: FASTA contains no IDs")
    return identifiers


def validate_candidate_bundle(unit: str, candidate_dir: Path, target_genome: Path) -> tuple[Path, Path, dict[str, object]]:
    manifest_path = candidate_dir / "run_manifest.json"
    candidates = candidate_dir / "candidates.tsv"
    queries = candidate_dir / "candidate_reference_cds.fasta"
    for path in (manifest_path, candidates, queries, candidate_dir / "checksums.tsv"):
        if not path.is_file() or path.stat().st_size <= 0:
            raise SearchError(f"{unit}: candidate bundle is incomplete")
    manifest = strict_json(manifest_path)
    if manifest.get("status") != "PASS" or manifest.get("unit") != unit:
        raise SearchError(f"{unit}: candidate manifest is not PASS/bound")
    inputs = manifest.get("inputs")
    outputs = manifest.get("outputs")
    if not isinstance(inputs, dict) or not isinstance(outputs, dict):
        raise SearchError(f"{unit}: candidate manifest lacks bindings")
    if inputs.get("target_genome") != binding(target_genome):
        raise SearchError(f"{unit}: target genome changed after candidate preparation")
    if outputs.get("candidates") != binding(candidates):
        raise SearchError(f"{unit}: candidate table binding mismatch")
    if outputs.get("candidate_reference_cds") != binding(queries):
        raise SearchError(f"{unit}: query FASTA binding mismatch")
    return candidates, queries, manifest


def read_candidates(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    required = {
        "unit", "reference_gene", "has_reference_cds", "target_chromosome",
        "target_interval_start_1based", "target_interval_end_1based", "callable",
        "callability_reason",
    }
    if not rows or not required.issubset(reader.fieldnames or []):
        raise SearchError(f"{path.name}: candidate table schema/rows invalid")
    genes = [row["reference_gene"] for row in rows]
    if len(genes) != len(set(genes)):
        raise SearchError(f"{path.name}: duplicate reference candidates")
    return rows


def numeric(value: str, label: str) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise SearchError(f"non-numeric {label}: {value!r}") from error
    if not math.isfinite(number):
        raise SearchError(f"non-finite {label}: {value!r}")
    return number


def parse_hits(path: Path, query_ids: set[str], subject_ids: set[str]) -> dict[str, list[dict[str, str]]]:
    hits: dict[str, list[dict[str, str]]] = defaultdict(list)
    if path.stat().st_size == 0:
        return hits
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, fieldnames=HIT_COLUMNS, delimiter="\t")
        for line_number, row in enumerate(reader, 1):
            if row["qseqid"] not in query_ids or row["sseqid"] not in subject_ids:
                raise SearchError(f"{path.name}:{line_number}: query/subject ID outside inputs")
            for field in ("pident", "evalue", "bitscore"):
                numeric(row[field], field)
            for field in ("length", "qstart", "qend", "sstart", "send", "qframe", "sframe"):
                try:
                    int(row[field])
                except ValueError as error:
                    raise SearchError(f"{path.name}:{line_number}: invalid integer {field}") from error
            hits[row["qseqid"]].append(row)
    return hits


def best_hit(rows: list[dict[str, str]]) -> dict[str, str] | None:
    if not rows:
        return None
    return min(
        rows,
        key=lambda row: (
            -numeric(row["bitscore"], "bitscore"),
            numeric(row["evalue"], "evalue"),
            -numeric(row["pident"], "pident"),
            row["sseqid"],
            min(int(row["sstart"]), int(row["send"])),
        ),
    )


def qualifies(row: dict[str, str], *, identity: float, bitscore: float, evalue: float) -> bool:
    return (
        numeric(row["pident"], "pident") >= identity
        and numeric(row["bitscore"], "bitscore") >= bitscore
        and numeric(row["evalue"], "evalue") < evalue
    )


def classify(
    *,
    unit: str,
    candidate_rows: list[dict[str, str]],
    hits: dict[str, list[dict[str, str]]],
    minimum_identity: float,
    minimum_bitscore: float,
    maximum_evalue: float,
) -> tuple[list[dict[str, str]], dict[str, int]]:
    states: list[dict[str, str]] = []
    counts: dict[str, int] = defaultdict(int)
    for candidate in candidate_rows:
        gene = candidate["reference_gene"]
        raw = hits.get(gene, [])
        qualifying = [
            row for row in raw
            if qualifies(row, identity=minimum_identity, bitscore=minimum_bitscore, evalue=maximum_evalue)
        ]
        local: list[dict[str, str]] = []
        if candidate["target_chromosome"] and candidate["target_interval_start_1based"]:
            start = int(candidate["target_interval_start_1based"])
            end = int(candidate["target_interval_end_1based"])
            local = [
                row for row in qualifying
                if row["sseqid"] == candidate["target_chromosome"]
                and max(min(int(row["sstart"]), int(row["send"])), start)
                <= min(max(int(row["sstart"]), int(row["send"])), end)
            ]
        best = best_hit(raw)
        best_qualifies = best is not None and qualifies(
            best, identity=minimum_identity, bitscore=minimum_bitscore, evalue=maximum_evalue
        )
        if candidate["callable"] != "true":
            primary = "uncertain_non_callable"
            positive = False
        elif local:
            primary = "uncertain_local_genomic_sequence_detected"
            positive = False
        else:
            primary = "positive_deleted"
            positive = True
        historical = "decayed" if best_qualifies else "deleted"
        counts[primary] += 1
        counts[f"historical_{historical}"] += 1
        states.append(
            {
                "unit": unit,
                "reference_gene": gene,
                "callable": candidate["callable"],
                "callability_reason": candidate["callability_reason"],
                "target_chromosome": candidate["target_chromosome"],
                "target_interval_start_1based": candidate["target_interval_start_1based"],
                "target_interval_end_1based": candidate["target_interval_end_1based"],
                "qualifying_genome_hit_count": str(len(qualifying)),
                "qualifying_local_hit_count": str(len(local)),
                "best_hit_subject": "" if best is None else best["sseqid"],
                "best_hit_percent_identity": "" if best is None else best["pident"],
                "best_hit_alignment_length": "" if best is None else best["length"],
                "best_hit_evalue": "" if best is None else best["evalue"],
                "best_hit_bitscore": "" if best is None else best["bitscore"],
                "primary_state": primary,
                "positive_loss": str(positive).lower(),
                "historical_reproduction_state": historical,
            }
        )
    return states, dict(sorted(counts.items()))


def run_command(command: list[str], *, stdout: Path | None = None, stderr: Path | None = None, stdin=None) -> None:
    stdout_handle = stdout.open("wb") if stdout else subprocess.DEVNULL
    stderr_handle = stderr.open("wb") if stderr else subprocess.DEVNULL
    try:
        completed = subprocess.run(command, stdin=stdin, stdout=stdout_handle, stderr=stderr_handle, check=False)
    finally:
        if stdout:
            stdout_handle.close()
        if stderr:
            stderr_handle.close()
    if completed.returncode:
        raise SearchError(f"command failed with exit code {completed.returncode}: {command[0]}")


def build_database(genome: Path, prefix: Path, makeblastdb: Path, stderr: Path) -> None:
    if genome.suffix.lower() == ".gz":
        decompressed = prefix.parent / "genome.input.fa"
        with gzip.open(genome, "rb") as source, decompressed.open("wb") as destination:
            shutil.copyfileobj(source, destination, length=1 << 20)
        try:
            run_command(
                [str(makeblastdb), "-in", str(decompressed), "-input_type", "fasta", "-dbtype", "nucl", "-parse_seqids", "-out", str(prefix)],
                stderr=stderr,
            )
        finally:
            decompressed.unlink(missing_ok=True)
    else:
        run_command(
            [str(makeblastdb), "-in", str(genome), "-input_type", "fasta", "-dbtype", "nucl", "-parse_seqids", "-out", str(prefix)],
            stderr=stderr,
        )


def run_unit(
    *,
    row: dict[str, str],
    root: Path,
    tools: dict[str, Path],
    threads: int,
    minimum_identity: float,
    minimum_bitscore: float,
    maximum_evalue: float,
    max_hsps: int,
) -> dict[str, object]:
    unit = row["unit"]
    genome = resolve(root, row["target_genome"], file=True)
    candidate_dir = resolve(root, row["candidate_dir"])
    output_dir = resolve(root, row["output_dir"])
    if output_dir.exists():
        raise SearchError(f"{unit}: refusing to overwrite output")
    candidates_path, queries_path, candidate_manifest = validate_candidate_bundle(
        unit, candidate_dir, genome
    )
    candidate_rows = read_candidates(candidates_path)
    query_ids = read_fasta_ids(queries_path)
    expected_query_ids = {
        row["reference_gene"] for row in candidate_rows if row["has_reference_cds"] == "true"
    }
    if query_ids != expected_query_ids:
        raise SearchError(f"{unit}: candidate table/query FASTA ID closure failed")
    subject_ids = read_fasta_ids(genome)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.staging.", dir=output_dir.parent))
    try:
        db_dir = staging / "_blastdb"
        db_dir.mkdir()
        db_prefix = db_dir / "genome"
        build_database(genome, db_prefix, tools["makeblastdb"], staging / "makeblastdb.stderr.log")
        db_ids_path = staging / "blastdb_ids.txt"
        run_command(
            [str(tools["blastdbcmd"]), "-db", str(db_prefix), "-entry", "all", "-outfmt", "%i"],
            stdout=db_ids_path,
            stderr=staging / "blastdbcmd.stderr.log",
        )
        db_ids = {line.strip().removeprefix("lcl|") for line in db_ids_path.read_text(encoding="utf-8").splitlines() if line.strip()}
        if db_ids != subject_ids:
            raise SearchError(f"{unit}: BLAST database sequence-ID closure failed")

        raw_hits = staging / "raw_hits.tsv"
        run_command(
            [
                str(tools["tblastx"]), "-query", str(queries_path), "-db", str(db_prefix),
                "-evalue", f"{maximum_evalue:.12g}", "-max_hsps", str(max_hsps),
                "-max_target_seqs", str(len(subject_ids)), "-num_threads", str(threads),
                "-outfmt", "6 " + " ".join(HIT_COLUMNS), "-out", str(raw_hits),
            ],
            stderr=staging / "tblastx.stderr.log",
        )
        hits = parse_hits(raw_hits, query_ids, subject_ids)
        state_rows, counts = classify(
            unit=unit,
            candidate_rows=candidate_rows,
            hits=hits,
            minimum_identity=minimum_identity,
            minimum_bitscore=minimum_bitscore,
            maximum_evalue=maximum_evalue,
        )
        state_path = staging / "loss_states.tsv"
        with state_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=STATE_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(state_rows)
        positives = staging / "positive_deleted_reference_genes.txt"
        positives.write_text(
            "".join(f"{row['reference_gene']}\n" for row in state_rows if row["positive_loss"] == "true"),
            encoding="utf-8",
        )
        manifest = {
            "schema_version": 1,
            "workflow": "callable_aware_translated_genome_search",
            "status": "PASS",
            "unit": unit,
            "finished_at_utc": utc_now(),
            "parameters": {
                "tblastx_evalue_command": maximum_evalue,
                "postfilter_evalue_operator": "<",
                "minimum_percent_identity": minimum_identity,
                "minimum_bitscore": minimum_bitscore,
                "minimum_alignment_length": None,
                "max_hsps_per_subject_sequence": max_hsps,
                "max_target_sequences": len(subject_ids),
                "threads": threads,
                "primary_positive_rule": "callable bilateral SynOrths locus with no qualifying local translated hit",
                "local_sequence_hit_rule": "uncertain, never pseudogenized without disruptive-mutation evidence",
                "historical_reproduction_rule": "best genome-wide hit passing identity/evalue/bitscore is decayed; otherwise deleted",
            },
            "tools": {role: binding(path) for role, path in tools.items()},
            "inputs": {
                "target_genome": binding(genome),
                "candidate_manifest": binding(candidate_dir / "run_manifest.json"),
                "candidates": binding(candidates_path),
                "candidate_reference_cds": binding(queries_path),
            },
            "candidate_manifest_metrics": candidate_manifest.get("metrics"),
            "metrics": {
                "candidate_rows": len(candidate_rows),
                "query_cds_records": len(query_ids),
                "subject_sequences": len(subject_ids),
                "raw_hit_rows": sum(len(items) for items in hits.values()),
                "state_counts": counts,
                "positive_deleted": sum(row["positive_loss"] == "true" for row in state_rows),
            },
        }
        manifest_path = staging / "run_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        shutil.rmtree(db_dir)
        checksums = staging / "checksums.tsv"
        with checksums.open("w", encoding="utf-8", newline="") as handle:
            handle.write("file\tbytes\tsha256\n")
            for path in sorted(staging.iterdir()):
                if path.is_file() and path.name != checksums.name:
                    item = binding(path, allow_empty=True)
                    handle.write(f"{path.name}\t{item['bytes']}\t{item['sha256']}\n")
        os.rename(staging, output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return {
        "unit": unit,
        "finished_at_utc": utc_now(),
        "manifest_sha256": sha256(output_dir / "run_manifest.json"),
        "positive_deleted": manifest["metrics"]["positive_deleted"],
    }


def write_state(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--blast-bin", required=True, type=Path)
    parser.add_argument("--queue-root", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=14)
    parser.add_argument("--minimum-identity", type=float, default=50.0)
    parser.add_argument("--minimum-bitscore", type=float, default=50.0)
    parser.add_argument("--maximum-evalue", type=float, default=1e-5)
    parser.add_argument("--max-hsps", type=int, default=20)
    parser.add_argument("--poll-seconds", type=int, default=30)
    args = parser.parse_args()
    lock_handle = None
    try:
        if not 1 <= args.threads <= 15 or args.max_hsps < 1 or args.poll_seconds < 1:
            raise SearchError("invalid worker/search parameter")
        root = args.data_root.resolve()
        rows = read_manifest(args.manifest)
        tools = {name: (args.blast_bin / name).resolve() for name in ("makeblastdb", "blastdbcmd", "tblastx")}
        for name, path in tools.items():
            if not path.is_file() or not os.access(path, os.X_OK):
                raise SearchError(f"missing executable {name}")
        queue_root = args.queue_root.resolve()
        queue_root.mkdir(parents=True, exist_ok=True)
        lock_handle = (queue_root / "controller.lock").open("w")
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise SearchError("another translated-search controller owns the queue") from error
        state_path = queue_root / "state.json"
        state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "callable_aware_translated_search_queue",
            "status": "running",
            "started_at_utc": utc_now(),
            "threads": args.threads,
            "pending": [row["unit"] for row in rows],
            "completed": [],
        }
        write_state(state_path, state)
        for row in rows:
            candidate_manifest = resolve(root, row["candidate_dir"]) / "run_manifest.json"
            while not candidate_manifest.is_file():
                state["waiting_for_candidate_unit"] = row["unit"]
                write_state(state_path, state)
                time.sleep(args.poll_seconds)
            state.pop("waiting_for_candidate_unit", None)
            state["active_unit"] = row["unit"]
            write_state(state_path, state)
            result = run_unit(
                row=row,
                root=root,
                tools=tools,
                threads=args.threads,
                minimum_identity=args.minimum_identity,
                minimum_bitscore=args.minimum_bitscore,
                maximum_evalue=args.maximum_evalue,
                max_hsps=args.max_hsps,
            )
            state["completed"].append(result)  # type: ignore[union-attr]
            state["pending"] = [item for item in state["pending"] if item != row["unit"]]  # type: ignore[union-attr]
            state.pop("active_unit", None)
            write_state(state_path, state)
        state["status"] = "PASS"
        state["finished_at_utc"] = utc_now()
        write_state(state_path, state)
        print(f"PASS\t{len(rows)} units")
        return 0
    except (OSError, UnicodeError, csv.Error, json.JSONDecodeError, SearchError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if lock_handle is not None:
            lock_handle.close()


if __name__ == "__main__":
    sys.exit(main())
