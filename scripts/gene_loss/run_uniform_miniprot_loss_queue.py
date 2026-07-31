#!/usr/bin/env python3
"""Classify all old/new loss candidates with one Miniprot evidence standard.

Exact SynOrths anchors are handled later as retained calls.  This queue treats
every non-anchor candidate identically: non-callable intervals remain
uncertain; a callable interval without a qualifying local protein-to-genome
alignment is deleted; a qualifying local alignment containing a Miniprot
frameshift or in-frame stop is pseudogenized; and an otherwise qualifying
local alignment remains uncertain because sequence alone is not retained or
pseudogene proof.
"""

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
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


class UniformLossError(RuntimeError):
    pass


MANIFEST_COLUMNS = ("unit", "target_genome", "candidate_dir", "output_dir")
STATE_COLUMNS = (
    "unit",
    "reference_gene",
    "callable",
    "callability_reason",
    "target_chromosome",
    "target_interval_start_1based",
    "target_interval_end_1based",
    "qualifying_local_alignment",
    "alignment_target_start_1based",
    "alignment_target_end_1based",
    "alignment_strand",
    "query_length_aa",
    "query_aligned_start_0based",
    "query_aligned_end_0based_exclusive",
    "query_coverage",
    "exact_alignment_identity",
    "alignment_score",
    "frameshift_events",
    "inframe_stop_codons",
    "classification",
    "positive_loss",
    "evidence_reason",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.resolve().open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    source = path.resolve()
    if not source.is_file() or (source.stat().st_size == 0 and not allow_empty):
        raise UniformLossError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def resolve(root: Path, value: str, *, file: bool = False) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise UniformLossError(f"unsafe data-root-relative path: {value!r}")
    path = (root / relative).absolute()
    if not path.is_relative_to(root):
        raise UniformLossError(f"path escapes data root: {value!r}")
    if file and (not path.is_file() or path.stat().st_size == 0):
        raise UniformLossError(f"missing input file: {value!r}")
    return path


def strict_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise UniformLossError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise UniformLossError(f"{path}: JSON root is not an object")
    return value


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MANIFEST_COLUMNS:
            raise UniformLossError("manifest columns differ from exact schema")
        rows = list(reader)
    units = [row["unit"] for row in rows]
    if not rows or len(units) != len(set(units)) or any(not value for row in rows for value in row.values()):
        raise UniformLossError("manifest is empty or has empty/duplicate unit rows")
    return rows


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8"
    )


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = ""
    pieces: list[str] = []
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if current:
                    sequence = "".join(pieces).replace("*", "").upper()
                    if not sequence:
                        raise UniformLossError(f"{path.name}: empty record {current!r}")
                    records[current] = sequence
                current = raw[1:].strip().split()[0]
                if not current or current in records:
                    raise UniformLossError(f"{path.name}:{line_number}: empty/duplicate FASTA ID")
                pieces = []
            elif raw.strip():
                if not current:
                    raise UniformLossError(f"{path.name}:{line_number}: sequence before header")
                pieces.append("".join(raw.split()))
    if current:
        sequence = "".join(pieces).replace("*", "").upper()
        if not sequence:
            raise UniformLossError(f"{path.name}: empty record {current!r}")
        records[current] = sequence
    if not records:
        raise UniformLossError(f"{path.name}: no FASTA records")
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for identifier in sorted(records):
            handle.write(f">{identifier}\n")
            sequence = records[identifier]
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def read_candidates(path: Path, unit: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    required = {
        "unit",
        "reference_gene",
        "callable",
        "callability_reason",
        "target_chromosome",
        "target_interval_start_1based",
        "target_interval_end_1based",
    }
    if not rows or not required.issubset(reader.fieldnames or []):
        raise UniformLossError(f"{path.name}: candidate schema/rows invalid")
    genes = [row["reference_gene"] for row in rows]
    if len(genes) != len(set(genes)) or any(row["unit"] != unit for row in rows):
        raise UniformLossError(f"{unit}: candidate unit/gene closure failed")
    return rows


def validate_candidates(unit: str, directory: Path, genome: Path) -> tuple[Path, dict[str, object]]:
    manifest_path = directory / "run_manifest.json"
    candidate_path = directory / "candidates.tsv"
    for path in (manifest_path, candidate_path, directory / "checksums.tsv"):
        if not path.is_file() or path.stat().st_size == 0:
            raise UniformLossError(f"{unit}: candidate bundle incomplete")
    manifest = strict_json(manifest_path)
    if manifest.get("status") != "PASS" or manifest.get("unit") != unit:
        raise UniformLossError(f"{unit}: candidate manifest is not PASS")
    inputs, outputs = manifest.get("inputs"), manifest.get("outputs")
    if not isinstance(inputs, dict) or not isinstance(outputs, dict):
        raise UniformLossError(f"{unit}: candidate bindings missing")
    if inputs.get("target_genome") != binding(genome) or outputs.get("candidates") != binding(candidate_path):
        raise UniformLossError(f"{unit}: candidate/genome binding mismatch")
    return candidate_path, manifest


def integer(value: str, context: str) -> int:
    try:
        result = int(value)
    except ValueError as error:
        raise UniformLossError(f"{context}: invalid integer {value!r}") from error
    return result


def parse_paf_line(raw: str, line_number: int) -> dict[str, object]:
    fields = raw.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise UniformLossError(f"PAF:{line_number}: fewer than 12 fields")
    values = [integer(fields[index], f"PAF:{line_number}:{index + 1}") for index in (1, 2, 3, 6, 7, 8, 9, 10, 11)]
    qlen, qstart, qend, tlen, tstart, tend, matches, aligned, mapq = values
    if qlen < 1 or not 0 <= qstart < qend <= qlen or tlen < 1 or not 0 <= tstart < tend <= tlen:
        raise UniformLossError(f"PAF:{line_number}: invalid query/target coordinates")
    if matches < 0 or aligned < 1 or matches > aligned or fields[4] not in {"+", "-"}:
        raise UniformLossError(f"PAF:{line_number}: invalid alignment metrics")
    tags: dict[str, tuple[str, str]] = {}
    for item in fields[12:]:
        pieces = item.split(":", 2)
        if len(pieces) != 3:
            raise UniformLossError(f"PAF:{line_number}: malformed tag {item!r}")
        tags[pieces[0]] = (pieces[1], pieces[2])
    for required in ("AS", "fs", "st"):
        if required not in tags or tags[required][0] != "i":
            raise UniformLossError(f"PAF:{line_number}: missing integer {required} tag")
    score = integer(tags["AS"][1], f"PAF:{line_number}:AS")
    frameshifts = integer(tags["fs"][1], f"PAF:{line_number}:fs")
    stops = integer(tags["st"][1], f"PAF:{line_number}:st")
    if frameshifts < 0 or stops < 0:
        raise UniformLossError(f"PAF:{line_number}: negative disruption count")
    return {
        "query": fields[0],
        "query_length": qlen,
        "query_start": qstart,
        "query_end": qend,
        "strand": fields[4],
        "target": fields[5],
        "target_length": tlen,
        "target_start_1based": tstart + 1,
        "target_end_1based": tend,
        "query_coverage": (qend - qstart) / qlen,
        "identity": matches / aligned,
        "alignment_score": score,
        "frameshifts": frameshifts,
        "stops": stops,
        "mapq": mapq,
    }


def read_paf(path: Path, query_ids: set[str]) -> dict[str, list[dict[str, object]]]:
    alignments: dict[str, list[dict[str, object]]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            row = parse_paf_line(raw, line_number)
            query = str(row["query"])
            if query not in query_ids:
                raise UniformLossError(f"PAF:{line_number}: query outside callable protein set")
            alignments[query].append(row)
    return alignments


def best_local_alignment(
    candidate: dict[str, str],
    alignments: list[dict[str, object]],
) -> dict[str, object] | None:
    if candidate["callable"] != "true":
        return None
    start = integer(candidate["target_interval_start_1based"], "candidate interval start")
    end = integer(candidate["target_interval_end_1based"], "candidate interval end")
    local = [
        row
        for row in alignments
        if row["target"] == candidate["target_chromosome"]
        and int(row["target_start_1based"]) >= start
        and int(row["target_end_1based"]) <= end
    ]
    if not local:
        return None
    return min(
        local,
        key=lambda row: (
            -int(row["alignment_score"]),
            -float(row["query_coverage"]),
            -float(row["identity"]),
            str(row["target"]),
            int(row["target_start_1based"]),
        ),
    )


def classify_candidate(
    candidate: dict[str, str],
    alignment: dict[str, object] | None,
    *,
    minimum_query_coverage: float,
    minimum_identity: float,
    minimum_alignment_score: int,
) -> dict[str, str]:
    qualifying = alignment is not None and (
        float(alignment["query_coverage"]) >= minimum_query_coverage
        and float(alignment["identity"]) >= minimum_identity
        and int(alignment["alignment_score"]) >= minimum_alignment_score
    )
    if candidate["callable"] != "true":
        classification, positive, reason = "uncertain", False, "non_callable_interval"
    elif not qualifying:
        classification, positive, reason = "deleted", True, "no_qualifying_local_protein_genome_alignment"
    elif int(alignment["frameshifts"]) > 0 or int(alignment["stops"]) > 0:
        classification, positive, reason = "pseudogenized", True, "explicit_frameshift_or_inframe_stop"
    else:
        classification, positive, reason = "uncertain", False, "local_sequence_without_disruptive_event"
    row = {
        "unit": candidate["unit"],
        "reference_gene": candidate["reference_gene"],
        "callable": candidate["callable"],
        "callability_reason": candidate["callability_reason"],
        "target_chromosome": candidate["target_chromosome"],
        "target_interval_start_1based": candidate["target_interval_start_1based"],
        "target_interval_end_1based": candidate["target_interval_end_1based"],
        "qualifying_local_alignment": str(qualifying).lower(),
        "classification": classification,
        "positive_loss": str(positive).lower(),
        "evidence_reason": reason,
    }
    alignment_fields = {
        "alignment_target_start_1based": "target_start_1based",
        "alignment_target_end_1based": "target_end_1based",
        "alignment_strand": "strand",
        "query_length_aa": "query_length",
        "query_aligned_start_0based": "query_start",
        "query_aligned_end_0based_exclusive": "query_end",
        "query_coverage": "query_coverage",
        "exact_alignment_identity": "identity",
        "alignment_score": "alignment_score",
        "frameshift_events": "frameshifts",
        "inframe_stop_codons": "stops",
    }
    for output, source in alignment_fields.items():
        value = "" if alignment is None else alignment[source]
        row[output] = f"{value:.12g}" if isinstance(value, float) else str(value)
    return row


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_miniprot(command: list[str], stdout: Path, stderr: Path) -> None:
    with stdout.open("wb") as out, stderr.open("wb") as err:
        completed = subprocess.run(command, stdout=out, stderr=err, check=False)
    if completed.returncode:
        raise UniformLossError(f"Miniprot failed with exit code {completed.returncode}")


def tool_version(path: Path) -> str:
    completed = subprocess.run([str(path), "--version"], capture_output=True, text=True, check=False)
    value = (completed.stdout or completed.stderr).strip().splitlines()
    if completed.returncode or len(value) != 1:
        raise UniformLossError("cannot obtain exact Miniprot version")
    return value[0]


def run_unit(
    row: dict[str, str],
    *,
    root: Path,
    reference_proteins: dict[str, str],
    reference_protein_path: Path,
    miniprot: Path,
    version: str,
    threads: int,
    minimum_query_coverage: float,
    minimum_identity: float,
    minimum_alignment_score: int,
) -> dict[str, object]:
    unit = row["unit"]
    genome = resolve(root, row["target_genome"], file=True)
    candidate_dir = resolve(root, row["candidate_dir"])
    output = resolve(root, row["output_dir"])
    candidate_path, candidate_manifest = validate_candidates(unit, candidate_dir, genome)
    candidates = read_candidates(candidate_path, unit)
    callable_genes = {item["reference_gene"] for item in candidates if item["callable"] == "true"}
    missing = callable_genes.difference(reference_proteins)
    if missing:
        raise UniformLossError(f"{unit}: {len(missing)} callable genes lack reference protein")
    parameter_record = {
        "minimum_query_coverage": minimum_query_coverage,
        "minimum_exact_alignment_identity": minimum_identity,
        "minimum_alignment_score": minimum_alignment_score,
        "positive_classes": ["deleted", "pseudogenized"],
        "disruptive_events": ["Miniprot fs tag >0", "Miniprot st tag >0"],
        "noncanonical_splice_sites_alone_are_not_positive": True,
        "threads": threads,
        "secondary_chain_ratio": 0.2,
        "secondary_output_score_ratio": 0.2,
        "maximum_alignments_per_query": 200,
    }
    if output.exists():
        manifest = strict_json(output / "run_manifest.json")
        if (
            manifest.get("status") == "PASS"
            and manifest.get("unit") == unit
            and manifest.get("parameters") == parameter_record
            and manifest.get("inputs", {}).get("target_genome") == binding(genome)
            and manifest.get("inputs", {}).get("candidate_manifest") == binding(candidate_dir / "run_manifest.json")
            and manifest.get("tools", {}).get("miniprot") == binding(miniprot)
        ):
            return {"unit": unit, "status": "reused_exact_pass", "metrics": manifest["metrics"]}
        raise UniformLossError(f"{unit}: existing output is not an exact reusable PASS")

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        query_path = staging / "callable_reference_proteins.faa"
        write_fasta(query_path, {gene: reference_proteins[gene] for gene in callable_genes})
        raw_paf = staging / "raw_alignments.paf"
        run_miniprot(
            [
                str(miniprot),
                "-t",
                str(threads),
                "-I",
                "-p",
                "0.2",
                "-N",
                "200",
                "--outn",
                "200",
                "--outs",
                "0.2",
                str(genome),
                str(query_path),
            ],
            raw_paf,
            staging / "miniprot.stderr.log",
        )
        alignments = read_paf(raw_paf, callable_genes)
        states = [
            classify_candidate(
                candidate,
                best_local_alignment(candidate, alignments.get(candidate["reference_gene"], [])),
                minimum_query_coverage=minimum_query_coverage,
                minimum_identity=minimum_identity,
                minimum_alignment_score=minimum_alignment_score,
            )
            for candidate in candidates
        ]
        states.sort(key=lambda item: item["reference_gene"])
        state_path = staging / "uniform_candidate_loss_states.tsv"
        write_tsv(state_path, states)
        with raw_paf.open("rb") as source, gzip.open(staging / "raw_alignments.paf.gz", "wb") as target:
            shutil.copyfileobj(source, target, length=1 << 20)
        raw_paf.unlink()
        counts = Counter(item["classification"] for item in states)
        manifest = {
            "schema_version": 1,
            "workflow": "uniform_old_new_miniprot_disruptive_loss_classification",
            "status": "PASS",
            "unit": unit,
            "finished_at_utc": utc_now(),
            "parameters": parameter_record,
            "tools": {"miniprot": binding(miniprot), "miniprot_version": version},
            "inputs": {
                "target_genome": binding(genome),
                "candidate_manifest": binding(candidate_dir / "run_manifest.json"),
                "candidates": binding(candidate_path),
                "reference_protein": binding(reference_protein_path),
            },
            "candidate_manifest_metrics": candidate_manifest.get("metrics"),
            "metrics": {
                "candidate_rows": len(candidates),
                "callable_rows": len(callable_genes),
                "queries_with_any_alignment": len(alignments),
                "classification_counts": dict(sorted(counts.items())),
                "positive_deleted": counts["deleted"],
                "positive_pseudogenized": counts["pseudogenized"],
            },
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        with (staging / "checksums.tsv").open("w", encoding="utf-8", newline="") as handle:
            handle.write("file\tbytes\tsha256\n")
            for path in sorted(staging.iterdir()):
                if path.is_file() and path.name != "checksums.tsv":
                    item = binding(path, allow_empty=True)
                    handle.write(f"{path.name}\t{item['bytes']}\t{item['sha256']}\n")
        os.rename(staging, output)
        return {"unit": unit, "status": "completed", "metrics": manifest["metrics"]}
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def write_state(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--reference-protein", required=True, type=Path)
    parser.add_argument("--miniprot", required=True, type=Path)
    parser.add_argument("--queue-root", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=14)
    parser.add_argument("--minimum-query-coverage", type=float, default=0.50)
    parser.add_argument("--minimum-identity", type=float, default=0.50)
    parser.add_argument("--minimum-alignment-score", type=int, default=50)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--stop-after", type=int, default=0, help="0 runs all units; positive values support audited pilots")
    args = parser.parse_args()
    lock_handle = None
    try:
        if not 1 <= args.threads <= 15 or not 0 < args.minimum_query_coverage <= 1 or not 0 < args.minimum_identity <= 1:
            raise UniformLossError("invalid worker/alignment thresholds")
        if args.minimum_alignment_score < 0 or args.poll_seconds < 1 or args.stop_after < 0:
            raise UniformLossError("invalid score/poll/stop parameter")
        root = args.data_root.resolve()
        rows = read_manifest(args.manifest)
        reference_protein = args.reference_protein.resolve()
        miniprot = args.miniprot.resolve()
        if not reference_protein.is_file() or not miniprot.is_file() or not os.access(miniprot, os.X_OK):
            raise UniformLossError("reference protein or Miniprot executable is unavailable")
        proteins = read_fasta(reference_protein)
        version = tool_version(miniprot)
        queue = args.queue_root.resolve()
        queue.mkdir(parents=True, exist_ok=True)
        lock_handle = (queue / "controller.lock").open("w")
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise UniformLossError("another Miniprot loss controller owns the queue") from error
        state: dict[str, object] = {
            "schema_version": 1,
            "workflow": "uniform_old_new_miniprot_loss_queue",
            "status": "running",
            "started_at_utc": utc_now(),
            "threads": args.threads,
            "pending": [row["unit"] for row in rows],
            "completed": [],
        }
        write_state(queue / "state.json", state)
        processed = 0
        for row in rows:
            candidate_manifest = resolve(root, row["candidate_dir"]) / "run_manifest.json"
            while not candidate_manifest.is_file():
                state["waiting_for_candidate_unit"] = row["unit"]
                write_state(queue / "state.json", state)
                time.sleep(args.poll_seconds)
            state.pop("waiting_for_candidate_unit", None)
            state["active_unit"] = row["unit"]
            write_state(queue / "state.json", state)
            result = run_unit(
                row,
                root=root,
                reference_proteins=proteins,
                reference_protein_path=reference_protein,
                miniprot=miniprot,
                version=version,
                threads=args.threads,
                minimum_query_coverage=args.minimum_query_coverage,
                minimum_identity=args.minimum_identity,
                minimum_alignment_score=args.minimum_alignment_score,
            )
            state["completed"].append(result)  # type: ignore[union-attr]
            state["pending"] = [item for item in state["pending"] if item != row["unit"]]  # type: ignore[union-attr]
            state.pop("active_unit", None)
            write_state(queue / "state.json", state)
            processed += 1
            if args.stop_after and processed >= args.stop_after:
                state["status"] = "PASS_PILOT_STOP"
                state["finished_at_utc"] = utc_now()
                write_state(queue / "state.json", state)
                print(f"PASS_PILOT_STOP\t{processed} units")
                return 0
        state["status"] = "PASS"
        state["finished_at_utc"] = utc_now()
        write_state(queue / "state.json", state)
        print(f"PASS\t{processed} units")
        return 0
    except (OSError, UnicodeError, csv.Error, json.JSONDecodeError, UniformLossError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if lock_handle is not None:
            lock_handle.close()


if __name__ == "__main__":
    sys.exit(main())
