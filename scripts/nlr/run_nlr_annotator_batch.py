#!/usr/bin/env python3
"""Run or safely resume NLR-Annotator with a hard worker-thread ceiling.

Without ``--execute`` this command only validates paths and prints the planned
commands.  Execution is deliberately sequential across manifest rows.  The
NLR-Annotator worker pool is capped at eight threads, and JVM processor/GC
settings are constrained so that this workflow can coexist with other server
jobs.  An interrupted controller can be continued only from an explicitly
named hidden temporary root.  Completed samples are revalidated before reuse;
ambiguous or corrupt state fails closed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


MAX_WORKER_THREADS = 8
SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
MEMORY = re.compile(r"^[1-9][0-9]*[mMgG]$")
SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
SELECTED_INPUT_FIELDS = [
    "sample_id", "species", "ploidy", "analysis_role", "input_scope",
    "relative_fasta", "expected_fasta_records",
]
CHECKSUM_OUTPUTS = ["nlr_calls.txt", "nlr_loci.gff", "stdout.log", "stderr.log"]


@dataclass(frozen=True)
class InputRow:
    sample_id: str
    species: str
    ploidy: str
    analysis_role: str
    input_scope: str
    fasta: Path
    relative_fasta: str
    expected_records: int | None


@dataclass(frozen=True)
class FastaAudit:
    sha256: str
    records: int
    total_bases: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--jar", type=Path, required=True)
    parser.add_argument("--motifs", type=Path, required=True)
    parser.add_argument("--store", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--java-bin", default="java")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--sequences-per-thread", type=int, default=1000)
    parser.add_argument("--max-heap", default="8G")
    parser.add_argument(
        "--role", choices=["all", "reference_callable", "target_repertoire"], default="all"
    )
    parser.add_argument("--sample-id", action="append", default=[],
                        help="Restrict to a sample ID; repeat to select more than one")
    parser.add_argument("--expected-jar-sha256")
    parser.add_argument("--expected-motifs-sha256")
    parser.add_argument("--expected-store-sha256")
    parser.add_argument("--execute", action="store_true",
                        help="Actually launch Java; otherwise print a no-write plan")
    parser.add_argument(
        "--resume-temp-root", type=Path,
        help=(
            "Continue an interrupted hidden .OUTPUT.tmp.* root. The root must be a "
            "direct sibling of --output-root and its selected_inputs.tsv is authoritative"
        ),
    )
    parser.add_argument(
        "--cleanup-partial", action="store_true",
        help=(
            "On resume, discard and rerun every unambiguously partial sample directory. "
            "Use only after confirming that no old controller/Java process is writing there"
        ),
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty TSV: {path}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def read_manifest(path: Path, input_root: Path) -> list[InputRow]:
    fields, rows = read_tsv(path)
    required = [
        "sample_id", "species", "ploidy", "analysis_role", "input_scope",
        "relative_fasta", "expected_fasta_records",
    ]
    missing = [field for field in required if field not in fields]
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
    parsed: list[InputRow] = []
    seen: set[str] = set()
    reference_rows = 0
    for line_no, row in enumerate(rows, 2):
        sample = row["sample_id"]
        if not sample or not SAFE_ID.fullmatch(sample):
            raise ValueError(f"Unsafe or empty sample_id at {path}:{line_no}: {sample!r}")
        if sample in seen:
            raise ValueError(f"Duplicate sample_id {sample!r} in {path}")
        seen.add(sample)
        role = row["analysis_role"]
        scope = row["input_scope"]
        if role == "reference_callable":
            reference_rows += 1
            if scope != "reference_transcript_cds":
                raise ValueError(f"Reference input must use reference_transcript_cds at {path}:{line_no}")
        elif role == "target_repertoire":
            if scope != "whole_genome":
                raise ValueError(f"Target repertoire input must use whole_genome at {path}:{line_no}")
        else:
            raise ValueError(f"Unknown analysis_role {role!r} at {path}:{line_no}")
        relative = Path(row["relative_fasta"])
        if not row["relative_fasta"] or relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"relative_fasta must be safe and relative at {path}:{line_no}")
        expected: int | None = None
        if row["expected_fasta_records"]:
            try:
                expected = int(row["expected_fasta_records"])
            except ValueError as exc:
                raise ValueError(f"Invalid expected_fasta_records at {path}:{line_no}") from exc
            if expected < 1:
                raise ValueError(f"expected_fasta_records must be positive at {path}:{line_no}")
        parsed.append(InputRow(
            sample_id=sample,
            species=row["species"],
            ploidy=row["ploidy"],
            analysis_role=role,
            input_scope=scope,
            fasta=input_root / relative,
            relative_fasta=row["relative_fasta"],
            expected_records=expected,
        ))
    if reference_rows != 1:
        raise ValueError(f"Manifest must contain exactly one reference_callable row; observed {reference_rows}")
    if not parsed:
        raise ValueError(f"No rows found in {path}")
    return parsed


def audit_fasta(path: Path) -> FastaAudit:
    digest = hashlib.sha256()
    records = 0
    total_bases = 0
    seen: set[str] = set()
    current_id: str | None = None
    current_bases = 0
    with path.open("rb") as handle:
        for line_no, raw in enumerate(handle, 1):
            digest.update(raw)
            if raw.startswith(b">"):
                if current_id is not None and current_bases == 0:
                    raise ValueError(f"Empty FASTA record {current_id!r} in {path}")
                try:
                    header = raw[1:].strip().decode("utf-8")
                except UnicodeDecodeError as exc:
                    raise ValueError(f"Non-UTF-8 FASTA header at {path}:{line_no}") from exc
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_no}")
                current_id = header.split()[0]
                if current_id in seen:
                    raise ValueError(f"Duplicate FASTA ID {current_id!r} in {path}")
                seen.add(current_id)
                records += 1
                current_bases = 0
            else:
                if current_id is None:
                    if raw.strip():
                        raise ValueError(f"Sequence data precedes the first header at {path}:{line_no}")
                    continue
                bases = len(b"".join(raw.split()))
                current_bases += bases
                total_bases += bases
    if current_id is not None and current_bases == 0:
        raise ValueError(f"Empty FASTA record {current_id!r} in {path}")
    if records == 0:
        raise ValueError(f"No FASTA records found in {path}")
    return FastaAudit(digest.hexdigest(), records, total_bases)


def build_command(
    args: argparse.Namespace,
    row: InputRow,
    sample_dir: Path,
) -> list[str]:
    java_tmp = sample_dir / "java_tmp"
    return [
        args.java_bin,
        f"-Xmx{args.max_heap}",
        f"-XX:ActiveProcessorCount={args.threads}",
        "-XX:+UseSerialGC",
        f"-Djava.io.tmpdir={java_tmp}",
        "-jar", str(args.jar),
        "-i", str(row.fasta),
        "-x", str(args.motifs),
        "-y", str(args.store),
        "-o", str(sample_dir / "nlr_calls.txt"),
        "-g", str(sample_dir / "nlr_loci.gff"),
        "-t", str(args.threads),
        "-n", str(args.sequences_per_thread),
    ]


def validate_nlr_output(path: Path) -> tuple[int, int, int]:
    if not path.is_file():
        raise ValueError(f"NLR-Annotator did not create {path}")
    rows = 0
    sequence_ids: set[str] = set()
    locus_ids: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_no, values in enumerate(reader, 1):
            if not values:
                continue
            if len(values) != 7:
                raise ValueError(f"Expected seven NLR output columns at {path}:{line_no}")
            if not values[0].strip() or not values[1].strip() or not values[2].strip():
                raise ValueError(f"Empty required NLR field at {path}:{line_no}")
            locus = values[1].strip()
            if locus in locus_ids:
                raise ValueError(f"Duplicate NLR locus ID {locus!r} in {path}")
            rows += 1
            sequence_ids.add(values[0].strip())
            locus_ids.add(locus)
    return rows, len(sequence_ids), len(locus_ids)


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def checked_tool_hash(path: Path, expected: str | None, label: str) -> str:
    observed = sha256_file(path)
    if expected:
        if not SHA256.fullmatch(expected):
            raise ValueError(f"Invalid expected SHA-256 for {label}: {expected!r}")
        if observed.lower() != expected.lower():
            raise ValueError(
                f"{label} checksum mismatch: observed {observed}, expected {expected.lower()}"
            )
    return observed


def input_row_record(row: InputRow) -> dict[str, object]:
    return {
        "sample_id": row.sample_id,
        "species": row.species,
        "ploidy": row.ploidy,
        "analysis_role": row.analysis_role,
        "input_scope": row.input_scope,
        "relative_fasta": row.relative_fasta,
        "expected_fasta_records": row.expected_records or "",
    }


def read_key_value(path: Path) -> dict[str, str]:
    fields, rows = read_tsv(path)
    if fields != ["key", "value"]:
        raise ValueError(f"Unexpected key/value header in {path}: {fields}")
    values: dict[str, str] = {}
    for line_no, row in enumerate(rows, 2):
        key = row["key"]
        if not key or key in values:
            raise ValueError(f"Empty or duplicate metadata key at {path}:{line_no}: {key!r}")
        values[key] = row["value"]
    return values


def validate_resume_root_path(resume_root: Path, output_root: Path) -> Path:
    if output_root.exists():
        raise ValueError(f"Output root already exists: {output_root}")
    if not resume_root.exists() or not resume_root.is_dir() or resume_root.is_symlink():
        raise ValueError(f"Resume root must be an existing non-symlink directory: {resume_root}")
    expected_prefix = f".{output_root.name}.tmp."
    if not resume_root.name.startswith(expected_prefix):
        raise ValueError(
            f"Resume root name must start with {expected_prefix!r}: {resume_root.name!r}"
        )
    if resume_root.parent.resolve() != output_root.parent.resolve():
        raise ValueError("Resume root must be a direct sibling of --output-root")
    return resume_root


def read_resume_snapshot(
    tmp_root: Path,
    manifest_rows: list[InputRow],
    role: str,
) -> list[InputRow]:
    snapshot = tmp_root / "selected_inputs.tsv"
    fields, records = read_tsv(snapshot)
    if fields != SELECTED_INPUT_FIELDS:
        raise ValueError(f"Unexpected selected-input header in {snapshot}: {fields}")
    current = {row.sample_id: row for row in manifest_rows}
    allowed = {
        row.sample_id for row in manifest_rows
        if role == "all" or row.analysis_role == role
    }
    selected: list[InputRow] = []
    seen: set[str] = set()
    for line_no, record in enumerate(records, 2):
        sample = record["sample_id"]
        if sample in seen:
            raise ValueError(f"Duplicate sample ID in resume snapshot at {snapshot}:{line_no}")
        seen.add(sample)
        if sample not in current:
            raise ValueError(f"Resume snapshot sample is absent from current manifest: {sample}")
        if sample not in allowed:
            raise ValueError(f"Resume snapshot sample is excluded by --role={role}: {sample}")
        row = current[sample]
        expected = {key: str(value) for key, value in input_row_record(row).items()}
        if record != expected:
            differences = [
                f"{key}: snapshot={record.get(key)!r}, current={expected.get(key)!r}"
                for key in SELECTED_INPUT_FIELDS if record.get(key) != expected.get(key)
            ]
            raise ValueError(
                f"Resume snapshot no longer matches the manifest for {sample}: "
                + "; ".join(differences)
            )
        selected.append(row)
    if not selected:
        raise ValueError(f"No rows in resume snapshot: {snapshot}")
    return selected


def metadata_int(metadata: dict[str, str], key: str, sample: str) -> int:
    try:
        return int(metadata[key])
    except KeyError as exc:
        raise ValueError(f"Completed sample {sample} lacks metadata key {key}") from exc
    except ValueError as exc:
        raise ValueError(f"Completed sample {sample} has non-integer metadata {key}") from exc


def validate_completed_sample(
    args: argparse.Namespace,
    row: InputRow,
    sample_dir: Path,
    tool_hashes: dict[str, str],
    *,
    expected_input_fasta_text: str | None = None,
    input_fasta_audit_path: Path | None = None,
) -> int:
    """Revalidate one completed directory and return its recorded worker count.

    ``expected_input_fasta_text`` preserves the immutable lexical path written
    by an older runner, while ``input_fasta_audit_path`` can independently
    select the already identity-checked current file for content auditing.
    Normal runner/resume validation leaves both unset.
    """
    if sample_dir.is_symlink() or not sample_dir.is_dir():
        raise ValueError(f"Completed sample path is not a normal directory: {sample_dir}")
    required = set(CHECKSUM_OUTPUTS + ["run_metadata.tsv", "output_checksums.tsv"])
    observed = {entry.name for entry in sample_dir.iterdir()}
    if observed != required:
        raise ValueError(
            f"Completed sample {row.sample_id} has unexpected directory contents: "
            f"missing={sorted(required - observed)}, extra={sorted(observed - required)}"
        )
    for name in required:
        path = sample_dir / name
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Completed sample file is not a regular non-symlink file: {path}")

    metadata = read_key_value(sample_dir / "run_metadata.tsv")
    fixed_expected = {
        "sample_id": row.sample_id,
        "species": row.species,
        "ploidy": row.ploidy,
        "analysis_role": row.analysis_role,
        "input_scope": row.input_scope,
        "input_fasta": (
            expected_input_fasta_text
            if expected_input_fasta_text is not None
            else str(row.fasta)
        ),
        "sequences_per_thread": str(args.sequences_per_thread),
        "java_max_heap": args.max_heap,
        "maximum_allowed_nlr_worker_threads": str(MAX_WORKER_THREADS),
        "jvm_gc": "UseSerialGC",
        "nlr_annotator_jar_sha256": tool_hashes["jar"],
        "motifs_sha256": tool_hashes["motifs"],
        "store_sha256": tool_hashes["store"],
        "completion_status": "complete",
    }
    for key, expected in fixed_expected.items():
        if metadata.get(key) != expected:
            raise ValueError(
                f"Completed sample {row.sample_id} metadata mismatch for {key}: "
                f"observed={metadata.get(key)!r}, expected={expected!r}"
            )
    for key in ["executed_command", "replay_command_at_final_output_path"]:
        if not metadata.get(key):
            raise ValueError(f"Completed sample {row.sample_id} lacks {key}")

    fasta_audit = audit_fasta(input_fasta_audit_path or row.fasta)
    if row.expected_records is not None and fasta_audit.records != row.expected_records:
        raise ValueError(
            f"{row.sample_id} FASTA has {fasta_audit.records} records, "
            f"expected {row.expected_records}"
        )
    audit_expected = {
        "input_fasta_sha256": fasta_audit.sha256,
        "input_fasta_records": str(fasta_audit.records),
        "input_fasta_total_bases": str(fasta_audit.total_bases),
    }
    for key, expected in audit_expected.items():
        if metadata.get(key) != expected:
            raise ValueError(
                f"Completed sample {row.sample_id} input audit mismatch for {key}: "
                f"observed={metadata.get(key)!r}, current={expected!r}"
            )

    output_rows, sequence_ids, locus_ids = validate_nlr_output(sample_dir / "nlr_calls.txt")
    output_expected = {
        "nlr_output_rows": output_rows,
        "nlr_output_sequence_ids": sequence_ids,
        "nlr_output_locus_ids": locus_ids,
    }
    for key, expected in output_expected.items():
        if metadata_int(metadata, key, row.sample_id) != expected:
            raise ValueError(f"Completed sample {row.sample_id} output audit mismatch for {key}")

    checksum_fields, checksum_rows = read_tsv(sample_dir / "output_checksums.tsv")
    if checksum_fields != ["path", "sha256"]:
        raise ValueError(f"Unexpected checksum header for completed sample {row.sample_id}")
    checksums: dict[str, str] = {}
    for record in checksum_rows:
        name, checksum = record["path"], record["sha256"].lower()
        if name in checksums or name not in CHECKSUM_OUTPUTS or not SHA256.fullmatch(checksum):
            raise ValueError(f"Invalid checksum row for completed sample {row.sample_id}: {record}")
        checksums[name] = checksum
    if set(checksums) != set(CHECKSUM_OUTPUTS):
        raise ValueError(f"Incomplete checksum inventory for completed sample {row.sample_id}")
    for name, expected in checksums.items():
        observed_hash = sha256_file(sample_dir / name)
        if observed_hash != expected:
            raise ValueError(
                f"Completed sample {row.sample_id} checksum mismatch for {name}: "
                f"observed={observed_hash}, recorded={expected}"
            )

    workers = metadata_int(metadata, "configured_nlr_worker_threads", row.sample_id)
    processor_cap = metadata_int(metadata, "jvm_processor_cap", row.sample_id)
    if not 1 <= workers <= MAX_WORKER_THREADS or processor_cap != workers:
        raise ValueError(f"Invalid recorded worker/JVM cap for completed sample {row.sample_id}")
    return workers


def inspect_resume_root(
    args: argparse.Namespace,
    tmp_root: Path,
    selected: list[InputRow],
    tool_hashes: dict[str, str],
) -> tuple[list[InputRow], list[InputRow], list[InputRow], dict[str, int]]:
    expected_ids = {row.sample_id for row in selected}
    allowed_files = {"selected_inputs.tsv", "resume_history.tsv", "batch_metadata.tsv"}
    for entry in tmp_root.iterdir():
        if entry.name in allowed_files or entry.name in expected_ids:
            continue
        raise ValueError(f"Unexpected entry in resume root: {entry}")

    completed: list[InputRow] = []
    missing: list[InputRow] = []
    partial: list[InputRow] = []
    worker_counts: dict[str, int] = {}
    for row in selected:
        sample_dir = tmp_root / row.sample_id
        if not sample_dir.exists():
            missing.append(row)
            continue
        if sample_dir.is_symlink() or not sample_dir.is_dir():
            raise ValueError(f"Sample state path is not a normal directory: {sample_dir}")
        has_metadata = (sample_dir / "run_metadata.tsv").exists()
        has_checksums = (sample_dir / "output_checksums.tsv").exists()
        if has_metadata and has_checksums:
            try:
                worker_counts[row.sample_id] = validate_completed_sample(
                    args, row, sample_dir, tool_hashes
                )
            except (OSError, UnicodeError, ValueError) as exc:
                raise ValueError(
                    f"Completed-looking sample {row.sample_id} failed revalidation and "
                    f"will not be cleaned automatically: {exc}"
                ) from exc
            completed.append(row)
        else:
            partial.append(row)
    batch_path = tmp_root / "batch_metadata.tsv"
    if batch_path.exists():
        if missing or partial:
            raise ValueError(
                "Resume root contains batch_metadata.tsv marked as published-ready but "
                "also has missing or partial samples"
            )
        batch = read_key_value(batch_path)
        batch_expected = {
            "selected_inputs": str(len(selected)),
            "maximum_allowed_nlr_worker_threads_per_process": str(MAX_WORKER_THREADS),
            "nlr_annotator_jar_sha256": tool_hashes["jar"],
            "motifs_sha256": tool_hashes["motifs"],
            "store_sha256": tool_hashes["store"],
            "completion_status": "complete",
        }
        for key, expected in batch_expected.items():
            if batch.get(key) != expected:
                raise ValueError(
                    f"Existing batch metadata mismatch for {key}: "
                    f"observed={batch.get(key)!r}, expected={expected!r}"
                )
    return completed, missing, partial, worker_counts


def run_one_sample(
    args: argparse.Namespace,
    row: InputRow,
    sample_dir: Path,
    tool_hashes: dict[str, str],
    index: int,
    total: int,
) -> None:
    sample_dir.mkdir()
    (sample_dir / "java_tmp").mkdir()
    fasta_audit = audit_fasta(row.fasta)
    if row.expected_records is not None and fasta_audit.records != row.expected_records:
        raise ValueError(
            f"{row.sample_id} FASTA has {fasta_audit.records} records, "
            f"expected {row.expected_records}"
        )
    command = build_command(args, row, sample_dir)
    print(
        f"[{index}/{total}] {row.sample_id}: launching one Java process "
        f"with {args.threads} NLR worker threads",
        flush=True,
    )
    environment = os.environ.copy()
    for variable in ["_JAVA_OPTIONS", "JAVA_TOOL_OPTIONS", "JDK_JAVA_OPTIONS"]:
        environment.pop(variable, None)
    for variable in [
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
    ]:
        environment[variable] = "1"
    with (sample_dir / "stdout.log").open("w", encoding="utf-8") as stdout, \
            (sample_dir / "stderr.log").open("w", encoding="utf-8") as stderr:
        completed = subprocess.run(command, stdout=stdout, stderr=stderr, env=environment, check=False)
    if completed.returncode != 0:
        raise ValueError(
            f"NLR-Annotator failed for {row.sample_id} with exit code "
            f"{completed.returncode}; inspect {sample_dir / 'stderr.log'}"
        )
    call_rows, sequence_ids, locus_ids = validate_nlr_output(sample_dir / "nlr_calls.txt")
    if not (sample_dir / "nlr_loci.gff").is_file():
        raise ValueError(f"NLR-Annotator did not create {sample_dir / 'nlr_loci.gff'}")
    shutil.rmtree(sample_dir / "java_tmp", ignore_errors=True)
    replay_command = build_command(args, row, args.output_root / row.sample_id)
    metadata = [
        ("sample_id", row.sample_id), ("species", row.species), ("ploidy", row.ploidy),
        ("analysis_role", row.analysis_role), ("input_scope", row.input_scope),
        ("input_fasta", str(row.fasta)), ("input_fasta_sha256", fasta_audit.sha256),
        ("input_fasta_records", fasta_audit.records),
        ("input_fasta_total_bases", fasta_audit.total_bases),
        ("nlr_output_rows", call_rows), ("nlr_output_sequence_ids", sequence_ids),
        ("nlr_output_locus_ids", locus_ids),
        ("configured_nlr_worker_threads", args.threads),
        ("maximum_allowed_nlr_worker_threads", MAX_WORKER_THREADS),
        ("sequences_per_thread", args.sequences_per_thread), ("java_max_heap", args.max_heap),
        ("jvm_processor_cap", args.threads), ("jvm_gc", "UseSerialGC"),
        ("nlr_annotator_jar_sha256", tool_hashes["jar"]),
        ("motifs_sha256", tool_hashes["motifs"]), ("store_sha256", tool_hashes["store"]),
        ("executed_command", shlex.join(command)),
        ("replay_command_at_final_output_path", shlex.join(replay_command)),
        ("completion_status", "complete"),
    ]
    write_tsv(
        sample_dir / "run_metadata.tsv", ["key", "value"],
        ({"key": key, "value": value} for key, value in metadata),
    )
    write_tsv(
        sample_dir / "output_checksums.tsv", ["path", "sha256"],
        ({"path": name, "sha256": sha256_file(sample_dir / name)} for name in CHECKSUM_OUTPUTS),
    )


def append_resume_history(path: Path, record: dict[str, object]) -> None:
    fields = [
        "timestamp_utc", "threads", "requested_sample_ids", "cleanup_partial",
        "completed_before", "missing_before", "partial_before", "discarded_partial",
        "samples_run", "outcome",
    ]
    previous: list[dict[str, str]] = []
    if path.exists():
        observed_fields, previous = read_tsv(path)
        if observed_fields != fields:
            raise ValueError(f"Unexpected resume-history header: {path}")
    temporary = path.with_name(path.name + ".new")
    write_tsv(temporary, fields, [*previous, record])
    os.replace(temporary, path)


def main() -> None:
    args = parse_args()
    if not 1 <= args.threads <= MAX_WORKER_THREADS:
        raise SystemExit(
            f"ERROR: --threads must be between 1 and {MAX_WORKER_THREADS}; "
            "the project contract stays below ten NLR worker threads"
        )
    if args.sequences_per_thread < 1:
        raise SystemExit("ERROR: --sequences-per-thread must be positive")
    if not MEMORY.fullmatch(args.max_heap):
        raise SystemExit("ERROR: --max-heap must look like 8000M or 8G")
    if not args.manifest.is_file() or args.manifest.stat().st_size == 0:
        raise SystemExit(f"ERROR: Missing or empty manifest: {args.manifest}")
    if not args.input_root.is_dir():
        raise SystemExit(f"ERROR: Missing input root: {args.input_root}")
    for path in [args.jar, args.motifs, args.store]:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: Missing or empty NLR-Annotator input: {path}")
    java_path = shutil.which(args.java_bin) if os.sep not in args.java_bin else args.java_bin
    if not java_path or not os.access(java_path, os.X_OK):
        raise SystemExit(f"ERROR: Java executable not found or not executable: {args.java_bin}")
    args.java_bin = str(java_path)

    try:
        rows = read_manifest(args.manifest, args.input_root)
        role_selected = [row for row in rows if args.role == "all" or row.analysis_role == args.role]
        requested = set(args.sample_id)
        if len(requested) != len(args.sample_id):
            raise ValueError("Duplicate --sample-id values are not allowed")
        tool_hashes = {
            "jar": checked_tool_hash(args.jar, args.expected_jar_sha256, "NLR-Annotator JAR"),
            "motifs": checked_tool_hash(args.motifs, args.expected_motifs_sha256, "motifs file"),
            "store": checked_tool_hash(args.store, args.expected_store_sha256, "store file"),
        }
        if args.resume_temp_root:
            tmp_root = validate_resume_root_path(args.resume_temp_root, args.output_root)
            selected = read_resume_snapshot(tmp_root, rows, args.role)
            unknown = requested - {row.sample_id for row in selected}
            if unknown:
                raise ValueError(
                    "Requested samples are absent from the frozen resume snapshot: "
                    + ", ".join(sorted(unknown))
                )
        else:
            if args.cleanup_partial:
                raise ValueError("--cleanup-partial requires --resume-temp-root")
            unknown = requested - {row.sample_id for row in role_selected}
            if unknown:
                raise ValueError(f"Unknown requested sample IDs: {', '.join(sorted(unknown))}")
            selected = [row for row in role_selected if not requested or row.sample_id in requested]
            tmp_root = None
        if not selected:
            raise ValueError("No manifest rows remain after role/sample filtering")
        for row in selected:
            if not row.fasta.is_file() or row.fasta.stat().st_size == 0:
                raise ValueError(f"Missing or empty FASTA for {row.sample_id}: {row.fasta}")
    except (OSError, UnicodeError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}")

    completed: list[InputRow] = []
    missing = list(selected)
    partial: list[InputRow] = []
    if tmp_root is not None:
        try:
            completed, missing, partial, _ = inspect_resume_root(
                args, tmp_root, selected, tool_hashes
            )
            partial_ids = {row.sample_id for row in partial}
            exact_partial_selection = bool(requested) and requested == partial_ids
            if requested and partial and not exact_partial_selection:
                raise ValueError(
                    "When --sample-id is used with partial state, it must name exactly all "
                    f"partial samples and no others; partial={sorted(partial_ids)}, "
                    f"requested={sorted(requested)}"
                )
            if partial and not (args.cleanup_partial or exact_partial_selection):
                raise ValueError(
                    "Unambiguously partial sample directories were found: "
                    f"{', '.join(row.sample_id for row in partial)}. Confirm that the old "
                    "controller and Java child are stopped, then use --cleanup-partial, or "
                    "use --sample-id once for exactly the listed current sample(s)"
                )
        except (OSError, UnicodeError, ValueError) as exc:
            raise SystemExit(f"ERROR: Resume validation failed closed: {exc}")

    unfinished_ids = {row.sample_id for row in [*missing, *partial]}
    run_candidates = [row for row in selected if row.sample_id in unfinished_ids]
    if requested and args.resume_temp_root:
        run_candidates = [row for row in run_candidates if row.sample_id in requested]

    if not args.execute:
        if args.resume_temp_root:
            print(
                f"RESUME PLAN ONLY: completed={len(completed)}, missing={len(missing)}, "
                f"partial={len(partial)}, selected_to_run={len(run_candidates)}; "
                f"worker threads for new runs={args.threads}"
            )
            if partial:
                print("WOULD DISCARD PARTIAL: " + ",".join(row.sample_id for row in partial))
        else:
            print(
                f"PLAN ONLY: {len(run_candidates)} NLR-Annotator process(es), sequential; "
                f"worker threads per process={args.threads} (hard maximum={MAX_WORKER_THREADS})"
            )
        plan_root = tmp_root if tmp_root is not None else args.output_root
        for row in run_candidates:
            print(shlex.join(build_command(args, row, plan_root / row.sample_id)))
        return

    fresh = tmp_root is None
    if fresh:
        if args.output_root.exists():
            raise SystemExit(f"ERROR: Output root already exists: {args.output_root}")
        output_parent = args.output_root.parent
        output_parent.mkdir(parents=True, exist_ok=True)
        tmp_root = Path(tempfile.mkdtemp(prefix=f".{args.output_root.name}.tmp.", dir=output_parent))
        write_tsv(
            tmp_root / "selected_inputs.tsv", SELECTED_INPUT_FIELDS,
            (input_row_record(row) for row in selected),
        )

    assert tmp_root is not None
    discarded = [row.sample_id for row in partial]
    samples_run: list[str] = []
    try:
        for row in partial:
            shutil.rmtree(tmp_root / row.sample_id)
        selected_positions = {row.sample_id: index for index, row in enumerate(selected, 1)}
        for row in run_candidates:
            run_one_sample(
                args, row, tmp_root / row.sample_id, tool_hashes,
                selected_positions[row.sample_id], len(selected),
            )
            samples_run.append(row.sample_id)

        final_completed, final_missing, final_partial, final_workers = inspect_resume_root(
            args, tmp_root, selected, tool_hashes
        )
        if final_partial:
            raise ValueError("Internal error: partial sample remained after successful execution")
        if final_missing:
            if not args.resume_temp_root:
                raise ValueError("Internal error: a fresh run left missing selected samples")
            append_resume_history(
                tmp_root / "resume_history.tsv",
                {
                    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                    "threads": args.threads,
                    "requested_sample_ids": ",".join(args.sample_id),
                    "cleanup_partial": str(args.cleanup_partial).upper(),
                    "completed_before": len(completed), "missing_before": len(missing),
                    "partial_before": len(partial), "discarded_partial": ",".join(discarded),
                    "samples_run": ",".join(samples_run), "outcome": "checkpoint_incomplete",
                },
            )
            print(
                f"Resume checkpoint retained at {tmp_root}: {len(final_completed)}/{len(selected)} "
                "samples complete; rerun without --sample-id to finish and publish"
            )
            return

        all_worker_counts = sorted(set(final_workers.values()))
        resumed = bool(args.resume_temp_root)
        if resumed:
            append_resume_history(
                tmp_root / "resume_history.tsv",
                {
                    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                    "threads": args.threads,
                    "requested_sample_ids": ",".join(args.sample_id),
                    "cleanup_partial": str(args.cleanup_partial).upper(),
                    "completed_before": len(completed), "missing_before": len(missing),
                    "partial_before": len(partial), "discarded_partial": ",".join(discarded),
                    "samples_run": ",".join(samples_run), "outcome": "published_complete",
                },
            )
        batch_metadata = [
            ("timestamp_utc", datetime.now(timezone.utc).isoformat()),
            ("command", " ".join(shlex.quote(value) for value in sys.argv)),
            ("python_version", platform.python_version()), ("selected_inputs", len(selected)),
            ("execution_order", "manifest order; one Java process at a time"),
            ("configured_nlr_worker_threads_per_process", args.threads),
            ("observed_completed_sample_worker_threads", ",".join(map(str, all_worker_counts))),
            ("mixed_worker_thread_counts", str(len(all_worker_counts) > 1).upper()),
            ("resumed_from_existing_temp_root", str(resumed).upper()),
            ("maximum_allowed_nlr_worker_threads_per_process", MAX_WORKER_THREADS),
            ("jvm_processor_cap_for_new_runs", args.threads), ("jvm_gc", "UseSerialGC"),
            ("nlr_annotator_jar_sha256", tool_hashes["jar"]),
            ("motifs_sha256", tool_hashes["motifs"]), ("store_sha256", tool_hashes["store"]),
            ("completion_status", "complete"),
        ]
        write_tsv(
            tmp_root / "batch_metadata.tsv", ["key", "value"],
            ({"key": key, "value": value} for key, value in batch_metadata),
        )
        os.replace(tmp_root, args.output_root)
    except Exception as exc:
        if args.resume_temp_root:
            raise SystemExit(
                f"ERROR: {exc}; resumable state retained in place at {tmp_root}. "
                "A newly partial sample will fail closed on the next attempt"
            )
        failure_path: Path | None = None
        try:
            (tmp_root / "FAILED.txt").write_text(f"{exc}\n", encoding="utf-8")
            failure_path = Path(str(args.output_root) + ".failed")
            if failure_path.exists():
                failure_path = Path(str(failure_path) + "." + datetime.now().strftime("%Y%m%d%H%M%S"))
            os.replace(tmp_root, failure_path)
        except Exception:
            failure_path = tmp_root
        raise SystemExit(f"ERROR: {exc}; partial audit retained at {failure_path}")

    print(f"Completed {len(selected)} sequential NLR-Annotator runs: {args.output_root}")


if __name__ == "__main__":
    main()
