#!/usr/bin/env python3
"""Validate and atomically merge disjoint NLR-Annotator batch roots.

This utility is for a controlled cutover in which an interrupted runner root
contains a prefix of completed samples and the remaining manifest rows were
completed in one or more disjoint continuation roots.  The original manifest
and current FASTA files remain authoritative.  Source roots are read-only:
the program copies only fully audited sample directories and never renames,
deletes, or edits a source path.

Without ``--execute`` the complete validation is performed but no output is
written.  Publication requires exactly one validated completed directory for
every manifest row and uses a temporary sibling followed by ``os.replace`` on
the output filesystem.
"""

from __future__ import annotations

import argparse
import csv
import os
import platform
import re
import shlex
import shutil
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

from run_nlr_annotator_batch import (
    CHECKSUM_OUTPUTS,
    MAX_WORKER_THREADS,
    MEMORY,
    SELECTED_INPUT_FIELDS,
    InputRow,
    checked_tool_hash,
    input_row_record,
    read_key_value,
    read_manifest,
    read_resume_snapshot,
    read_tsv,
    sha256_file,
    validate_completed_sample,
    write_tsv,
)


SAFE_LABEL = re.compile(r"^[A-Za-z0-9_.-]+$")
SOURCE_OPTIONAL_FILES = {"batch_metadata.tsv", "resume_history.tsv"}
MERGE_PROVENANCE_FIELDS = [
    "manifest_order",
    "sample_id",
    "species",
    "ploidy",
    "analysis_role",
    "input_scope",
    "relative_fasta",
    "input_fasta_sha256",
    "input_fasta_records",
    "input_fasta_total_bases",
    "nlr_output_rows",
    "nlr_output_sequence_ids",
    "nlr_output_locus_ids",
    "source_root_index",
    "source_root_label",
    "source_root",
    "source_selected_inputs_sha256",
    "source_batch_status",
    "configured_nlr_worker_threads",
    "worker_lane_provenance",
    "jvm_processor_cap",
    "nlr_annotator_jar_sha256",
    "motifs_sha256",
    "store_sha256",
    "nlr_calls_sha256",
    "nlr_loci_gff_sha256",
    "stdout_log_sha256",
    "stderr_log_sha256",
]


@dataclass(frozen=True)
class SourceSpec:
    index: int
    label: str
    root: Path


@dataclass(frozen=True)
class Candidate:
    row: InputRow
    sample_dir: Path
    source: SourceSpec
    source_selected_inputs_sha256: str
    source_batch_status: str


@dataclass(frozen=True)
class ValidatedSample:
    candidate: Candidate
    metadata: dict[str, str]
    checksums: dict[str, str]
    workers: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--jar", type=Path, required=True)
    parser.add_argument("--motifs", type=Path, required=True)
    parser.add_argument("--store", type=Path, required=True)
    parser.add_argument(
        "--source",
        action="append",
        required=True,
        metavar="LABEL=ROOT",
        help=(
            "Read-only runner root with an explicit provenance label; repeat as "
            "original14=ROOT, lane1=ROOT, and lane2=ROOT"
        ),
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--expected-samples", type=int, default=22)
    parser.add_argument(
        "--expected-source-worker-threads",
        default="original14=4,lane1=5,lane2=5",
        help="Comma-separated exact per-sample NLR worker counts by source label",
    )
    parser.add_argument(
        "--expected-source-sample-counts",
        default="original14=14,lane1=4,lane2=4",
        help="Comma-separated exact completed-sample counts by source label",
    )
    parser.add_argument(
        "--continuation-aggregate-workers",
        type=int,
        default=10,
        help=(
            "Recorded aggregate scientific concurrency of lane1 plus lane2; this is "
            "never interpreted as a per-sample NLR worker count"
        ),
    )
    parser.add_argument("--expected-jar-sha256")
    parser.add_argument("--expected-motifs-sha256")
    parser.add_argument("--expected-store-sha256")
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Publish the validated merge; otherwise perform a no-write audit",
    )
    return parser.parse_args()


def parse_expected_source_counts(value: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for token in value.split(","):
        token = token.strip()
        if not token or "=" not in token:
            raise ValueError(
                "--expected-source-sample-counts must contain LABEL=COUNT entries"
            )
        label, raw_count = token.split("=", 1)
        if not SAFE_LABEL.fullmatch(label) or label in counts:
            raise ValueError(f"Unsafe or duplicate source-count label: {label!r}")
        try:
            count = int(raw_count)
        except ValueError as exc:
            raise ValueError(f"Invalid sample count for source {label!r}: {raw_count!r}") from exc
        if count < 1:
            raise ValueError(f"Expected source sample count must be positive: {label}={count}")
        counts[label] = count
    if not counts:
        raise ValueError("No expected source sample counts were declared")
    return counts


def parse_expected_source_workers(value: str) -> dict[str, int]:
    workers = parse_expected_source_counts(value)
    for label, count in workers.items():
        if count > MAX_WORKER_THREADS:
            raise ValueError(
                f"Expected per-sample worker count exceeds the runner maximum: {label}={count}"
            )
    return workers


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def build_source_specs(args: argparse.Namespace) -> list[SourceSpec]:
    roots: set[Path] = set()
    labels: set[str] = set()
    specs: list[SourceSpec] = []
    output = args.output_root.resolve(strict=False)
    for index, value in enumerate(args.source, 1):
        if "=" not in value:
            raise ValueError(f"--source must use LABEL=ROOT syntax: {value!r}")
        label, raw_path = value.split("=", 1)
        if not SAFE_LABEL.fullmatch(label):
            raise ValueError(f"Unsafe source label: {label!r}")
        if label in labels:
            raise ValueError(f"Duplicate source label: {label!r}")
        labels.add(label)
        if not raw_path:
            raise ValueError(f"Empty source root for label {label!r}")
        raw_root = Path(raw_path)
        if raw_root.is_symlink() or not raw_root.is_dir():
            raise ValueError(f"Source root must be an existing non-symlink directory: {raw_root}")
        root = raw_root.resolve()
        if root in roots:
            raise ValueError(f"Duplicate source root: {root}")
        roots.add(root)
        if output == root or is_relative_to(output, root):
            raise ValueError(f"Output root must not equal or be nested below a source root: {root}")
        specs.append(SourceSpec(index=index, label=label, root=root))
    return specs


def validate_regular_file(path: Path, description: str) -> None:
    if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"{description} must be a non-empty regular non-symlink file: {path}")


def parse_positive_int(value: str | None, description: str) -> int:
    try:
        parsed = int(value or "")
    except ValueError as exc:
        raise ValueError(f"{description} is not an integer: {value!r}") from exc
    if parsed < 1:
        raise ValueError(f"{description} must be positive: {parsed}")
    return parsed


def read_sample_checksums(sample_dir: Path) -> dict[str, str]:
    fields, rows = read_tsv(sample_dir / "output_checksums.tsv")
    if fields != ["path", "sha256"]:
        raise ValueError(f"Unexpected output-checksum header: {sample_dir}")
    checksums = {row["path"]: row["sha256"].lower() for row in rows}
    if len(checksums) != len(rows) or set(checksums) != set(CHECKSUM_OUTPUTS):
        raise ValueError(f"Invalid four-file checksum inventory: {sample_dir}")
    return checksums


def validate_source_batch_metadata(
    source: SourceSpec,
    snapshot_count: int,
    completed_ids: set[str],
    partial_ids: set[str],
    sample_workers: set[int],
    tool_hashes: dict[str, str],
) -> str:
    path = source.root / "batch_metadata.tsv"
    if not path.exists():
        if source.label == "original14":
            return "interrupted_no_batch_metadata"
        raise ValueError(
            "Continuation source lacks completed batch metadata: "
            f"label={source.label}, root={source.root}"
        )
    validate_regular_file(path, "Source batch metadata")
    metadata = read_key_value(path)
    if metadata.get("completion_status") != "complete":
        raise ValueError(f"Source batch is not marked complete: {source.root}")
    recorded_selected = parse_positive_int(
        metadata.get("selected_inputs"), f"selected_inputs in {path}"
    )
    if recorded_selected != snapshot_count:
        raise ValueError(
            f"Source batch selected-input count mismatch at {source.root}: "
            f"metadata={recorded_selected}, snapshot={snapshot_count}"
        )
    if len(completed_ids) != snapshot_count or partial_ids:
        raise ValueError(
            f"Source marked complete but does not contain exactly its selected inputs: "
            f"{source.root}"
        )
    for key, expected in [
        ("nlr_annotator_jar_sha256", tool_hashes["jar"]),
        ("motifs_sha256", tool_hashes["motifs"]),
        ("store_sha256", tool_hashes["store"]),
    ]:
        if metadata.get(key) != expected:
            raise ValueError(
                f"Source batch tool checksum mismatch for {key} at {source.root}"
            )
    worker_inventory_key = "observed_completed_sample_worker_threads"
    if worker_inventory_key in metadata:
        observed_text = metadata[worker_inventory_key]
        if not observed_text:
            raise ValueError(f"Source batch has an empty worker inventory at {source.root}")
        try:
            recorded_workers = {int(value) for value in observed_text.split(",") if value}
        except ValueError as exc:
            raise ValueError(f"Invalid source worker inventory at {source.root}") from exc
        if recorded_workers != sample_workers:
            raise ValueError(
                f"Source batch worker inventory mismatch at {source.root}: "
                f"metadata={sorted(recorded_workers)}, samples={sorted(sample_workers)}"
            )
    else:
        # The two production continuation controllers predate the explicit
        # observed-worker inventory field.  Accept only that missing-field
        # schema after independently validated per-sample metadata agrees with
        # all older batch fields and the recorded command.
        if source.label not in {"lane1", "lane2"} or len(sample_workers) != 1:
            raise ValueError(f"Source batch lacks worker inventory at {source.root}")
        worker = next(iter(sample_workers))
        for key, expected in [
            ("configured_nlr_worker_threads_per_process", worker),
            ("jvm_processor_cap", worker),
            ("maximum_allowed_nlr_worker_threads_per_process", MAX_WORKER_THREADS),
        ]:
            if parse_positive_int(metadata.get(key), f"{key} in {path}") != expected:
                raise ValueError(
                    f"Legacy source batch worker mismatch for {key} at {source.root}"
                )
        if metadata.get("jvm_gc") != "UseSerialGC":
            raise ValueError(f"Legacy source batch JVM GC mismatch at {source.root}")
        try:
            command = shlex.split(metadata.get("command", ""))
        except ValueError as exc:
            raise ValueError(f"Invalid legacy source command at {source.root}") from exc
        thread_values = [
            command[index + 1]
            for index, token in enumerate(command[:-1])
            if token == "--threads"
        ]
        sample_values = [
            command[index + 1]
            for index, token in enumerate(command[:-1])
            if token == "--sample-id"
        ]
        if thread_values != [str(worker)] or "--execute" not in command:
            raise ValueError(f"Legacy source command worker/execute mismatch at {source.root}")
        if len(sample_values) != len(set(sample_values)) or set(sample_values) != completed_ids:
            raise ValueError(f"Legacy source command sample set mismatch at {source.root}")
    return "complete"


def inventory_sources(
    specs: list[SourceSpec],
    manifest_rows: list[InputRow],
) -> tuple[dict[str, Candidate], dict[int, tuple[int, set[str], set[str]]]]:
    manifest_by_id = {row.sample_id: row for row in manifest_rows}
    candidates: dict[str, Candidate] = {}
    source_state: dict[int, tuple[int, set[str], set[str]]] = {}

    for source in specs:
        snapshot = source.root / "selected_inputs.tsv"
        validate_regular_file(snapshot, "Source selected-input snapshot")
        selected = read_resume_snapshot(source.root, manifest_rows, "all")
        selected_ids = {row.sample_id for row in selected}
        if len(selected_ids) != len(selected):
            raise ValueError(f"Duplicate selected-input sample in {snapshot}")

        allowed_files = {"selected_inputs.tsv", *SOURCE_OPTIONAL_FILES}
        completed_ids: set[str] = set()
        partial_ids: set[str] = set()
        for entry in source.root.iterdir():
            if entry.name in allowed_files:
                if entry.is_symlink() or not entry.is_file():
                    raise ValueError(f"Unexpected source audit-file type: {entry}")
                continue
            if entry.name not in manifest_by_id or entry.name not in selected_ids:
                raise ValueError(f"Unexpected entry in source root: {entry}")
            if entry.is_symlink() or not entry.is_dir():
                raise ValueError(f"Source sample path is not a regular directory: {entry}")
            has_metadata = (entry / "run_metadata.tsv").exists()
            has_checksums = (entry / "output_checksums.tsv").exists()
            if has_metadata != has_checksums:
                raise ValueError(
                    f"Ambiguous completed-looking sample has only one audit file: {entry}"
                )
            if has_metadata:
                completed_ids.add(entry.name)
                if entry.name in candidates:
                    previous = candidates[entry.name].source.root
                    raise ValueError(
                        f"Duplicate completed sample {entry.name!r} in {previous} and {source.root}"
                    )
                candidates[entry.name] = Candidate(
                    row=manifest_by_id[entry.name],
                    sample_dir=entry,
                    source=source,
                    source_selected_inputs_sha256=sha256_file(snapshot),
                    source_batch_status="pending_validation",
                )
            else:
                partial_ids.add(entry.name)
        source_state[source.index] = (len(selected), completed_ids, partial_ids)

    expected_ids = set(manifest_by_id)
    observed_ids = set(candidates)
    missing = sorted(expected_ids - observed_ids)
    extra = sorted(observed_ids - expected_ids)
    if missing or extra:
        raise ValueError(
            "Completed sample coverage does not equal the authoritative manifest: "
            f"missing={missing}, extra={extra}"
        )
    return candidates, source_state


def validate_candidates(
    specs: list[SourceSpec],
    manifest_rows: list[InputRow],
    candidates: dict[str, Candidate],
    source_state: dict[int, tuple[int, set[str], set[str]]],
    tool_hashes: dict[str, str],
    expected_workers_by_source: dict[str, int],
    invocation_cwd: Path,
    input_root_lexical: Path,
) -> list[ValidatedSample]:
    validated: list[ValidatedSample] = []
    workers_by_source: dict[int, set[int]] = {source.index: set() for source in specs}

    for row in manifest_rows:
        candidate = candidates[row.sample_id]
        metadata = read_key_value(candidate.sample_dir / "run_metadata.tsv")
        recorded_fasta_text = metadata.get("input_fasta", "")
        if not recorded_fasta_text:
            raise ValueError(f"Sample {row.sample_id} lacks recorded input_fasta")
        recorded_fasta = Path(recorded_fasta_text)
        if recorded_fasta.is_absolute():
            recorded_for_resolution = recorded_fasta
        else:
            # Production runners were launched from the repository root with a
            # relative --input-root, so their immutable metadata intentionally
            # records paths such as data/linked/.../sample.fa.  Accept only the
            # exact safe lexical path reconstructed from the current manifest;
            # resolution, FASTA hash/record/base audits below still prove file
            # identity and content.  This does not permit path traversal or a
            # different relative alias that merely resolves to the same file.
            if ".." in recorded_fasta.parts:
                raise ValueError(
                    f"Sample {row.sample_id} recorded an unsafe relative input_fasta"
                )
            if input_root_lexical.is_absolute() or ".." in input_root_lexical.parts:
                raise ValueError(
                    "A relative recorded FASTA requires the original --input-root to be "
                    f"safe and relative: {input_root_lexical}"
                )
            expected_relative = input_root_lexical / Path(row.relative_fasta)
            expected_relative_text = expected_relative.as_posix()
            if recorded_fasta_text != expected_relative_text:
                raise ValueError(
                    f"Sample {row.sample_id} recorded a relative input FASTA that "
                    "does not exactly match the current manifest path relative to the "
                    f"invocation base: recorded={recorded_fasta_text!r}, "
                    f"expected={expected_relative_text!r}"
                )
            recorded_for_resolution = invocation_cwd / recorded_fasta
        try:
            recorded_resolved = recorded_for_resolution.resolve(strict=True)
            current_resolved = row.fasta.resolve(strict=True)
        except OSError as exc:
            raise ValueError(f"Could not resolve current FASTA for {row.sample_id}: {exc}") from exc
        if recorded_resolved != current_resolved:
            raise ValueError(
                f"Sample {row.sample_id} recorded a different input FASTA: "
                f"recorded={recorded_resolved}, current={current_resolved}"
            )
        sequences_per_thread = parse_positive_int(
            metadata.get("sequences_per_thread"),
            f"sequences_per_thread for {row.sample_id}",
        )
        max_heap = metadata.get("java_max_heap", "")
        if not MEMORY.fullmatch(max_heap):
            raise ValueError(f"Invalid java_max_heap for {row.sample_id}: {max_heap!r}")
        validation_args = SimpleNamespace(
            sequences_per_thread=sequences_per_thread,
            max_heap=max_heap,
        )
        # Preserve the exact lexical path written by the runner for its
        # metadata comparison after proving that it resolves to the current
        # manifest FASTA.  This handles /var versus /private/var aliases on
        # macOS without weakening file-identity or content checks.
        workers = validate_completed_sample(
            validation_args,
            row,
            candidate.sample_dir,
            tool_hashes,
            expected_input_fasta_text=recorded_fasta_text,
            input_fasta_audit_path=current_resolved,
        )
        expected_workers = expected_workers_by_source[candidate.source.label]
        if workers != expected_workers:
            raise ValueError(
                f"Sample {row.sample_id} used {workers} NLR workers; production merge "
                f"requires {expected_workers} for source {candidate.source.label}"
            )
        checksums = read_sample_checksums(candidate.sample_dir)
        workers_by_source[candidate.source.index].add(workers)
        validated.append(ValidatedSample(candidate, metadata, checksums, workers))

    source_status: dict[int, str] = {}
    for source in specs:
        snapshot_count, completed_ids, partial_ids = source_state[source.index]
        source_status[source.index] = validate_source_batch_metadata(
            source,
            snapshot_count,
            completed_ids,
            partial_ids,
            workers_by_source[source.index],
            tool_hashes,
        )

    return [
        ValidatedSample(
            Candidate(
                row=item.candidate.row,
                sample_dir=item.candidate.sample_dir,
                source=item.candidate.source,
                source_selected_inputs_sha256=item.candidate.source_selected_inputs_sha256,
                source_batch_status=source_status[item.candidate.source.index],
            ),
            item.metadata,
            item.checksums,
            item.workers,
        )
        for item in validated
    ]


def provenance_record(index: int, item: ValidatedSample) -> dict[str, object]:
    row = item.candidate.row
    metadata = item.metadata
    return {
        "manifest_order": index,
        "sample_id": row.sample_id,
        "species": row.species,
        "ploidy": row.ploidy,
        "analysis_role": row.analysis_role,
        "input_scope": row.input_scope,
        "relative_fasta": row.relative_fasta,
        "input_fasta_sha256": metadata["input_fasta_sha256"],
        "input_fasta_records": metadata["input_fasta_records"],
        "input_fasta_total_bases": metadata["input_fasta_total_bases"],
        "nlr_output_rows": metadata["nlr_output_rows"],
        "nlr_output_sequence_ids": metadata["nlr_output_sequence_ids"],
        "nlr_output_locus_ids": metadata["nlr_output_locus_ids"],
        "source_root_index": item.candidate.source.index,
        "source_root_label": item.candidate.source.label,
        "source_root": item.candidate.source.root,
        "source_selected_inputs_sha256": item.candidate.source_selected_inputs_sha256,
        "source_batch_status": item.candidate.source_batch_status,
        "configured_nlr_worker_threads": item.workers,
        "worker_lane_provenance": f"{item.workers}-worker",
        "jvm_processor_cap": metadata["jvm_processor_cap"],
        "nlr_annotator_jar_sha256": metadata["nlr_annotator_jar_sha256"],
        "motifs_sha256": metadata["motifs_sha256"],
        "store_sha256": metadata["store_sha256"],
        "nlr_calls_sha256": item.checksums["nlr_calls.txt"],
        "nlr_loci_gff_sha256": item.checksums["nlr_loci.gff"],
        "stdout_log_sha256": item.checksums["stdout.log"],
        "stderr_log_sha256": item.checksums["stderr.log"],
    }


def checksum_inventory(root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in sorted(root.rglob("*"), key=lambda value: value.as_posix()):
        if path.name == "merge_checksums.tsv" and path.parent == root:
            continue
        if path.is_symlink():
            raise ValueError(f"Symlink encountered in canonical merge staging: {path}")
        if path.is_file():
            rows.append({
                "path": path.relative_to(root).as_posix(),
                "sha256": sha256_file(path),
            })
        elif not path.is_dir():
            raise ValueError(f"Unsupported filesystem entry in merge staging: {path}")
    return rows


def verify_checksum_inventory(root: Path) -> None:
    fields, rows = read_tsv(root / "merge_checksums.tsv")
    if fields != ["path", "sha256"]:
        raise ValueError("Unexpected merge checksum header")
    seen: set[str] = set()
    for row in rows:
        relative = Path(row["path"])
        if (
            not row["path"]
            or relative.is_absolute()
            or ".." in relative.parts
            or row["path"] in seen
        ):
            raise ValueError(f"Invalid merge checksum path: {row['path']!r}")
        seen.add(row["path"])
        path = root / relative
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Merge checksum target is not a regular file: {path}")
        if sha256_file(path) != row["sha256"]:
            raise ValueError(f"Merge staging checksum mismatch: {path}")
    expected = {row["path"] for row in checksum_inventory(root)}
    if seen != expected:
        raise ValueError(
            f"Merge checksum coverage mismatch: missing={sorted(expected-seen)}, "
            f"extra={sorted(seen-expected)}"
        )


def publish(
    args: argparse.Namespace,
    specs: list[SourceSpec],
    manifest_rows: list[InputRow],
    validated: list[ValidatedSample],
    tool_hashes: dict[str, str],
) -> None:
    if args.output_root.exists() or args.output_root.is_symlink():
        raise ValueError(f"Output root already exists: {args.output_root}")
    parent = args.output_root.parent
    parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{args.output_root.name}.merge.tmp.", dir=parent))
    try:
        if staging.stat().st_dev != parent.stat().st_dev:
            raise ValueError("Merge staging is not on the output filesystem")
        for item in validated:
            destination = staging / item.candidate.row.sample_id
            shutil.copytree(item.candidate.sample_dir, destination, symlinks=False)
            for name in [*CHECKSUM_OUTPUTS, "run_metadata.tsv", "output_checksums.tsv"]:
                source_hash = sha256_file(item.candidate.sample_dir / name)
                copied_hash = sha256_file(destination / name)
                if source_hash != copied_hash:
                    raise ValueError(
                        f"Copy verification failed for {item.candidate.row.sample_id}/{name}"
                    )

        write_tsv(
            staging / "selected_inputs.tsv",
            SELECTED_INPUT_FIELDS,
            (input_row_record(row) for row in manifest_rows),
        )
        write_tsv(
            staging / "merge_provenance.tsv",
            MERGE_PROVENANCE_FIELDS,
            (
                provenance_record(index, item)
                for index, item in enumerate(validated, 1)
            ),
        )

        worker_counts: dict[int, int] = {}
        for item in validated:
            worker_counts[item.workers] = worker_counts.get(item.workers, 0) + 1
        source_counts: dict[str, int] = {source.label: 0 for source in specs}
        for item in validated:
            source_counts[item.candidate.source.label] += 1
        metadata_rows = [
            ("timestamp_utc", datetime.now(timezone.utc).isoformat()),
            ("command", shlex.join(sys.argv)),
            ("python_version", platform.python_version()),
            ("merge_mode", "validated_copy_from_disjoint_runner_roots"),
            ("manifest", args.manifest.resolve()),
            ("manifest_sha256", sha256_file(args.manifest)),
            ("input_root", args.input_root),
            ("selected_inputs", len(manifest_rows)),
            ("source_roots", len(specs)),
            (
                "source_sample_counts",
                ";".join(f"{label}:{source_counts[label]}" for label in source_counts),
            ),
            (
                "per_sample_worker_thread_counts",
                ";".join(f"{workers}:{worker_counts[workers]}" for workers in sorted(worker_counts)),
            ),
            (
                "expected_per_sample_workers_by_source",
                args.expected_source_worker_threads,
            ),
            ("continuation_lane_labels", "lane1,lane2"),
            ("continuation_lanes_run_concurrently", "TRUE"),
            ("continuation_aggregate_scientific_workers", args.continuation_aggregate_workers),
            (
                "aggregate_concurrency_semantics",
                "lane1 plus lane2 only; never a per-sample NLR-Annotator worker count",
            ),
            ("nlr_annotator_jar_sha256", tool_hashes["jar"]),
            ("motifs_sha256", tool_hashes["motifs"]),
            ("store_sha256", tool_hashes["store"]),
            ("source_roots_modified", "FALSE"),
            ("publication_method", "same-filesystem temporary sibling then os.replace"),
            ("checksum_inventory", "merge_checksums.tsv; excludes itself"),
            ("completion_status", "complete"),
        ]
        write_tsv(
            staging / "batch_metadata.tsv",
            ["key", "value"],
            ({"key": key, "value": value} for key, value in metadata_rows),
        )
        write_tsv(
            staging / "merge_checksums.tsv",
            ["path", "sha256"],
            checksum_inventory(staging),
        )
        verify_checksum_inventory(staging)

        expected_top = {
            "selected_inputs.tsv", "batch_metadata.tsv", "merge_provenance.tsv",
            "merge_checksums.tsv", *(row.sample_id for row in manifest_rows),
        }
        observed_top = {entry.name for entry in staging.iterdir()}
        if observed_top != expected_top:
            raise ValueError(
                f"Canonical staging content mismatch: missing={sorted(expected_top-observed_top)}, "
                f"extra={sorted(observed_top-expected_top)}"
            )
        if args.output_root.exists() or args.output_root.is_symlink():
            raise ValueError(f"Output root appeared during validation: {args.output_root}")
        os.replace(staging, args.output_root)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise


def main() -> None:
    args = parse_args()
    try:
        invocation_cwd = Path.cwd().resolve(strict=True)
        input_root_lexical = args.input_root
        if args.expected_samples < 1:
            raise ValueError("--expected-samples must be positive")
        if not args.manifest.is_file() or args.manifest.stat().st_size == 0:
            raise ValueError(f"Missing or empty manifest: {args.manifest}")
        if not args.input_root.is_dir():
            raise ValueError(f"Missing input root: {args.input_root}")
        if args.output_root.exists() or args.output_root.is_symlink():
            raise ValueError(f"Output root already exists: {args.output_root}")
        for path, label in [
            (args.jar, "NLR-Annotator JAR"),
            (args.motifs, "motifs file"),
            (args.store, "store file"),
        ]:
            validate_regular_file(path, label)

        args.manifest = args.manifest.resolve()
        args.input_root = args.input_root.resolve()
        args.jar = args.jar.resolve()
        args.motifs = args.motifs.resolve()
        args.store = args.store.resolve()
        args.output_root = args.output_root.resolve(strict=False)
        if not 1 <= args.continuation_aggregate_workers <= 10:
            raise ValueError("--continuation-aggregate-workers must be at most ten")
        expected_source_counts = parse_expected_source_counts(args.expected_source_sample_counts)
        expected_source_workers = parse_expected_source_workers(
            args.expected_source_worker_threads
        )
        source_specs = build_source_specs(args)
        source_labels = {source.label for source in source_specs}
        if source_labels != set(expected_source_counts):
            raise ValueError(
                "Source labels do not equal the declared source-count contract: "
                f"observed={sorted(source_labels)}, expected={sorted(expected_source_counts)}"
            )
        if source_labels != set(expected_source_workers):
            raise ValueError(
                "Source labels do not equal the declared worker-count contract: "
                f"observed={sorted(source_labels)}, expected={sorted(expected_source_workers)}"
            )
        continuation_labels = {"lane1", "lane2"}
        if not continuation_labels <= source_labels:
            raise ValueError("Production merge requires explicit lane1 and lane2 sources")
        expected_aggregate = sum(
            expected_source_workers[label] for label in sorted(continuation_labels)
        )
        if args.continuation_aggregate_workers != expected_aggregate:
            raise ValueError(
                "Continuation aggregate concurrency does not reconcile with two concurrent "
                f"lanes: recorded={args.continuation_aggregate_workers}, "
                f"expected={expected_aggregate}"
            )
        manifest_rows = read_manifest(args.manifest, args.input_root)
        if len(manifest_rows) != args.expected_samples:
            raise ValueError(
                f"Authoritative manifest has {len(manifest_rows)} rows; "
                f"expected {args.expected_samples}"
            )
        for row in manifest_rows:
            validate_regular_file(row.fasta, f"Current FASTA for {row.sample_id}")
        tool_hashes = {
            "jar": checked_tool_hash(args.jar, args.expected_jar_sha256, "NLR-Annotator JAR"),
            "motifs": checked_tool_hash(args.motifs, args.expected_motifs_sha256, "motifs file"),
            "store": checked_tool_hash(args.store, args.expected_store_sha256, "store file"),
        }
        candidates, source_state = inventory_sources(source_specs, manifest_rows)
        observed_source_counts = {
            source.label: len(source_state[source.index][1]) for source in source_specs
        }
        if observed_source_counts != expected_source_counts:
            raise ValueError(
                "Completed sample counts by source do not match the production contract: "
                f"observed={observed_source_counts}, expected={expected_source_counts}"
            )
        validated = validate_candidates(
            source_specs,
            manifest_rows,
            candidates,
            source_state,
            tool_hashes,
            expected_source_workers,
            invocation_cwd,
            input_root_lexical,
        )
        if len(validated) != args.expected_samples:
            raise ValueError("Internal validation error: canonical sample count changed")
        if args.execute:
            publish(args, source_specs, manifest_rows, validated, tool_hashes)
    except (OSError, UnicodeError, csv.Error, ValueError) as exc:
        raise SystemExit(f"ERROR: NLR batch merge failed closed: {exc}")

    worker_summary: dict[int, int] = {}
    for item in validated:
        worker_summary[item.workers] = worker_summary.get(item.workers, 0) + 1
    summary = ", ".join(
        f"{count} sample(s) at {workers} workers"
        for workers, count in sorted(worker_summary.items())
    )
    if args.execute:
        print(
            f"Published {len(validated)} validated NLR samples atomically to "
            f"{args.output_root} ({summary})"
        )
    else:
        print(
            f"MERGE PLAN ONLY: validated {len(validated)} unique complete samples from "
            f"{len(source_specs)} read-only roots ({summary}); add --execute to publish"
        )


if __name__ == "__main__":
    main()
