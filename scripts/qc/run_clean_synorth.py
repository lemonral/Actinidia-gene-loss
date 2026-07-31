#!/usr/bin/env python3
"""Run one SynOrths comparison from an isolated, provenance-recorded copy.

The installed SynOrths package is treated as read-only.  This wrapper copies
only ``SynOrths`` and ``bin/`` into a new per-run directory, changes the single
literal ``-num_threads 30`` in the copied ``bin/initialBlast.pl``, and executes
the copied program from that isolated directory.  Absolute input and output
paths prevent the run directory from changing which biological assets are
used, while the isolated current working directory contains SynOrths scratch
files and avoids reusing cached BLAST intermediates from another comparison.

Example
-------
python scripts/assembly_qc/run_clean_synorth.py \
  --synorth-dir /opt/SynOrths_V1.5 \
  --query-protein target.primary.faa \
  --query-coords target.coords \
  --reference-protein csc.primary.faa \
  --reference-coords csc.coords \
  --output-dir results/synorth/clean/target \
  --output-name target_vs_csc.synorths.txt \
  --m 20 --n 100 --r 0.2 --blast-threads 2

SynOrths 1.5 starts three BLASTP searches concurrently in
``bin/initialBlast.pl``.  A value of 2 therefore limits the BLAST phase to at
most six requested worker threads for one isolated run.

The output directory receives the requested SynOrths result, raw stdout and
stderr logs, a provenance JSON document, a status JSON document, and a hidden
``.clean_synorth_runs`` directory containing the isolated tool/work copy.
``--validate-only`` still makes and patches the isolated copy and writes both
JSON documents, but it does not invoke SynOrths or create the requested result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PATCH_TOKEN = "-num_threads 30"
SAFE_OUTPUT_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,254}$")


class CleanSynorthError(RuntimeError):
    """Raised when inputs, isolation, patching, or execution are unsafe."""


@dataclass(frozen=True)
class RunPaths:
    """Resolved paths for one isolated SynOrths invocation."""

    source_dir: Path
    source_executable: Path
    source_bin: Path
    output_dir: Path
    output_path: Path
    stdout_log: Path
    stderr_log: Path
    provenance_json: Path
    status_json: Path
    isolated_dir: Path
    copied_executable: Path
    copied_bin: Path
    copied_initial_blast: Path


def utc_now() -> str:
    """Return a timezone-explicit ISO 8601 timestamp."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    """Calculate SHA-256 without loading a potentially large file at once."""
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        raise CleanSynorthError(f"Cannot hash {path}: {error}") from error
    return digest.hexdigest()


def file_record(path: Path, role: str, relative_to: Path | None = None) -> dict[str, Any]:
    """Build one checksum record for provenance."""
    try:
        size = path.stat().st_size
    except OSError as error:
        raise CleanSynorthError(f"Cannot stat {path}: {error}") from error
    record: dict[str, Any] = {
        "role": role,
        "path": str(path),
        "size_bytes": size,
        "sha256": sha256_file(path),
    }
    if relative_to is not None:
        record["relative_path"] = str(path.relative_to(relative_to))
    return record


def tool_records(root: Path, executable: Path, bin_dir: Path, role: str) -> list[dict[str, Any]]:
    """Hash the executable and every regular file under the package bin tree."""
    records = [file_record(executable, role, root)]
    try:
        bin_files = sorted(
            (path for path in bin_dir.rglob("*") if path.is_file()),
            key=lambda path: str(path.relative_to(root)),
        )
    except OSError as error:
        raise CleanSynorthError(f"Cannot inspect SynOrths bin directory {bin_dir}: {error}") from error
    records.extend(file_record(path, role, root) for path in bin_files)
    return records


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Atomically write a deterministic, human-readable JSON document."""
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except OSError as error:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass
        raise CleanSynorthError(f"Cannot write JSON file {path}: {error}") from error


def require_nonempty_file(raw_path: Path, label: str, *, executable: bool = False) -> Path:
    """Resolve and validate one required regular, non-empty file."""
    path = raw_path.expanduser().resolve()
    if not path.is_file():
        raise CleanSynorthError(f"{label} is not a regular file: {path}")
    try:
        size = path.stat().st_size
    except OSError as error:
        raise CleanSynorthError(f"Cannot stat {label} {path}: {error}") from error
    if size == 0:
        raise CleanSynorthError(f"{label} is empty: {path}")
    if executable and not os.access(path, os.X_OK):
        raise CleanSynorthError(f"{label} is not executable: {path}")
    return path


def is_within(path: Path, directory: Path) -> bool:
    """Return whether a resolved path is the directory or one of its children."""
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def validate_positive(value: int, option: str) -> None:
    if value < 1:
        raise CleanSynorthError(f"{option} must be at least 1")


def resolve_paths(args: argparse.Namespace) -> tuple[RunPaths, dict[str, Path]]:
    """Validate source/input/output paths and allocate a unique isolated run."""
    validate_positive(args.m, "--m")
    validate_positive(args.n, "--n")
    validate_positive(args.blast_threads, "--blast-threads")
    if args.blast_threads > 2:
        raise CleanSynorthError(
            "--blast-threads must not exceed 2 in this project; SynOrths starts "
            "three concurrent BLASTP searches, so this caps the BLAST phase at six workers"
        )
    if not math.isfinite(args.r) or args.r <= 0:
        raise CleanSynorthError("--r must be a finite number greater than 0")
    if not SAFE_OUTPUT_RE.fullmatch(args.output_name):
        raise CleanSynorthError(
            "--output-name must be a plain 1-255 character filename made from "
            "ASCII letters, digits, dots, underscores, or hyphens and must start "
            "with a letter or digit"
        )

    source_dir = args.synorth_dir.expanduser().resolve()
    if not source_dir.is_dir():
        raise CleanSynorthError(f"SynOrths installation directory does not exist: {source_dir}")
    source_executable = require_nonempty_file(
        source_dir / "SynOrths", "SynOrths executable", executable=True
    )
    source_bin = (source_dir / "bin").resolve()
    if not source_bin.is_dir():
        raise CleanSynorthError(f"SynOrths bin directory does not exist: {source_bin}")
    if not is_within(source_executable, source_dir) or not is_within(source_bin, source_dir):
        raise CleanSynorthError(
            "SynOrths executable and bin directory must resolve inside --synorth-dir; "
            "external package symlinks are not accepted for an isolated auditable copy"
        )
    require_nonempty_file(source_bin / "initialBlast.pl", "SynOrths initialBlast.pl")

    inputs = {
        "query_protein": require_nonempty_file(args.query_protein, "Query protein FASTA"),
        "query_coords": require_nonempty_file(args.query_coords, "Query coordinate table"),
        "reference_protein": require_nonempty_file(
            args.reference_protein, "Reference protein FASTA"
        ),
        "reference_coords": require_nonempty_file(
            args.reference_coords, "Reference coordinate table"
        ),
    }

    output_dir = args.output_dir.expanduser().resolve()
    if is_within(output_dir, source_dir):
        raise CleanSynorthError(
            f"Output directory must not be the SynOrths installation or a child of it: {output_dir}"
        )
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise CleanSynorthError(f"Cannot create output directory {output_dir}: {error}") from error
    if not output_dir.is_dir():
        raise CleanSynorthError(f"Output directory is not a directory: {output_dir}")

    output_path = output_dir / args.output_name
    stdout_log = output_dir / f"{args.output_name}.stdout.log"
    stderr_log = output_dir / f"{args.output_name}.stderr.log"
    provenance_json = output_dir / f"{args.output_name}.provenance.json"
    status_json = output_dir / f"{args.output_name}.status.json"
    protected_files = set(inputs.values()) | {
        source_executable,
        (source_bin / "initialBlast.pl").resolve(),
    }
    generated_files = (output_path, stdout_log, stderr_log, provenance_json, status_json)
    collisions = [path for path in generated_files if path.resolve() in protected_files]
    if collisions:
        raise CleanSynorthError(
            "A generated result, log, or metadata path would overwrite a source or input file: "
            + ", ".join(str(path) for path in collisions)
        )
    if output_path.exists() and output_path.is_dir():
        raise CleanSynorthError(f"Requested output path is a directory: {output_path}")
    if output_path.exists() and output_path.stat().st_size > 0 and not args.force:
        raise CleanSynorthError(
            f"Refusing to replace existing non-empty output without --force: {output_path}"
        )

    run_id = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{os.getpid()}-{uuid.uuid4().hex[:8]}"
    isolated_dir = output_dir / ".clean_synorth_runs" / run_id
    copied_executable = isolated_dir / "SynOrths"
    copied_bin = isolated_dir / "bin"
    return (
        RunPaths(
            source_dir=source_dir,
            source_executable=source_executable,
            source_bin=source_bin,
            output_dir=output_dir,
            output_path=output_path,
            stdout_log=stdout_log,
            stderr_log=stderr_log,
            provenance_json=provenance_json,
            status_json=status_json,
            isolated_dir=isolated_dir,
            copied_executable=copied_executable,
            copied_bin=copied_bin,
            copied_initial_blast=copied_bin / "initialBlast.pl",
        ),
        inputs,
    )


def copy_and_patch(paths: RunPaths, blast_threads: int) -> dict[str, Any]:
    """Make the minimal isolated package and patch its one thread literal."""
    try:
        paths.isolated_dir.mkdir(parents=True, exist_ok=False)
        shutil.copy2(paths.source_executable, paths.copied_executable)
        shutil.copytree(paths.source_bin, paths.copied_bin, symlinks=False)
    except OSError as error:
        raise CleanSynorthError(f"Cannot create isolated SynOrths copy: {error}") from error

    try:
        original_text = paths.copied_initial_blast.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise CleanSynorthError(
            f"Cannot read copied initialBlast.pl as UTF-8: {paths.copied_initial_blast}: {error}"
        ) from error
    occurrences = original_text.count(PATCH_TOKEN)
    if occurrences != 1:
        raise CleanSynorthError(
            f"Expected exactly one literal {PATCH_TOKEN!r} in copied initialBlast.pl; "
            f"found {occurrences}. The source package was not modified."
        )
    replacement = f"-num_threads {blast_threads}"
    patched_text = original_text.replace(PATCH_TOKEN, replacement, 1)
    try:
        paths.copied_initial_blast.write_text(patched_text, encoding="utf-8")
    except OSError as error:
        raise CleanSynorthError(
            f"Cannot patch copied initialBlast.pl {paths.copied_initial_blast}: {error}"
        ) from error
    return {
        "file": str(paths.copied_initial_blast),
        "search_literal": PATCH_TOKEN,
        "replacement_literal": replacement,
        "replacement_count": 1,
    }


def build_command(
    paths: RunPaths, inputs: dict[str, Path], args: argparse.Namespace
) -> tuple[str, ...]:
    """Create the exact no-shell SynOrths command."""
    return (
        str(paths.copied_executable),
        "-a",
        str(inputs["query_protein"]),
        "-b",
        str(inputs["reference_protein"]),
        "-p",
        str(inputs["query_coords"]),
        "-q",
        str(inputs["reference_coords"]),
        "-m",
        str(args.m),
        "-n",
        str(args.n),
        "-r",
        format(args.r, ".15g"),
        "-o",
        str(paths.output_path),
    )


def status_payload(
    *,
    status: str,
    started_at: str,
    finished_at: str | None,
    paths: RunPaths,
    message: str,
    exit_code: int | None,
    output_record: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the compact machine-readable run status document."""
    return {
        "schema_version": 1,
        "status": status,
        "message": message,
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "exit_code": exit_code,
        "output": output_record,
        "output_path": str(paths.output_path),
        "stdout_log": str(paths.stdout_log),
        "stderr_log": str(paths.stderr_log),
        "provenance_json": str(paths.provenance_json),
        "isolated_work_dir": str(paths.isolated_dir),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--synorth-dir",
        required=True,
        type=Path,
        help="Read-only installed package directory containing SynOrths and bin/",
    )
    parser.add_argument(
        "--blast-bin",
        type=Path,
        help=(
            "directory containing executable blastp and makeblastdb; when supplied, "
            "the directory and executable hashes are recorded and PATH is set explicitly"
        ),
    )
    parser.add_argument("--query-protein", required=True, type=Path)
    parser.add_argument("--query-coords", required=True, type=Path)
    parser.add_argument("--reference-protein", required=True, type=Path)
    parser.add_argument("--reference-coords", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--output-name",
        required=True,
        help="Plain output filename passed to SynOrths -o (no directory components)",
    )
    parser.add_argument("--m", type=int, default=20, help="SynOrths -m value (default: 20)")
    parser.add_argument("--n", type=int, default=100, help="SynOrths -n value (default: 100)")
    parser.add_argument("--r", type=float, default=0.2, help="SynOrths -r value (default: 0.2)")
    parser.add_argument(
        "--blast-threads",
        type=int,
        default=2,
        help=(
            "Threads requested by each of the three concurrent BLASTP searches "
            "in initialBlast.pl (default: 2; at most six requested BLAST workers)"
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Permit replacement of an existing non-empty requested output",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate, copy, patch, hash, and write metadata without invoking SynOrths",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    """Prepare and optionally execute a clean SynOrths run."""
    started_at = utc_now()
    paths, inputs = resolve_paths(args)
    if args.blast_bin is None:
        resolved_blastp = shutil.which("blastp")
        resolved_makeblastdb = shutil.which("makeblastdb")
        if resolved_blastp is None or resolved_makeblastdb is None:
            raise CleanSynorthError(
                "blastp and makeblastdb are not both available on PATH; supply --blast-bin"
            )
        blast_bin = Path(resolved_blastp).resolve().parent
        blast_tools = {
            "blastp": require_nonempty_file(Path(resolved_blastp), "blastp", executable=True),
            "makeblastdb": require_nonempty_file(
                Path(resolved_makeblastdb), "makeblastdb", executable=True
            ),
        }
    else:
        blast_bin = args.blast_bin.expanduser().resolve()
        if not blast_bin.is_dir():
            raise CleanSynorthError(f"--blast-bin is not a directory: {blast_bin}")
        blast_tools = {
            "blastp": require_nonempty_file(blast_bin / "blastp", "blastp", executable=True),
            "makeblastdb": require_nonempty_file(
                blast_bin / "makeblastdb", "makeblastdb", executable=True
            ),
        }
    blast_tool_before = [file_record(path, role) for role, path in blast_tools.items()]
    # Existing successful metadata should not be replaced by a refused rerun;
    # resolve_paths performs that collision check before this point.
    write_json(
        paths.status_json,
        status_payload(
            status="preparing",
            started_at=started_at,
            finished_at=None,
            paths=paths,
            message="Validating inputs and preparing an isolated SynOrths copy",
            exit_code=None,
        ),
    )

    source_tool_before = tool_records(
        paths.source_dir, paths.source_executable, paths.source_bin, "source_tool"
    )
    patch_record = copy_and_patch(paths, args.blast_threads)
    copied_tool = tool_records(
        paths.isolated_dir, paths.copied_executable, paths.copied_bin, "isolated_tool"
    )
    input_records = [file_record(path, role) for role, path in inputs.items()]
    command = build_command(paths, inputs, args)

    provenance: dict[str, Any] = {
        "schema_version": 1,
        "status": "validated" if args.validate_only else "prepared",
        "created_at_utc": started_at,
        "finished_at_utc": None,
        "runner": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256_file(Path(__file__).resolve()),
            "python_executable": sys.executable,
            "python_version": sys.version,
        },
        "source_synorth_directory": str(paths.source_dir),
        "blast_bin": str(blast_bin),
        "external_blast_tools_before": blast_tool_before,
        "isolated_work_directory": str(paths.isolated_dir),
        "source_tool_files_before": source_tool_before,
        "isolated_tool_files_after_patch": copied_tool,
        "inputs": input_records,
        "patch": patch_record,
        "parameters": {
            "m": args.m,
            "n": args.n,
            "r": args.r,
            "blast_threads": args.blast_threads,
        },
        "working_directory": str(paths.isolated_dir),
        "command": list(command),
        "command_shell_display_only": shlex.join(command),
        "stdout_log": str(paths.stdout_log),
        "stderr_log": str(paths.stderr_log),
        "requested_output": str(paths.output_path),
        "output": None,
    }

    # Hash the source a second time after patching and fail if any source file
    # changed.  This is a guardrail and a provenance assertion, not merely a
    # consequence of having copied to a different path.
    source_tool_after = tool_records(
        paths.source_dir, paths.source_executable, paths.source_bin, "source_tool"
    )
    provenance["source_tool_files_after"] = source_tool_after
    provenance["source_tool_unchanged"] = source_tool_before == source_tool_after
    if source_tool_before != source_tool_after:
        provenance["status"] = "failed_source_changed"
        provenance["finished_at_utc"] = utc_now()
        write_json(paths.provenance_json, provenance)
        raise CleanSynorthError(
            "The installed SynOrths executable or bin tree changed during preparation; aborting"
        )

    if args.validate_only:
        finished_at = utc_now()
        provenance["finished_at_utc"] = finished_at
        write_json(paths.provenance_json, provenance)
        write_json(
            paths.status_json,
            status_payload(
                status="validated",
                started_at=started_at,
                finished_at=finished_at,
                paths=paths,
                message="Inputs and the isolated thread patch were validated; SynOrths was not run",
                exit_code=None,
            ),
        )
        print(f"Validated clean SynOrths run; provenance: {paths.provenance_json}")
        return 0

    # Remove an empty placeholder, or a non-empty file explicitly authorized
    # by --force, before launch so a failed command cannot be mistaken for a
    # newly generated result.
    if paths.output_path.exists():
        try:
            paths.output_path.unlink()
        except OSError as error:
            raise CleanSynorthError(f"Cannot remove prior output {paths.output_path}: {error}") from error

    write_json(paths.provenance_json, provenance)
    write_json(
        paths.status_json,
        status_payload(
            status="running",
            started_at=started_at,
            finished_at=None,
            paths=paths,
            message="SynOrths is running in the isolated working directory",
            exit_code=None,
        ),
    )

    environment = os.environ.copy()
    environment["PATH"] = os.pathsep.join(
        [
            str(paths.isolated_dir),
            str(paths.copied_bin),
            str(blast_bin),
            environment.get("PATH", ""),
        ]
    )
    try:
        with paths.stdout_log.open("w", encoding="utf-8") as stdout_handle, paths.stderr_log.open(
            "w", encoding="utf-8"
        ) as stderr_handle:
            completed = subprocess.run(
                command,
                cwd=paths.isolated_dir,
                env=environment,
                stdout=stdout_handle,
                stderr=stderr_handle,
                text=True,
                check=False,
            )
    except OSError as error:
        finished_at = utc_now()
        provenance["status"] = "failed_runner"
        provenance["finished_at_utc"] = finished_at
        provenance["runner_error"] = str(error)
        write_json(paths.provenance_json, provenance)
        write_json(
            paths.status_json,
            status_payload(
                status="failed_runner",
                started_at=started_at,
                finished_at=finished_at,
                paths=paths,
                message=f"Could not start copied SynOrths executable: {error}",
                exit_code=None,
            ),
        )
        print(f"ERROR: could not start copied SynOrths executable: {error}", file=sys.stderr)
        return 1

    finished_at = utc_now()
    provenance["exit_code"] = completed.returncode
    provenance["finished_at_utc"] = finished_at
    if completed.returncode != 0:
        provenance["status"] = "failed"
        write_json(paths.provenance_json, provenance)
        write_json(
            paths.status_json,
            status_payload(
                status="failed",
                started_at=started_at,
                finished_at=finished_at,
                paths=paths,
                message="SynOrths returned a non-zero exit code; inspect stdout and stderr logs",
                exit_code=completed.returncode,
            ),
        )
        print(
            f"ERROR: SynOrths exited with code {completed.returncode}; inspect {paths.stderr_log}",
            file=sys.stderr,
        )
        return 1

    if not paths.output_path.is_file() or paths.output_path.stat().st_size == 0:
        provenance["status"] = "failed_no_output"
        write_json(paths.provenance_json, provenance)
        write_json(
            paths.status_json,
            status_payload(
                status="failed_no_output",
                started_at=started_at,
                finished_at=finished_at,
                paths=paths,
                message="SynOrths exited successfully but did not create a non-empty requested output",
                exit_code=completed.returncode,
            ),
        )
        print(
            f"ERROR: SynOrths exited successfully but output is missing or empty: {paths.output_path}",
            file=sys.stderr,
        )
        return 1

    output_record = file_record(paths.output_path, "synorth_output")
    source_tool_final = tool_records(
        paths.source_dir, paths.source_executable, paths.source_bin, "source_tool"
    )
    blast_tool_final = [file_record(path, role) for role, path in blast_tools.items()]
    provenance["external_blast_tools_final"] = blast_tool_final
    provenance["external_blast_tools_unchanged"] = blast_tool_final == blast_tool_before
    if blast_tool_final != blast_tool_before:
        provenance["status"] = "failed_external_tool_changed"
        write_json(paths.provenance_json, provenance)
        write_json(
            paths.status_json,
            status_payload(
                status="failed_external_tool_changed",
                started_at=started_at,
                finished_at=finished_at,
                paths=paths,
                message="blastp or makeblastdb changed during execution",
                exit_code=completed.returncode,
                output_record=output_record,
            ),
        )
        print("ERROR: blastp or makeblastdb changed during execution", file=sys.stderr)
        return 1
    if source_tool_final != source_tool_before:
        provenance["status"] = "failed_source_changed"
        provenance["source_tool_files_final"] = source_tool_final
        provenance["source_tool_unchanged"] = False
        write_json(paths.provenance_json, provenance)
        write_json(
            paths.status_json,
            status_payload(
                status="failed_source_changed",
                started_at=started_at,
                finished_at=finished_at,
                paths=paths,
                message="Installed SynOrths files changed during execution; result requires investigation",
                exit_code=completed.returncode,
                output_record=output_record,
            ),
        )
        print("ERROR: installed SynOrths files changed during execution", file=sys.stderr)
        return 1

    provenance["status"] = "completed"
    provenance["source_tool_files_final"] = source_tool_final
    provenance["source_tool_unchanged"] = True
    provenance["output"] = output_record
    write_json(paths.provenance_json, provenance)
    write_json(
        paths.status_json,
        status_payload(
            status="completed",
            started_at=started_at,
            finished_at=finished_at,
            paths=paths,
            message="SynOrths completed and created a non-empty output",
            exit_code=completed.returncode,
            output_record=output_record,
        ),
    )
    print(f"Completed clean SynOrths run: {paths.output_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except CleanSynorthError as error:
        # A preparing status exists only after path/collision validation.  If
        # possible, replace it with a failed-validation status without hiding
        # a previously completed run rejected by the collision check.
        output_dir = args.output_dir.expanduser().resolve()
        status_path = output_dir / f"{args.output_name}.status.json"
        if status_path.exists():
            try:
                existing = json.loads(status_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                existing = {}
            if existing.get("status") in {"preparing", "running"}:
                existing.update(
                    {
                        "status": "failed_validation",
                        "message": str(error),
                        "finished_at_utc": utc_now(),
                    }
                )
                try:
                    write_json(status_path, existing)
                except CleanSynorthError:
                    pass
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
