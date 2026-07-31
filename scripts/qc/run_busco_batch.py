#!/usr/bin/env python3
"""Run BUSCO reproducibly for all selected rows in the assembly-QC manifest.

This standard-library-only wrapper provides bounded parallelism, one stdout
and stderr log per sample, safe output names, restart/skip behavior, and a
deterministically ordered ``batch_status.tsv``.  A single invocation runs one
BUSCO mode (``genome`` or ``proteins``); use separate output directories for
the two modes.  Gzip compression is detected from magic bytes.  Because BUSCO
5.8 does not reliably decode compressed FASTA itself, gzip inputs are
stream-decompressed into atomic, content-addressed plain-FASTA stages.  Source
and staged SHA-256 provenance is validated before reuse.

The manifest schema is shared with ``basic_stats.py``::

    sample  current_or_alternative  accession  genome  gff  protein  source_url

Example
-------
python scripts/assembly_qc/run_busco_batch.py \
  --manifest config/assembly_qc_manifest.tsv \
  --output-dir results/assembly_qc/busco/genome \
  --mode genome --lineage data/busco/embryophyta_odb10 \
  --busco "$HOME/.local/bin/busco" --jobs 1 --cpus-per-job 8

BUSCO is run with ``--offline`` because ``--lineage`` must be a local lineage
directory and with ``--opt-out-run-stats`` to disable BUSCO's separate
anonymous run-statistics submission.  ``--offline`` alone is not a telemetry
control.  The wrapper requires the opt-out flag exactly once when planning a
job and revalidates it immediately before execution; the exact command is
retained in the status table and both per-sample logs for an executed job.  A
sample is skipped only when it has a parseable specific short summary and
passes the available engine-integrity checks.  In particular, a
genome run with Miniprot logs is complete only when the log corresponding to
the summary's lineage contains exactly one normal Miniprot completion marker.
An incomplete existing run is resumed with BUSCO ``--restart``.
For a compressed input, restart and completed-summary reuse additionally
require an exact run-input binding.  An unbound or checksum-mismatched legacy
run fails closed and must use a new output directory or an explicit reviewed
``--force`` rerun.  Plain FASTA continues to follow the original behavior.
``--validate-only`` checks inputs and writes the planned status table without
invoking BUSCO.
"""

from __future__ import annotations

import argparse
import codecs
import csv
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from collect_busco import BuscoParseError, parse_short_summary


REQUIRED_MANIFEST_COLUMNS = (
    "sample",
    "current_or_alternative",
    "accession",
    "genome",
    "gff",
    "protein",
    "source_url",
)
SAFE_SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
STATUS_COLUMNS = (
    "sample",
    "mode",
    "input_path",
    "lineage_path",
    "status",
    "exit_code",
    "summary_path",
    "run_dir",
    "stdout_log",
    "stderr_log",
    "message",
    "command",
)
HIDDEN_THREAD_ENVIRONMENT = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
)
GZIP_MAGIC = b"\x1f\x8b"
HASH_CHUNK_BYTES = 1024 * 1024
STAGING_PROVENANCE_SCHEMA_VERSION = 1
RUN_INPUT_BINDING_SCHEMA_VERSION = 1
BUSCO_RUN_STATS_OPT_OUT_FLAG = "--opt-out-run-stats"
STAGING_PROVENANCE_KEYS = {
    "schema_version",
    "sample",
    "mode",
    "source_path",
    "source_bytes",
    "source_sha256",
    "compression_detection",
    "staged_filename",
    "staged_bytes",
    "staged_sha256",
}


class BatchInputError(RuntimeError):
    """Raised for an invalid manifest, executable, lineage, or input asset."""


@dataclass(frozen=True)
class Job:
    """One validated BUSCO sample job."""

    index: int
    sample: str
    mode: str
    # ``input_path`` remains the manifest-declared source so the public status
    # schema and uncompressed-input behavior do not change.  BUSCO receives
    # ``effective_input_path``; the two paths differ only for gzip input.
    input_path: Path
    effective_input_path: Path
    lineage_path: Path
    run_dir: Path
    stdout_log: Path
    stderr_log: Path
    staged_input: "StagedInput | None"
    input_binding_path: Path
    command: tuple[str, ...]


@dataclass(frozen=True)
class RunIntegrity:
    """Result of inspecting one BUSCO sample directory for safe reuse."""

    summary_path: Path | None
    complete: bool
    message: str


@dataclass(frozen=True)
class StagedInput:
    """Checksum-bound plain FASTA materialized from one gzip source."""

    source_path: Path
    source_bytes: int
    source_sha256: str
    staged_path: Path
    staged_bytes: int
    staged_sha256: str
    provenance_path: Path


def resolve_asset(raw_path: str, manifest_dir: Path) -> Path:
    """Resolve absolute or manifest-relative paths without requiring symlinks."""
    candidate = Path(raw_path).expanduser()
    if not candidate.is_absolute():
        candidate = manifest_dir / candidate
    return candidate.resolve()


def resolve_executable(value: str) -> Path:
    """Resolve an executable path or PATH command and verify execute access."""
    if os.sep in value or (os.altsep and os.altsep in value):
        path = Path(value).expanduser().resolve()
        found = str(path)
    else:
        found = shutil.which(value) or ""
        path = Path(found) if found else Path(value)
    if not found or not path.is_file() or not os.access(path, os.X_OK):
        raise BatchInputError(f"BUSCO executable is missing or not executable: {value}")
    return path.resolve()


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Write deterministic JSON with an fsync-before-replace commit.

    The temporary file is created in the destination directory, so
    ``os.replace`` is atomic on the target filesystem.  A failed write never
    leaves a provenance file that could authorize reuse of a partial FASTA.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, raw_temporary = tempfile.mkstemp(
        prefix=f".{path.name}.tmp.",
        dir=path.parent,
    )
    temporary = Path(raw_temporary)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    except (OSError, TypeError, ValueError) as error:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass
        raise BatchInputError(f"Cannot atomically write JSON provenance {path}: {error}") from error


def stable_sha256(path: Path, label: str) -> tuple[int, str]:
    """Hash one file and reject a file that changed while it was read."""
    try:
        before = path.stat()
        digest = hashlib.sha256()
        observed_bytes = 0
        with path.open("rb") as handle:
            while True:
                chunk = handle.read(HASH_CHUNK_BYTES)
                if not chunk:
                    break
                observed_bytes += len(chunk)
                digest.update(chunk)
        after = path.stat()
    except OSError as error:
        raise BatchInputError(f"Cannot checksum {label} {path}: {error}") from error

    before_signature = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    after_signature = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )
    if before_signature != after_signature or observed_bytes != after.st_size:
        raise BatchInputError(
            f"{label.capitalize()} changed while it was checksummed; refusing an unstable input: {path}"
        )
    return observed_bytes, digest.hexdigest()


def has_gzip_magic(path: Path) -> bool:
    """Detect gzip by magic bytes rather than a potentially misleading suffix."""
    try:
        with path.open("rb") as handle:
            return handle.read(len(GZIP_MAGIC)) == GZIP_MAGIC
    except OSError as error:
        raise BatchInputError(f"Cannot inspect BUSCO input compression for {path}: {error}") from error


def require_sha256(value: object, field: str, path: Path) -> str:
    """Return one lowercase SHA-256 string from a provenance document."""
    text = str(value)
    if not re.fullmatch(r"[0-9a-f]{64}", text):
        raise BatchInputError(f"Invalid {field} in staging provenance {path}")
    return text


def canonical_json(payload: Mapping[str, Any]) -> str:
    """Return a type-sensitive canonical representation for exact comparison."""
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def read_staging_provenance(path: Path) -> dict[str, Any]:
    """Read and structurally validate a staging provenance document."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise BatchInputError(f"Cannot validate existing staging provenance {path}: {error}") from error
    if not isinstance(payload, dict):
        raise BatchInputError(f"Staging provenance is not a JSON object: {path}")
    if set(payload) != STAGING_PROVENANCE_KEYS:
        raise BatchInputError(
            f"Staging provenance has an unexpected field set in {path}"
        )
    if (
        type(payload.get("schema_version")) is not int
        or payload.get("schema_version") != STAGING_PROVENANCE_SCHEMA_VERSION
    ):
        raise BatchInputError(
            f"Unsupported staging provenance schema in {path}: {payload.get('schema_version')!r}"
        )
    sample = payload.get("sample")
    mode = payload.get("mode")
    if not isinstance(sample, str) or not SAFE_SAMPLE_RE.fullmatch(sample):
        raise BatchInputError(f"Invalid sample in staging provenance {path}")
    if mode not in {"genome", "proteins"}:
        raise BatchInputError(f"Invalid mode in staging provenance {path}")
    require_sha256(payload.get("source_sha256"), "source_sha256", path)
    staged_sha256 = require_sha256(payload.get("staged_sha256"), "staged_sha256", path)
    expected_filename = f"{sample}.{mode}.{staged_sha256}.fasta"
    staged_filename = payload.get("staged_filename")
    if (
        staged_filename != expected_filename
        or Path(str(staged_filename)).name != staged_filename
    ):
        raise BatchInputError(
            f"Staging provenance filename/checksum mismatch in {path}; refusing ambiguous reuse"
        )
    for field in ("source_bytes", "staged_bytes"):
        value = payload.get(field)
        if type(value) is not int or value < 1:
            raise BatchInputError(f"Invalid {field} in staging provenance {path}")
    source_path = payload.get("source_path")
    if not isinstance(source_path, str) or not source_path or not Path(source_path).is_absolute():
        raise BatchInputError(f"Invalid source_path in staging provenance {path}")
    if payload.get("compression_detection") != "gzip_magic_bytes":
        raise BatchInputError(f"Invalid compression_detection in staging provenance {path}")
    return payload


def _decompress_gzip_to_temporary(source: Path, temporary: Path) -> tuple[int, str]:
    """Stream one gzip FASTA into a temporary plain file and validate text/FASTA form."""
    digest = hashlib.sha256()
    observed_bytes = 0
    first_non_whitespace: int | None = None
    decoder = codecs.getincrementaldecoder("utf-8")("strict")
    try:
        with gzip.open(source, "rb") as input_handle, temporary.open("wb") as output_handle:
            while True:
                chunk = input_handle.read(HASH_CHUNK_BYTES)
                if not chunk:
                    break
                decoder.decode(chunk, final=False)
                if first_non_whitespace is None:
                    candidate = chunk
                    if observed_bytes == 0 and candidate.startswith(codecs.BOM_UTF8):
                        candidate = candidate[len(codecs.BOM_UTF8) :]
                    stripped = candidate.lstrip(b" \t\r\n")
                    if stripped:
                        first_non_whitespace = stripped[0]
                output_handle.write(chunk)
                digest.update(chunk)
                observed_bytes += len(chunk)
            decoder.decode(b"", final=True)
            output_handle.flush()
            os.fsync(output_handle.fileno())
    except (OSError, EOFError, UnicodeError, gzip.BadGzipFile) as error:
        raise BatchInputError(f"Cannot safely decompress gzip FASTA {source}: {error}") from error

    if observed_bytes == 0:
        raise BatchInputError(f"Gzip BUSCO input expands to an empty file: {source}")
    if first_non_whitespace != ord(">"):
        raise BatchInputError(
            f"Gzip BUSCO input does not expand to a FASTA beginning with '>': {source}"
        )
    return observed_bytes, digest.hexdigest()


def stage_gzip_input(
    source: Path,
    output_dir: Path,
    sample: str,
    mode: str,
) -> StagedInput:
    """Materialize and checksum-bind one gzip FASTA for BUSCO.

    A valid existing stage is reused only after both the compressed source and
    the staged plain FASTA pass SHA-256 validation.  Corrupt or malformed
    provenance fails closed.  A genuinely changed source creates a new
    content-addressed FASTA; older content is retained rather than silently
    deleted while an older BUSCO run may still refer to it.
    """
    staging_dir = output_dir / "staged_inputs" / mode
    staging_dir.mkdir(parents=True, exist_ok=True)
    provenance_path = staging_dir / f"{sample}.provenance.json"
    source_bytes, source_sha256 = stable_sha256(source, "gzip source")

    if provenance_path.exists():
        payload = read_staging_provenance(provenance_path)
        if payload["sample"] != sample or payload["mode"] != mode:
            raise BatchInputError(
                f"Staging provenance identity does not match its sample/mode path: {provenance_path}"
            )
        same_source = (
            payload.get("sample") == sample
            and payload.get("mode") == mode
            and payload.get("source_path") == str(source)
            and payload.get("source_bytes") == source_bytes
            and payload.get("source_sha256") == source_sha256
        )
        if same_source:
            staged_path = staging_dir / str(payload["staged_filename"])
            if not staged_path.is_file():
                raise BatchInputError(
                    f"Staging provenance points to a missing plain FASTA: {staged_path}"
                )
            staged_bytes, staged_sha256 = stable_sha256(staged_path, "staged plain FASTA")
            if (
                staged_bytes != payload["staged_bytes"]
                or staged_sha256 != payload["staged_sha256"]
            ):
                raise BatchInputError(
                    f"Staged plain FASTA failed provenance checksum validation: {staged_path}"
                )
            return StagedInput(
                source_path=source,
                source_bytes=source_bytes,
                source_sha256=source_sha256,
                staged_path=staged_path,
                staged_bytes=staged_bytes,
                staged_sha256=staged_sha256,
                provenance_path=provenance_path,
            )

    descriptor, raw_temporary = tempfile.mkstemp(
        prefix=f".{sample}.{mode}.decompress.tmp.",
        dir=staging_dir,
    )
    os.close(descriptor)
    temporary = Path(raw_temporary)
    try:
        staged_bytes, staged_sha256 = _decompress_gzip_to_temporary(source, temporary)
        # Hash the compressed file again after decompression.  A source that
        # changed between the initial fingerprint and the stream read cannot
        # be safely associated with this staged payload.
        final_source_bytes, final_source_sha256 = stable_sha256(source, "gzip source")
        if (final_source_bytes, final_source_sha256) != (source_bytes, source_sha256):
            raise BatchInputError(
                f"Gzip source changed while it was being staged; refusing mixed provenance: {source}"
            )

        staged_filename = f"{sample}.{mode}.{staged_sha256}.fasta"
        staged_path = staging_dir / staged_filename
        if staged_path.exists():
            existing_bytes, existing_sha256 = stable_sha256(
                staged_path, "existing content-addressed plain FASTA"
            )
            if (existing_bytes, existing_sha256) != (staged_bytes, staged_sha256):
                raise BatchInputError(
                    f"Content-addressed staging target is corrupt or colliding: {staged_path}"
                )
            temporary.unlink()
        else:
            os.chmod(temporary, 0o644)
            os.replace(temporary, staged_path)

        payload: dict[str, Any] = {
            "schema_version": STAGING_PROVENANCE_SCHEMA_VERSION,
            "sample": sample,
            "mode": mode,
            "source_path": str(source),
            "source_bytes": source_bytes,
            "source_sha256": source_sha256,
            "compression_detection": "gzip_magic_bytes",
            "staged_filename": staged_filename,
            "staged_bytes": staged_bytes,
            "staged_sha256": staged_sha256,
        }
        atomic_write_json(provenance_path, payload)
        return StagedInput(
            source_path=source,
            source_bytes=source_bytes,
            source_sha256=source_sha256,
            staged_path=staged_path,
            staged_bytes=staged_bytes,
            staged_sha256=staged_sha256,
            provenance_path=provenance_path,
        )
    finally:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass


def prepare_busco_input(
    source: Path,
    output_dir: Path,
    sample: str,
    mode: str,
) -> StagedInput | None:
    """Return gzip staging metadata, or ``None`` for unchanged plain input."""
    if not has_gzip_magic(source):
        return None
    return stage_gzip_input(source, output_dir, sample, mode)


def validate_staged_input_identity(staged: StagedInput, sample: str, mode: str) -> None:
    """Revalidate source, stage, and provenance immediately before BUSCO."""
    payload = read_staging_provenance(staged.provenance_path)
    expected_payload: dict[str, Any] = {
        "schema_version": STAGING_PROVENANCE_SCHEMA_VERSION,
        "sample": sample,
        "mode": mode,
        "source_path": str(staged.source_path),
        "source_bytes": staged.source_bytes,
        "source_sha256": staged.source_sha256,
        "compression_detection": "gzip_magic_bytes",
        "staged_filename": staged.staged_path.name,
        "staged_bytes": staged.staged_bytes,
        "staged_sha256": staged.staged_sha256,
    }
    if canonical_json(payload) != canonical_json(expected_payload):
        raise BatchInputError(
            f"Staging provenance changed after job planning: {staged.provenance_path}"
        )
    source_bytes, source_sha256 = stable_sha256(staged.source_path, "gzip source")
    if (source_bytes, source_sha256) != (staged.source_bytes, staged.source_sha256):
        raise BatchInputError(
            f"Compressed BUSCO source changed after job planning: {staged.source_path}"
        )
    staged_bytes, staged_sha256 = stable_sha256(staged.staged_path, "staged plain FASTA")
    if (staged_bytes, staged_sha256) != (staged.staged_bytes, staged.staged_sha256):
        raise BatchInputError(
            f"Staged BUSCO FASTA changed after job planning: {staged.staged_path}"
        )


def run_input_binding_payload(job: Job) -> dict[str, Any]:
    """Return the checksum-bound input identity for one compressed-input run."""
    if job.staged_input is None:
        raise BatchInputError("Internal error: plain input has no run-input binding")
    staged = job.staged_input
    return {
        "schema_version": RUN_INPUT_BINDING_SCHEMA_VERSION,
        "sample": job.sample,
        "mode": job.mode,
        "source_path": str(staged.source_path),
        "source_bytes": staged.source_bytes,
        "source_sha256": staged.source_sha256,
        "effective_input_path": str(staged.staged_path),
        "effective_input_bytes": staged.staged_bytes,
        "effective_input_sha256": staged.staged_sha256,
        "staging_provenance_path": str(staged.provenance_path),
    }


def validate_run_input_binding(job: Job) -> tuple[bool, str]:
    """Require an exact existing binding before restarting/reusing gzip runs."""
    if job.staged_input is None:
        return True, "Plain input does not require a staging binding"
    path = job.input_binding_path
    if not path.is_file():
        return False, f"run-input binding is missing: {path}"
    try:
        observed = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        return False, f"run-input binding cannot be validated: {error}"
    if not isinstance(observed, dict):
        return False, "run-input binding is not a JSON object"
    expected = run_input_binding_payload(job)
    if canonical_json(observed) != canonical_json(expected):
        return False, "run-input binding does not match the current compressed/staged checksums"
    return True, "Compressed source and staged FASTA match the run-input binding"


def summary_matches_effective_input(summary_path: Path, effective_input: Path) -> tuple[bool, str]:
    """Bind a reusable BUSCO summary to the exact content-addressed input path."""
    try:
        summary = parse_short_summary(summary_path)
    except BuscoParseError as error:
        return False, f"BUSCO summary is not parseable: {error}"
    if not summary.input_path:
        return False, "BUSCO summary lacks input-path provenance"
    try:
        observed = Path(summary.input_path).expanduser().resolve()
    except (OSError, RuntimeError) as error:
        return False, f"BUSCO summary input path cannot be resolved: {error}"
    expected = effective_input.resolve()
    if observed != expected:
        return False, f"BUSCO summary input {observed} does not match staged input {expected}"
    return True, "BUSCO summary is bound to the exact staged FASTA"


def read_manifest(path: Path, mode: str, selected_samples: set[str]) -> list[tuple[str, Path]]:
    """Read and validate BUSCO sample names and mode-specific input paths."""
    path = Path(path).expanduser().resolve()
    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise BatchInputError(f"Cannot open manifest {path}: {error}") from error

    input_column = "genome" if mode == "genome" else "protein"
    rows: list[tuple[str, Path]] = []
    seen: set[str] = set()
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise BatchInputError(f"Manifest has no header: {path}")
        missing = [column for column in REQUIRED_MANIFEST_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise BatchInputError(
                f"Manifest {path} is missing required columns: {', '.join(missing)}"
            )
        for line_number, row in enumerate(reader, start=2):
            if not row or not any((value or "").strip() for value in row.values()):
                continue
            sample = (row.get("sample") or "").strip()
            if not SAFE_SAMPLE_RE.fullmatch(sample):
                raise BatchInputError(
                    f"Manifest line {line_number} has unsafe sample name {sample!r}; "
                    "use 1-128 ASCII letters, digits, dots, underscores, or hyphens, "
                    "starting with a letter or digit"
                )
            if sample in seen:
                raise BatchInputError(f"Duplicate sample in manifest: {sample}")
            seen.add(sample)
            if selected_samples and sample not in selected_samples:
                continue
            raw_input = (row.get(input_column) or "").strip()
            if not raw_input:
                raise BatchInputError(
                    f"Manifest line {line_number} has no {input_column} path for {sample}"
                )
            input_path = resolve_asset(raw_input, path.parent)
            if not input_path.is_file():
                raise BatchInputError(f"{sample}: {input_column} input is not a file: {input_path}")
            rows.append((sample, input_path))

    if selected_samples:
        unknown = sorted(selected_samples - seen)
        if unknown:
            raise BatchInputError(f"Requested samples are absent from the manifest: {', '.join(unknown)}")
    if not rows:
        raise BatchInputError("No manifest rows were selected")
    return rows


def miniprot_completion_count(path: Path) -> int:
    """Count normal Miniprot final-runtime records without loading a large log."""
    count = 0
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if re.match(r"^\[M::main\]\s+Real time:", line):
                    count += 1
    except OSError:
        return -1
    return count


def summary_integrity(run_dir: Path, mode: str, summary_path: Path) -> RunIntegrity:
    """Validate engine evidence associated with one parseable BUSCO summary.

    BUSCO 5 can write a syntactically complete short summary after Miniprot's
    alignment output is interrupted.  The partial predictions then yield a
    deceptively low but parseable result.  Miniprot writes ``[M::main] Real
    time:`` exactly once after a normal alignment exit, so require that marker
    in every matching Miniprot stderr log.  Recursive discovery covers both
    BUSCO layouts observed in this project: logs directly below the sample run
    and logs below ``run_<dataset>/logs``.

    Protein mode is deliberately exempt.  A genome run with no Miniprot logs
    is also accepted so that BUSCO versions using another gene-prediction
    engine are not rejected.
    """
    summary_path = summary_path.resolve()
    if mode != "genome":
        return RunIntegrity(summary_path, True, "Parseable protein BUSCO summary")

    try:
        summary = parse_short_summary(summary_path)
    except BuscoParseError as error:
        return RunIntegrity(summary_path, False, f"BUSCO summary is not parseable: {error}")

    all_logs = sorted(
        run_dir.rglob("miniprot_align_*_err.log"),
        key=lambda item: str(item),
    )
    if not all_logs:
        return RunIntegrity(
            summary_path,
            True,
            "Parseable genome BUSCO summary; no Miniprot alignment log was present",
        )

    expected_name = f"miniprot_align_{summary.dataset}_err.log"
    relevant_logs = [path for path in all_logs if path.name == expected_name]
    if not relevant_logs:
        found = ", ".join(path.name for path in all_logs)
        return RunIntegrity(
            summary_path,
            False,
            f"Miniprot integrity failed: no {expected_name} log matches summary dataset "
            f"{summary.dataset}; found {found}",
        )

    problems: list[str] = []
    for path in relevant_logs:
        count = miniprot_completion_count(path)
        try:
            display_path = str(path.relative_to(run_dir))
        except ValueError:
            display_path = str(path)
        if count < 0:
            problems.append(f"{display_path} could not be read")
        elif count != 1:
            problems.append(
                f"{display_path} has {count} '[M::main] Real time:' completion markers; expected 1"
            )
    if problems:
        return RunIntegrity(
            summary_path,
            False,
            "Miniprot integrity failed: " + "; ".join(problems),
        )
    return RunIntegrity(
        summary_path,
        True,
        f"Miniprot integrity passed in {len(relevant_logs)} matching log(s)",
    )


def inspect_run_integrity(run_dir: Path, mode: str) -> RunIntegrity:
    """Return the first parseable, integrity-valid summary below ``run_dir``.

    If summaries are parseable but all fail integrity, retain the first path
    and diagnostic message so callers can distinguish a corrupted/partial run
    from a run that has not produced any summary.
    """
    if not run_dir.exists():
        return RunIntegrity(None, False, "BUSCO run directory does not exist")
    first_invalid: RunIntegrity | None = None
    for path in sorted(run_dir.rglob("short_summary.specific.*.txt"), key=lambda item: str(item)):
        try:
            parse_short_summary(path)
        except BuscoParseError:
            continue
        inspected = summary_integrity(run_dir, mode, path)
        if inspected.complete:
            return inspected
        if first_invalid is None:
            first_invalid = inspected
    if first_invalid is not None:
        return first_invalid
    return RunIntegrity(None, False, "No parseable specific BUSCO summary was found")


def command_text(command: tuple[str, ...]) -> str:
    """Return a shell-readable command for provenance only."""
    return shlex.join(command)


def validate_busco_privacy_command(command: tuple[str, ...]) -> None:
    """Fail closed unless BUSCO run-stat collection is explicitly disabled."""
    occurrences = command.count(BUSCO_RUN_STATS_OPT_OUT_FLAG)
    if occurrences != 1:
        raise BatchInputError(
            "BUSCO privacy contract requires exactly one "
            f"{BUSCO_RUN_STATS_OPT_OUT_FLAG}; observed {occurrences}"
        )


def result_row(
    job: Job,
    status: str,
    *,
    exit_code: int | str = "",
    summary_path: Path | None = None,
    message: str = "",
) -> dict[str, str]:
    """Build one deterministic batch-status row."""
    return {
        "sample": job.sample,
        "mode": job.mode,
        "input_path": str(job.input_path),
        "lineage_path": str(job.lineage_path),
        "status": status,
        "exit_code": str(exit_code),
        "summary_path": str(summary_path) if summary_path else "",
        "run_dir": str(job.run_dir),
        "stdout_log": str(job.stdout_log),
        "stderr_log": str(job.stderr_log),
        "message": message,
        "command": command_text(job.command),
    }


def run_job(job: Job) -> dict[str, str]:
    """Execute one BUSCO command and capture its two streams separately."""
    try:
        # Revalidate at the execution boundary as well as during planning so a
        # malformed or programmatically modified Job can never start BUSCO.
        validate_busco_privacy_command(job.command)
    except BatchInputError as error:
        return result_row(job, "failed_runner", message=str(error))
    job.stdout_log.parent.mkdir(parents=True, exist_ok=True)
    if job.staged_input is not None:
        try:
            # Jobs can wait in the bounded executor while another sample is
            # running.  Recheck all three identities at the last possible
            # point so a source/stage edit during that wait cannot be hidden
            # behind the earlier planning check.
            validate_staged_input_identity(job.staged_input, job.sample, job.mode)
            # This binding is committed immediately before BUSCO starts.  A
            # partial run may therefore restart only against the same source
            # and staged checksums.  Completed-run reuse additionally checks
            # the input path recorded inside BUSCO's own short summary.
            atomic_write_json(job.input_binding_path, run_input_binding_payload(job))
        except BatchInputError as error:
            return result_row(
                job,
                "failed_runner",
                message=f"Could not validate/bind compressed input: {error}",
            )
    # BUSCO console scripts commonly use ``#!/usr/bin/env python3``.  Prepend
    # the selected executable's environment directory so a non-interactive
    # SSH launch cannot silently resolve that shebang to the system Python
    # instead of the Python that owns the BUSCO installation.
    environment = os.environ.copy()
    executable_dir = str(Path(job.command[0]).parent)
    existing_path = environment.get("PATH", "")
    environment["PATH"] = (
        f"{executable_dir}{os.pathsep}{existing_path}" if existing_path else executable_dir
    )
    # BUSCO's explicit ``--cpu`` value controls its declared tool workers.
    # Numerical Python libraries can otherwise create additional implicit
    # BLAS/OpenMP pools, so keep those hidden pools serial.  This preserves the
    # project-wide eight-CPU ceiling even when the caller's shell exports a
    # larger thread count.
    for variable in HIDDEN_THREAD_ENVIRONMENT:
        environment[variable] = "1"
    try:
        with job.stdout_log.open("w", encoding="utf-8") as stdout_handle, job.stderr_log.open(
            "w", encoding="utf-8"
        ) as stderr_handle:
            stdout_handle.write(f"# command: {command_text(job.command)}\n")
            stderr_handle.write(f"# command: {command_text(job.command)}\n")
            stdout_handle.flush()
            stderr_handle.flush()
            completed = subprocess.run(
                job.command,
                stdout=stdout_handle,
                stderr=stderr_handle,
                text=True,
                check=False,
                env=environment,
            )
    except OSError as error:
        return result_row(job, "failed_runner", message=f"Could not start BUSCO: {error}")

    integrity = inspect_run_integrity(job.run_dir, job.mode)
    summary = integrity.summary_path
    if completed.returncode != 0:
        return result_row(
            job,
            "failed",
            exit_code=completed.returncode,
            summary_path=summary,
            message="BUSCO returned a non-zero exit code; inspect the per-sample logs",
        )
    if summary is None:
        return result_row(
            job,
            "failed_no_summary",
            exit_code=completed.returncode,
            message="BUSCO exited successfully but no complete specific short summary was found",
        )
    if not integrity.complete:
        return result_row(
            job,
            "failed_integrity",
            exit_code=completed.returncode,
            summary_path=summary,
            message=integrity.message,
        )
    if job.staged_input is not None:
        binding_valid, binding_message = validate_run_input_binding(job)
        summary_input_valid, summary_input_message = summary_matches_effective_input(
            summary,
            job.effective_input_path,
        )
        if not binding_valid or not summary_input_valid:
            return result_row(
                job,
                "failed_integrity",
                exit_code=completed.returncode,
                summary_path=summary,
                message=(
                    "Compressed-input provenance failed after BUSCO: "
                    f"{binding_message}; {summary_input_message}"
                ),
            )
    return result_row(job, "completed", exit_code=completed.returncode, summary_path=summary)


def write_status(path: Path, rows: list[dict[str, str]]) -> None:
    """Atomically replace the status TSV in manifest order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def build_jobs(args: argparse.Namespace) -> tuple[list[Job], list[dict[str, str]], list[Job]]:
    """Validate inputs and partition selected rows into immediate and runnable jobs."""
    if args.jobs < 1:
        raise BatchInputError("--jobs must be at least 1")
    if args.cpus_per_job < 1:
        raise BatchInputError("--cpus-per-job must be at least 1")
    if args.jobs * args.cpus_per_job > 8:
        raise BatchInputError(
            "The project CPU cap is 8: --jobs multiplied by --cpus-per-job must not exceed 8"
        )

    busco = resolve_executable(args.busco)
    lineage = Path(args.lineage).expanduser().resolve()
    if not lineage.is_dir():
        raise BatchInputError(f"Local BUSCO lineage directory does not exist: {lineage}")
    selected_samples = set(args.sample or [])
    rows = read_manifest(args.manifest, args.mode, selected_samples)

    output_dir = Path(args.output_dir).expanduser().resolve()
    runs_dir = output_dir / "runs"
    logs_dir = output_dir / "logs"
    jobs: list[Job] = []
    immediate_rows: list[dict[str, str]] = []
    runnable_jobs: list[Job] = []

    for index, (sample, input_path) in enumerate(rows):
        staged_input = prepare_busco_input(
            input_path,
            output_dir,
            sample,
            args.mode,
        )
        effective_input_path = (
            staged_input.staged_path if staged_input is not None else input_path
        )
        run_dir = runs_dir / sample
        input_binding_path = output_dir / "run_input_bindings" / f"{sample}.{args.mode}.json"
        command = [
            str(busco),
            "--in",
            str(effective_input_path),
            "--lineage_dataset",
            str(lineage),
            "--mode",
            args.mode,
            "--out",
            sample,
            "--out_path",
            str(runs_dir),
            "--cpu",
            str(args.cpus_per_job),
            "--offline",
            BUSCO_RUN_STATS_OPT_OUT_FLAG,
        ]
        existing_integrity = inspect_run_integrity(run_dir, args.mode)
        existing_summary = existing_integrity.summary_path

        job = Job(
            index=index,
            sample=sample,
            mode=args.mode,
            input_path=input_path,
            effective_input_path=effective_input_path,
            lineage_path=lineage,
            run_dir=run_dir,
            stdout_log=logs_dir / f"{sample}.stdout.log",
            stderr_log=logs_dir / f"{sample}.stderr.log",
            staged_input=staged_input,
            input_binding_path=input_binding_path,
            command=tuple(command),
        )

        binding_valid, binding_message = validate_run_input_binding(job)
        summary_input_valid = True
        summary_input_message = "Plain input preserves legacy summary reuse behavior"
        if staged_input is not None and existing_summary is not None:
            summary_input_valid, summary_input_message = summary_matches_effective_input(
                existing_summary,
                effective_input_path,
            )
        gzip_run_reusable = (
            staged_input is None or (binding_valid and summary_input_valid)
        )

        if staged_input is not None and run_dir.exists() and not gzip_run_reusable and not args.force:
            raise BatchInputError(
                f"{sample}: existing BUSCO run cannot be safely reused or restarted for the "
                f"compressed input ({binding_message}; {summary_input_message}). Use a new "
                "--output-dir, or review the old run and rerun explicitly with --force."
            )

        mutable_command = list(job.command)
        if args.force:
            mutable_command.append("--force")
        elif run_dir.exists() and not existing_integrity.complete:
            # For gzip input the binding check above proves that this partial
            # directory belongs to the exact content-addressed staged FASTA.
            mutable_command.append("--restart")
        validate_busco_privacy_command(tuple(mutable_command))
        job = Job(
            index=job.index,
            sample=job.sample,
            mode=job.mode,
            input_path=job.input_path,
            effective_input_path=job.effective_input_path,
            lineage_path=job.lineage_path,
            run_dir=job.run_dir,
            stdout_log=job.stdout_log,
            stderr_log=job.stderr_log,
            staged_input=job.staged_input,
            input_binding_path=job.input_binding_path,
            command=tuple(mutable_command),
        )
        jobs.append(job)
        if args.validate_only:
            message = "Inputs validated; BUSCO was not run"
            if staged_input is not None:
                message += "; gzip FASTA was staged and checksum-validated"
            if run_dir.exists() and not existing_integrity.complete:
                message += f"; existing run would restart: {existing_integrity.message}"
            immediate_rows.append(result_row(job, "validated", message=message))
        elif existing_integrity.complete and gzip_run_reusable and not args.force:
            retained_message = (
                f"Existing complete specific short summary retained; {existing_integrity.message}"
            )
            if staged_input is not None:
                retained_message += f"; {binding_message}; {summary_input_message}"
            immediate_rows.append(
                result_row(
                    job,
                    "skipped_complete",
                    summary_path=existing_summary,
                    message=retained_message,
                )
            )
        else:
            runnable_jobs.append(job)

    return jobs, immediate_rows, runnable_jobs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path, help="Seven-column assembly QC TSV")
    parser.add_argument("--output-dir", required=True, type=Path, help="Batch output directory")
    parser.add_argument("--mode", required=True, choices=("genome", "proteins"))
    parser.add_argument("--lineage", required=True, type=Path, help="Local unpacked BUSCO lineage directory")
    parser.add_argument("--busco", default="busco", help="BUSCO executable path or PATH command")
    parser.add_argument("--jobs", type=int, default=1, help="Concurrent BUSCO processes (default: 1)")
    parser.add_argument("--cpus-per-job", type=int, default=4, help="BUSCO CPUs per process (default: 4)")
    parser.add_argument(
        "--sample",
        action="append",
        help="Run only this exact sample; repeat to select multiple samples",
    )
    parser.add_argument("--force", action="store_true", help="Rerun even when a complete summary exists")
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate all selected inputs and write status without invoking BUSCO",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    status_path = Path(args.output_dir).expanduser().resolve() / "batch_status.tsv"
    try:
        jobs, immediate_rows, runnable_jobs = build_jobs(args)
    except BatchInputError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    results_by_index: dict[int, dict[str, str]] = {}
    immediate_by_sample = {row["sample"]: row for row in immediate_rows}
    for job in jobs:
        if job.sample in immediate_by_sample:
            results_by_index[job.index] = immediate_by_sample[job.sample]

    # Write a useful status file before expensive work starts.  Final writes
    # always return to manifest order, so parallel completion order cannot
    # change the table.
    for job in runnable_jobs:
        results_by_index[job.index] = result_row(job, "pending")
    write_status(status_path, [results_by_index[index] for index in sorted(results_by_index)])

    if runnable_jobs:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            future_to_job = {executor.submit(run_job, job): job for job in runnable_jobs}
            for future in as_completed(future_to_job):
                job = future_to_job[future]
                try:
                    results_by_index[job.index] = future.result()
                except Exception as error:  # Defensive: retain a row for every selected sample.
                    results_by_index[job.index] = result_row(
                        job,
                        "failed_runner",
                        message=f"Unexpected runner error: {type(error).__name__}: {error}",
                    )
                write_status(
                    status_path,
                    [results_by_index[index] for index in sorted(results_by_index)],
                )

    final_rows = [results_by_index[index] for index in sorted(results_by_index)]
    failures = [row for row in final_rows if row["status"].startswith("failed")]
    print(f"Wrote {len(final_rows)} BUSCO job statuses to {status_path}")
    if failures:
        print(f"ERROR: {len(failures)} BUSCO jobs failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
