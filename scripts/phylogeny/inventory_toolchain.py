#!/usr/bin/env python3
"""Audit registered scientific executables without publishing their paths.

The external registry is intentionally path-bearing and must stay outside the
repository. The emitted TSV contains only identifiers, sizes, hashes, a
restricted matched probe token, and pass/fail states. Exact-version banners
are distinguished from program-identity-only smoke tests. Probe commands are
read from the reviewed toolchain manifest and are executed without a shell in
an empty temporary working directory and a minimized environment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
import unicodedata
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TOOLCHAIN = PROJECT_ROOT / "config" / "phylogeny" / "toolchain.tsv"
REGISTRY_FIELDS = ("tool_id", "executable_id", "executable_path")
OUTPUT_FIELDS = (
    "audit_timestamp_utc",
    "toolchain_manifest_sha256",
    "registry_sha256",
    "tool_id",
    "executable_id",
    "manifest_version",
    "expected_primary_executable",
    "probe_evidence_level",
    "allowed_probe_exit_codes_json",
    "executable_size_bytes",
    "executable_sha256",
    "post_probe_executable_sha256",
    "executable_identity_stable",
    "probe_exit_code",
    "probe_output_bytes",
    "probe_output_sha256",
    "normalized_probe_output_sha256",
    "stripped_sgr_sequence_count",
    "matched_probe_token",
    "probe_match",
    "version_match",
    "inventory_status",
)
SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9_.+-]+$")
SAFE_MATCH_TOKEN = re.compile(r"^[A-Za-z0-9_.+-]{1,80}$")
ANSI_SGR = re.compile(r"\x1b\[(?:[0-9]{1,3}(?:;[0-9]{1,3})*)?m")
MAX_PROBE_OUTPUT_BYTES = 1_048_576
PROBE_EVIDENCE_LEVELS = {"exact_version_banner", "program_identity_smoke_only"}
PASS_STATUSES = {"PASS_EXACT_VERSION_MATCH", "PASS_PROGRAM_IDENTITY_SMOKE"}
INHERITED_RUNTIME_ENVIRONMENT = (
    "PATH",
    "LD_LIBRARY_PATH",
    "DYLD_LIBRARY_PATH",
    "SYSTEMROOT",
    "WINDIR",
    "PATHEXT",
)


class InventoryError(ValueError):
    """Raised for invalid manifests, registries, or output destinations."""


class ProbeOutputError(ValueError):
    """Raised when probe bytes cannot form a safe normalized text view."""

    def __init__(self, status: str) -> None:
        super().__init__(status)
        self.status = status


@dataclass(frozen=True)
class ToolSpec:
    tool_id: str
    version: str
    executable: str
    probe_args: tuple[str, ...]
    allowed_probe_exit_codes: tuple[int, ...]
    probe_regex: re.Pattern[str]
    probe_evidence_level: str


@dataclass(frozen=True)
class RegistryEntry:
    tool_id: str
    executable_id: str
    executable_path: Path


@dataclass(frozen=True)
class FileIdentity:
    device: int
    inode: int
    mode: int
    size: int
    modified_ns: int
    sha256: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def snapshot_regular_file(path: Path) -> FileIdentity:
    """Hash one opened file and reject mutation during the hash operation."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise InventoryError("registered executable is not a regular file")
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
        after = os.fstat(handle.fileno())
    before_identity = (
        before.st_dev,
        before.st_ino,
        before.st_mode,
        before.st_size,
        before.st_mtime_ns,
    )
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_mode,
        after.st_size,
        after.st_mtime_ns,
    )
    if before_identity != after_identity:
        raise InventoryError("registered executable changed while it was hashed")
    return FileIdentity(
        device=after.st_dev,
        inode=after.st_ino,
        mode=after.st_mode,
        size=after.st_size,
        modified_ns=after.st_mtime_ns,
        sha256=digest.hexdigest(),
    )


def normalize_probe_output(raw_output: bytes) -> tuple[str, int, str]:
    """Return a safe probe view after stripping only ANSI SGR sequences.

    Raw bytes remain the authoritative audit artifact. Matching is performed
    only on strict UTF-8 text with CR/CRLF normalized to LF and standard SGR
    color/style sequences removed. Any other control or escape character is a
    hard failure rather than an invitation to broaden the version regex.
    """
    try:
        decoded = raw_output.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise ProbeOutputError("FAIL_PROBE_OUTPUT_ENCODING") from error

    stripped_sgr_count = 0

    def remove_sgr(match: re.Match[str]) -> str:
        nonlocal stripped_sgr_count
        stripped_sgr_count += 1
        return ""

    normalized = ANSI_SGR.sub(remove_sgr, decoded)
    normalized = normalized.replace("\r\n", "\n").replace("\r", "\n")
    for character in normalized:
        if character in {"\n", "\t"}:
            continue
        if unicodedata.category(character).startswith("C"):
            raise ProbeOutputError("FAIL_PROBE_OUTPUT_CONTROL_SEQUENCE")

    normalized_sha256 = hashlib.sha256(normalized.encode("utf-8")).hexdigest()
    return normalized, stripped_sgr_count, normalized_sha256


def build_probe_environment(temporary_directory: str) -> dict[str, str]:
    """Build a deterministic runtime environment without credentials/proxies."""
    environment = {
        name: os.environ[name]
        for name in INHERITED_RUNTIME_ENVIRONMENT
        if os.environ.get(name)
    }
    environment.setdefault("PATH", os.defpath)
    environment.update(
        {
            "HOME": temporary_directory,
            "TMPDIR": temporary_directory,
            "LC_ALL": "C",
            "LANG": "C",
            "TERM": "dumb",
            "TZ": "UTC",
            "NO_COLOR": "1",
            "PYTHONNOUSERSITE": "1",
            "R_ENVIRON_USER": os.devnull,
            "R_PROFILE_USER": os.devnull,
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "BLIS_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    return environment


def is_within(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def load_toolchain(path: Path) -> dict[str, ToolSpec]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "tool_id",
                "version",
                "executable",
                "probe_args_json",
                "allowed_probe_exit_codes_json",
                "probe_regex",
                "probe_evidence_level",
            }
            if not required.issubset(reader.fieldnames or ()):
                raise InventoryError("toolchain manifest is missing required columns")
            rows = list(reader)
    except OSError as error:
        raise InventoryError("toolchain manifest cannot be read") from error

    specs: dict[str, ToolSpec] = {}
    for row_number, row in enumerate(rows, start=2):
        tool_id = row["tool_id"]
        if not SAFE_IDENTIFIER.fullmatch(tool_id):
            raise InventoryError(f"toolchain row {row_number} has an unsafe tool_id")
        if tool_id in specs:
            raise InventoryError(f"toolchain row {row_number} duplicates tool_id")
        try:
            parsed_args = json.loads(row["probe_args_json"])
        except json.JSONDecodeError as error:
            raise InventoryError(
                f"toolchain row {row_number} has invalid probe-argument JSON"
            ) from error
        if not isinstance(parsed_args, list) or not all(
            isinstance(argument, str) for argument in parsed_args
        ):
            raise InventoryError(
                f"toolchain row {row_number} probe arguments must be a JSON string list"
            )
        try:
            parsed_exit_codes = json.loads(row["allowed_probe_exit_codes_json"])
        except json.JSONDecodeError as error:
            raise InventoryError(
                f"toolchain row {row_number} has invalid allowed-exit-code JSON"
            ) from error
        canonical_exit_codes = json.dumps(parsed_exit_codes, separators=(",", ":"))
        if (
            not isinstance(parsed_exit_codes, list)
            or not parsed_exit_codes
            or any(type(code) is not int for code in parsed_exit_codes)
            or any(code < 0 or code > 255 for code in parsed_exit_codes)
            or parsed_exit_codes != sorted(set(parsed_exit_codes))
            or row["allowed_probe_exit_codes_json"] != canonical_exit_codes
        ):
            raise InventoryError(
                f"toolchain row {row_number} allowed exit codes must be a "
                "canonical, nonempty, sorted, unique JSON integer list in 0..255"
            )
        try:
            probe_regex = re.compile(row["probe_regex"])
        except re.error as error:
            raise InventoryError(
                f"toolchain row {row_number} has an invalid probe regex"
            ) from error
        probe_evidence_level = row["probe_evidence_level"]
        if probe_evidence_level not in PROBE_EVIDENCE_LEVELS:
            raise InventoryError(
                f"toolchain row {row_number} has an invalid probe evidence level"
            )
        specs[tool_id] = ToolSpec(
            tool_id=tool_id,
            version=row["version"],
            executable=row["executable"],
            probe_args=tuple(parsed_args),
            allowed_probe_exit_codes=tuple(parsed_exit_codes),
            probe_regex=probe_regex,
            probe_evidence_level=probe_evidence_level,
        )
    if not specs:
        raise InventoryError("toolchain manifest has no tool rows")
    return specs


def load_registry(path: Path, specs: dict[str, ToolSpec]) -> list[RegistryEntry]:
    resolved_registry = path.expanduser().resolve()
    if is_within(resolved_registry, PROJECT_ROOT):
        raise InventoryError("the path-bearing registry must remain outside the repository")
    try:
        with resolved_registry.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not set(REGISTRY_FIELDS).issubset(reader.fieldnames or ()):
                raise InventoryError("registry is missing required columns")
            rows = list(reader)
    except OSError as error:
        raise InventoryError("registry cannot be read") from error

    entries: list[RegistryEntry] = []
    observed_keys: set[tuple[str, str]] = set()
    for row_number, row in enumerate(rows, start=2):
        tool_id = row["tool_id"]
        executable_id = row["executable_id"]
        if tool_id not in specs:
            raise InventoryError(f"registry row {row_number} names an unknown tool_id")
        if not SAFE_IDENTIFIER.fullmatch(executable_id):
            raise InventoryError(
                f"registry row {row_number} has an unsafe executable_id"
            )
        key = (tool_id, executable_id)
        if key in observed_keys:
            raise InventoryError(
                f"registry row {row_number} duplicates tool_id and executable_id"
            )
        observed_keys.add(key)
        executable_path = Path(row["executable_path"]).expanduser()
        if not executable_path.is_absolute():
            raise InventoryError(
                f"registry row {row_number} executable_path must be absolute"
            )
        entries.append(
            RegistryEntry(
                tool_id=tool_id,
                executable_id=executable_id,
                executable_path=executable_path,
            )
        )
    if not entries:
        raise InventoryError("registry has no executable rows")
    return entries


def base_result(
    *,
    timestamp: str,
    toolchain_sha256: str,
    registry_sha256: str,
    entry: RegistryEntry,
    spec: ToolSpec,
) -> dict[str, str]:
    result = {field: "NA" for field in OUTPUT_FIELDS}
    result.update(
        {
            "audit_timestamp_utc": timestamp,
            "toolchain_manifest_sha256": toolchain_sha256,
            "registry_sha256": registry_sha256,
            "tool_id": entry.tool_id,
            "executable_id": entry.executable_id,
            "manifest_version": spec.version,
            "expected_primary_executable": spec.executable,
            "probe_evidence_level": spec.probe_evidence_level,
            "allowed_probe_exit_codes_json": json.dumps(
                spec.allowed_probe_exit_codes, separators=(",", ":")
            ),
            "executable_identity_stable": "false",
            "probe_match": "false",
            "version_match": (
                "false"
                if spec.probe_evidence_level == "exact_version_banner"
                else "not_tested"
            ),
        }
    )
    return result


def audit_entry(
    entry: RegistryEntry,
    spec: ToolSpec,
    *,
    timestamp: str,
    toolchain_sha256: str,
    registry_sha256: str,
    timeout_seconds: int,
) -> dict[str, str]:
    result = base_result(
        timestamp=timestamp,
        toolchain_sha256=toolchain_sha256,
        registry_sha256=registry_sha256,
        entry=entry,
        spec=spec,
    )
    path = entry.executable_path
    try:
        resolved_path = path.resolve(strict=True)
    except (OSError, RuntimeError):
        result["inventory_status"] = "FAIL_MISSING"
        return result
    if not resolved_path.is_file():
        result["inventory_status"] = "FAIL_NOT_REGULAR_FILE"
        return result
    if not os.access(resolved_path, os.X_OK):
        result["inventory_status"] = "FAIL_NOT_EXECUTABLE"
        return result

    try:
        pre_identity = snapshot_regular_file(resolved_path)
    except (InventoryError, OSError):
        result["inventory_status"] = "FAIL_HASH_ERROR"
        return result
    result["executable_size_bytes"] = str(pre_identity.size)
    result["executable_sha256"] = pre_identity.sha256

    probe_failure: str | None = None
    probe_text = ""
    completed_returncode: int | None = None
    with tempfile.TemporaryDirectory(prefix="toolchain-probe-") as temporary:
        probe_output = Path(temporary) / "probe-output.bin"
        try:
            with probe_output.open("wb") as handle:
                completed = subprocess.run(
                    [str(resolved_path), *spec.probe_args],
                    cwd=temporary,
                    stdin=subprocess.DEVNULL,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    timeout=timeout_seconds,
                    check=False,
                    env=build_probe_environment(temporary),
                )
            completed_returncode = completed.returncode
            result["probe_exit_code"] = str(completed.returncode)
        except subprocess.TimeoutExpired:
            probe_failure = "FAIL_TIMEOUT"
        except OSError:
            probe_failure = "FAIL_PROBE_ERROR"

        try:
            output_size = probe_output.stat().st_size
            result["probe_output_bytes"] = str(output_size)
            result["probe_output_sha256"] = sha256_file(probe_output)
        except OSError:
            probe_failure = probe_failure or "FAIL_PROBE_OUTPUT_ERROR"
        else:
            if output_size > MAX_PROBE_OUTPUT_BYTES:
                probe_failure = probe_failure or "FAIL_PROBE_OUTPUT_TOO_LARGE"
            else:
                try:
                    raw_probe_output = probe_output.read_bytes()
                except OSError:
                    probe_failure = probe_failure or "FAIL_PROBE_OUTPUT_ERROR"
                else:
                    raw_sha256 = hashlib.sha256(raw_probe_output).hexdigest()
                    if (
                        len(raw_probe_output) != output_size
                        or raw_sha256 != result["probe_output_sha256"]
                    ):
                        probe_failure = probe_failure or "FAIL_PROBE_OUTPUT_CHANGED"
                    else:
                        try:
                            (
                                probe_text,
                                stripped_sgr_count,
                                normalized_sha256,
                            ) = normalize_probe_output(raw_probe_output)
                        except ProbeOutputError as error:
                            probe_failure = probe_failure or error.status
                        else:
                            result["stripped_sgr_sequence_count"] = str(
                                stripped_sgr_count
                            )
                            result["normalized_probe_output_sha256"] = (
                                normalized_sha256
                            )

    try:
        post_resolved_path = path.resolve(strict=True)
        post_identity = snapshot_regular_file(post_resolved_path)
    except (InventoryError, OSError, RuntimeError):
        result["inventory_status"] = "FAIL_EXECUTABLE_CHANGED"
        return result
    result["post_probe_executable_sha256"] = post_identity.sha256
    if (
        post_resolved_path != resolved_path
        or post_identity != pre_identity
        or not os.access(post_resolved_path, os.X_OK)
    ):
        result["inventory_status"] = "FAIL_EXECUTABLE_CHANGED"
        return result
    result["executable_identity_stable"] = "true"

    if probe_failure is not None:
        result["inventory_status"] = probe_failure
        return result
    if completed_returncode is None:
        result["inventory_status"] = "FAIL_PROBE_ERROR"
        return result

    match = spec.probe_regex.search(probe_text)
    if match:
        token = match.group(0)
        if SAFE_MATCH_TOKEN.fullmatch(token):
            result["matched_probe_token"] = token
    if completed_returncode not in spec.allowed_probe_exit_codes:
        result["inventory_status"] = "FAIL_PROBE_EXIT"
    elif match is None:
        if spec.probe_evidence_level == "exact_version_banner":
            result["inventory_status"] = "FAIL_VERSION_MISMATCH"
        else:
            result["inventory_status"] = "FAIL_PROGRAM_IDENTITY_MISMATCH"
    elif spec.probe_evidence_level == "exact_version_banner":
        result["probe_match"] = "true"
        result["version_match"] = "true"
        result["inventory_status"] = "PASS_EXACT_VERSION_MATCH"
    else:
        result["probe_match"] = "true"
        result["version_match"] = "not_tested"
        result["inventory_status"] = "PASS_PROGRAM_IDENTITY_SMOKE"
    return result


def render_tsv(rows: Sequence[dict[str, str]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=OUTPUT_FIELDS,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def write_output(content: str, output: Path | None) -> None:
    if output is None:
        sys.stdout.write(content)
        return
    try:
        with output.open("x", encoding="utf-8", newline="") as handle:
            handle.write(content)
    except FileExistsError as error:
        raise InventoryError("output already exists; refusing to overwrite") from error
    except OSError as error:
        raise InventoryError("output cannot be created") from error


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--toolchain", type=Path, default=DEFAULT_TOOLCHAIN)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=30,
        help="per-probe timeout in seconds (1-120; default: 30)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="return nonzero if any executable does not pass its declared probe",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if not 1 <= args.timeout_seconds <= 120:
        print("ERROR: --timeout-seconds must be between 1 and 120", file=sys.stderr)
        return 2
    if args.output is not None and args.output.exists():
        print("ERROR: output already exists; refusing to overwrite", file=sys.stderr)
        return 2
    try:
        toolchain_path = args.toolchain.expanduser().resolve(strict=True)
        registry_path = args.registry.expanduser().resolve(strict=True)
        specs = load_toolchain(toolchain_path)
        entries = load_registry(registry_path, specs)
        timestamp = (
            datetime.now(timezone.utc)
            .replace(microsecond=0)
            .isoformat()
            .replace("+00:00", "Z")
        )
        toolchain_sha256 = sha256_file(toolchain_path)
        registry_sha256 = sha256_file(registry_path)
        results = [
            audit_entry(
                entry,
                specs[entry.tool_id],
                timestamp=timestamp,
                toolchain_sha256=toolchain_sha256,
                registry_sha256=registry_sha256,
                timeout_seconds=args.timeout_seconds,
            )
            for entry in entries
        ]
        write_output(render_tsv(results), args.output)
    except (InventoryError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    if args.strict and any(
        row["inventory_status"] not in PASS_STATUSES for row in results
    ):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
