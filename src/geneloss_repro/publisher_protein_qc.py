"""Fail-closed compatibility QC for derived and publisher protein FASTA files.

The comparison is deliberately narrower than a general sequence-alignment
tool.  Record identifiers must reconcile exactly.  After removal of at most
one terminal ``*`` from each sequence, residues must be identical except that
an ``X`` in the publisher sequence may match one unambiguous, canonical amino
acid in the derived sequence.  Every use of that one-way wildcard is recorded.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from .io_utils import natural_key


WORKFLOW_VERSION = "1.0.0"
GZIP_MAGIC = b"\x1f\x8b"
SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
# B, J, O, U, and Z may occur in publisher files.  They can pass only by an
# exact character match; none is accepted as a wildcard.
PROTEIN_ALPHABET = CANONICAL_AMINO_ACIDS | frozenset("XBJOUZ*")


class PublishedProteinCompatibilityError(RuntimeError):
    """Raised when a publisher-protein compatibility bundle cannot pass."""


@dataclass(frozen=True)
class ProteinCompatibilityResult:
    """Counts from one atomically published compatibility audit."""

    output_dir: Path
    record_count: int
    exact_record_count: int
    normalized_exact_record_count: int
    terminal_stop_normalized_record_count: int
    publisher_x_wildcard_record_count: int
    publisher_x_wildcard_position_count: int


@dataclass(frozen=True)
class _ProteinRecord:
    identifier: str
    sequence: str
    terminal_stop: bool

    @property
    def normalized_sequence(self) -> str:
        return self.sequence[:-1] if self.terminal_stop else self.sequence


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == GZIP_MAGIC
    except OSError as error:
        raise PublishedProteinCompatibilityError(
            f"Cannot inspect protein FASTA {path.name}: {error}"
        ) from error


def _open_text(path: Path) -> TextIO:
    try:
        if _is_gzip(path):
            return gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="")
        return path.open("rt", encoding="utf-8", errors="strict", newline="")
    except OSError as error:
        raise PublishedProteinCompatibilityError(
            f"Cannot open protein FASTA {path.name}: {error}"
        ) from error


def _checksum(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                size += len(block)
                digest.update(block)
    except OSError as error:
        raise PublishedProteinCompatibilityError(
            f"Cannot checksum protein FASTA {path.name}: {error}"
        ) from error
    return size, digest.hexdigest()


def _input_signature(path: Path) -> tuple[int, int, int, int]:
    try:
        status = path.stat()
    except OSError as error:
        raise PublishedProteinCompatibilityError(
            f"Cannot stat protein FASTA {path.name}: {error}"
        ) from error
    return status.st_dev, status.st_ino, status.st_size, status.st_mtime_ns


def _validate_sequence(sequence: str, identifier: str, role: str) -> _ProteinRecord:
    if not sequence:
        raise PublishedProteinCompatibilityError(
            f"{role} protein record {identifier!r} has an empty sequence"
        )
    invalid = sorted(set(sequence) - PROTEIN_ALPHABET)
    if invalid:
        raise PublishedProteinCompatibilityError(
            f"{role} protein record {identifier!r} contains unsupported symbols: "
            + ",".join(invalid[:10])
        )
    stop_positions = [index for index, residue in enumerate(sequence) if residue == "*"]
    if stop_positions and stop_positions != [len(sequence) - 1]:
        raise PublishedProteinCompatibilityError(
            f"{role} protein record {identifier!r} contains an internal or repeated stop codon"
        )
    normalized = sequence[:-1] if stop_positions else sequence
    if not normalized:
        raise PublishedProteinCompatibilityError(
            f"{role} protein record {identifier!r} contains no amino acids"
        )
    return _ProteinRecord(identifier, sequence, bool(stop_positions))


def read_protein_fasta(path: str | Path, role: str) -> dict[str, _ProteinRecord]:
    """Read a protein FASTA without silently normalizing case or whitespace."""
    source = Path(path).expanduser()
    records: dict[str, _ProteinRecord] = {}
    identifier: str | None = None
    parts: list[str] = []

    def finish() -> None:
        nonlocal identifier, parts
        if identifier is None:
            return
        records[identifier] = _validate_sequence("".join(parts), identifier, role)

    try:
        with _open_text(source) as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith(">"):
                    finish()
                    fields = line[1:].split()
                    if not fields:
                        raise PublishedProteinCompatibilityError(
                            f"{role} protein FASTA {source.name}:{line_number} has an empty header"
                        )
                    identifier = fields[0]
                    if any(character.isspace() or ord(character) < 32 for character in identifier):
                        raise PublishedProteinCompatibilityError(
                            f"{role} protein FASTA {source.name}:{line_number} has an invalid ID"
                        )
                    if identifier in records:
                        raise PublishedProteinCompatibilityError(
                            f"{role} protein FASTA {source.name}:{line_number} repeats ID "
                            f"{identifier!r}"
                        )
                    parts = []
                elif identifier is None:
                    raise PublishedProteinCompatibilityError(
                        f"{role} protein FASTA {source.name}:{line_number} has sequence before a header"
                    )
                else:
                    if any(character.isspace() for character in line):
                        raise PublishedProteinCompatibilityError(
                            f"{role} protein record {identifier!r} contains embedded whitespace"
                        )
                    parts.append(line)
        finish()
    except (UnicodeError, EOFError, gzip.BadGzipFile) as error:
        raise PublishedProteinCompatibilityError(
            f"Cannot read {role} protein FASTA {source.name}: {error}"
        ) from error
    if not records:
        raise PublishedProteinCompatibilityError(
            f"{role} protein FASTA {source.name} contains no records"
        )
    return records


def _write_tsv(path: Path, rows: list[dict[str, object]], columns: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _compare_records(
    derived: dict[str, _ProteinRecord],
    publisher: dict[str, _ProteinRecord],
) -> list[dict[str, object]]:
    derived_ids = set(derived)
    publisher_ids = set(publisher)
    if derived_ids != publisher_ids:
        derived_only = sorted(derived_ids - publisher_ids, key=natural_key)
        publisher_only = sorted(publisher_ids - derived_ids, key=natural_key)
        details = []
        if derived_only:
            details.append(
                f"derived_only={len(derived_only)} ({','.join(derived_only[:5])})"
            )
        if publisher_only:
            details.append(
                f"publisher_only={len(publisher_only)} ({','.join(publisher_only[:5])})"
            )
        raise PublishedProteinCompatibilityError(
            "Protein ID sets are not exactly equal: " + "; ".join(details)
        )

    rows: list[dict[str, object]] = []
    failures: list[str] = []
    for identifier in sorted(derived_ids, key=natural_key):
        derived_record = derived[identifier]
        publisher_record = publisher[identifier]
        derived_sequence = derived_record.normalized_sequence
        publisher_sequence = publisher_record.normalized_sequence
        stop_difference = derived_record.terminal_stop != publisher_record.terminal_stop
        wildcard_positions: list[int] = []
        mismatch_positions: list[int] = []
        if len(derived_sequence) != len(publisher_sequence):
            status = "FAIL_LENGTH_MISMATCH"
            failures.append(
                f"{identifier}:{status}:{len(derived_sequence)}!={len(publisher_sequence)}"
            )
        else:
            for position, (derived_residue, publisher_residue) in enumerate(
                zip(derived_sequence, publisher_sequence), start=1
            ):
                if derived_residue == publisher_residue:
                    continue
                if (
                    publisher_residue == "X"
                    and derived_residue in CANONICAL_AMINO_ACIDS
                ):
                    wildcard_positions.append(position)
                else:
                    mismatch_positions.append(position)
            if mismatch_positions:
                status = "FAIL_RESIDUE_MISMATCH"
                failures.append(
                    f"{identifier}:{status}:positions="
                    + ",".join(str(value) for value in mismatch_positions[:5])
                )
            elif wildcard_positions and stop_difference:
                status = "PASS_TERMINAL_STOP_AND_PUBLISHER_X_WILDCARD"
            elif wildcard_positions:
                status = "PASS_PUBLISHER_X_WILDCARD"
            elif stop_difference:
                status = "PASS_TERMINAL_STOP_NORMALIZED"
            else:
                status = "PASS_EXACT"
        rows.append(
            {
                "protein_id": identifier,
                "derived_raw_length": len(derived_record.sequence),
                "publisher_raw_length": len(publisher_record.sequence),
                "derived_normalized_length": len(derived_sequence),
                "publisher_normalized_length": len(publisher_sequence),
                "derived_terminal_stop": str(derived_record.terminal_stop).lower(),
                "publisher_terminal_stop": str(publisher_record.terminal_stop).lower(),
                "terminal_stop_difference": str(stop_difference).lower(),
                "publisher_X_wildcard_count": len(wildcard_positions),
                "publisher_X_wildcard_positions_1based": ",".join(
                    str(value) for value in wildcard_positions
                ),
                "nonpermitted_mismatch_count": len(mismatch_positions),
                "nonpermitted_mismatch_positions_1based": ",".join(
                    str(value) for value in mismatch_positions
                ),
                "status": status,
            }
        )
    if failures:
        raise PublishedProteinCompatibilityError(
            f"Protein compatibility failed for {len(failures)} records ("
            + "; ".join(failures[:5])
            + "); no output was published"
        )
    return rows


def _validate_inputs(
    derived: Path, publisher: Path, output: Path, sample_id: str
) -> None:
    if not SAFE_SAMPLE_ID.fullmatch(sample_id):
        raise PublishedProteinCompatibilityError(
            "sample_id must start with an alphanumeric character and contain only letters, "
            "numbers, periods, underscores, or hyphens"
        )
    for role, path in (("derived", derived), ("publisher", publisher)):
        if not path.is_file() or path.stat().st_size == 0:
            raise PublishedProteinCompatibilityError(
                f"The {role} protein FASTA is missing or empty: {path.name}"
            )
    if derived.resolve() == publisher.resolve():
        raise PublishedProteinCompatibilityError(
            "Derived and publisher protein inputs must be different files"
        )
    if output.exists() or output.is_symlink():
        raise PublishedProteinCompatibilityError(
            f"Output directory already exists; refusing overwrite: {output.name}"
        )


def audit_published_protein_compatibility(
    derived_proteins: str | Path,
    publisher_proteins: str | Path,
    output_dir: str | Path,
    sample_id: str,
) -> ProteinCompatibilityResult:
    """Compare two protein sets and atomically publish an audit only on PASS."""
    derived_path = Path(derived_proteins).expanduser()
    publisher_path = Path(publisher_proteins).expanduser()
    output = Path(output_dir).expanduser()
    _validate_inputs(derived_path, publisher_path, output, sample_id)

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        input_signatures = {
            derived_path: _input_signature(derived_path),
            publisher_path: _input_signature(publisher_path),
        }
        derived_bytes, derived_sha256 = _checksum(derived_path)
        publisher_bytes, publisher_sha256 = _checksum(publisher_path)
        derived_records = read_protein_fasta(derived_path, "derived")
        publisher_records = read_protein_fasta(publisher_path, "publisher")
        record_rows = _compare_records(derived_records, publisher_records)

        exact_count = sum(row["status"] == "PASS_EXACT" for row in record_rows)
        normalized_exact_count = sum(
            int(row["publisher_X_wildcard_count"]) == 0 for row in record_rows
        )
        stop_count = sum(
            str(row["terminal_stop_difference"]) == "true" for row in record_rows
        )
        wildcard_record_count = sum(
            int(row["publisher_X_wildcard_count"]) > 0 for row in record_rows
        )
        wildcard_position_count = sum(
            int(row["publisher_X_wildcard_count"]) for row in record_rows
        )
        summary_rows = [
            {
                "sample_id": sample_id,
                "status": "PASS",
                "publication_gate": "PASS",
                "workflow_version": WORKFLOW_VERSION,
                "derived_record_count": len(derived_records),
                "publisher_record_count": len(publisher_records),
                "exact_ID_set": "true",
                "exact_record_count": exact_count,
                "normalized_exact_record_count": normalized_exact_count,
                "terminal_stop_normalized_record_count": stop_count,
                "publisher_X_wildcard_record_count": wildcard_record_count,
                "publisher_X_wildcard_position_count": wildcard_position_count,
                "nonpermitted_mismatch_record_count": 0,
                "derived_input_bytes": derived_bytes,
                "derived_input_sha256": derived_sha256,
                "publisher_input_bytes": publisher_bytes,
                "publisher_input_sha256": publisher_sha256,
            }
        ]
        summary_name = f"{sample_id}.published_protein_compatibility.summary.tsv"
        records_name = f"{sample_id}.published_protein_compatibility.records.tsv"
        _write_tsv(
            staging / summary_name,
            summary_rows,
            (
                "sample_id", "status", "publication_gate", "workflow_version",
                "derived_record_count", "publisher_record_count", "exact_ID_set",
                "exact_record_count", "normalized_exact_record_count",
                "terminal_stop_normalized_record_count",
                "publisher_X_wildcard_record_count", "publisher_X_wildcard_position_count",
                "nonpermitted_mismatch_record_count", "derived_input_bytes",
                "derived_input_sha256", "publisher_input_bytes", "publisher_input_sha256",
            ),
        )
        _write_tsv(
            staging / records_name,
            record_rows,
            (
                "protein_id", "derived_raw_length", "publisher_raw_length",
                "derived_normalized_length", "publisher_normalized_length",
                "derived_terminal_stop", "publisher_terminal_stop",
                "terminal_stop_difference", "publisher_X_wildcard_count",
                "publisher_X_wildcard_positions_1based", "nonpermitted_mismatch_count",
                "nonpermitted_mismatch_positions_1based", "status",
            ),
        )
        manifest = {
            "schema_version": 1,
            "workflow": "published_protein_compatibility",
            "workflow_version": WORKFLOW_VERSION,
            "status": "PASS",
            "publication_gate": "PASS",
            "sample_id": sample_id,
            "execution": {"processes": 1, "worker_threads": 0},
            "inputs": [
                {
                    "role": "derived_primary_proteins",
                    "file_name": derived_path.name,
                    "bytes": derived_bytes,
                    "sha256": derived_sha256,
                },
                {
                    "role": "publisher_proteins",
                    "file_name": publisher_path.name,
                    "bytes": publisher_bytes,
                    "sha256": publisher_sha256,
                },
            ],
            "policy": {
                "identifier_policy": "exact_first_token_ID_set_equality",
                "case_policy": "no_case_normalization",
                "terminal_stop_policy": "remove_at_most_one_terminal_asterisk_per_sequence",
                "internal_stop_policy": "reject_record",
                "publisher_X_policy": (
                    "one_way_wildcard_only_against_one_unambiguous_canonical_derived_amino_acid"
                ),
                "all_other_residue_differences": "reject_complete_run",
            },
            "counts": {
                "compared_records": len(record_rows),
                "exact_records": exact_count,
                "normalized_exact_records": normalized_exact_count,
                "terminal_stop_normalized_records": stop_count,
                "publisher_X_wildcard_records": wildcard_record_count,
                "publisher_X_wildcard_positions": wildcard_position_count,
                "nonpermitted_mismatch_records": 0,
            },
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

        for path, signature in input_signatures.items():
            if _input_signature(path) != signature:
                raise PublishedProteinCompatibilityError(
                    f"Input {path.name} changed during comparison; no output was published"
                )
        checksum_rows = []
        for path in sorted(staging.iterdir(), key=lambda item: item.name):
            if path.is_file() and path.name != "checksums.tsv":
                size, digest = _checksum(path)
                checksum_rows.append({"file": path.name, "bytes": size, "sha256": digest})
        _write_tsv(
            staging / "checksums.tsv",
            checksum_rows,
            ("file", "bytes", "sha256"),
        )
        if output.exists() or output.is_symlink():
            raise PublishedProteinCompatibilityError(
                f"Output directory appeared during the run; refusing overwrite: {output.name}"
            )
        os.replace(staging, output)
        return ProteinCompatibilityResult(
            output_dir=output,
            record_count=len(record_rows),
            exact_record_count=exact_count,
            normalized_exact_record_count=normalized_exact_count,
            terminal_stop_normalized_record_count=stop_count,
            publisher_x_wildcard_record_count=wildcard_record_count,
            publisher_x_wildcard_position_count=wildcard_position_count,
        )
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
