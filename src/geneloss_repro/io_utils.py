"""Small, dependency-free I/O and validation helpers.

The analysis tables are written as UTF-8 TSV rather than Excel by default.
TSV keeps a stable schema, is suitable for Git, and avoids an implicit Excel
engine dependency.  Excel can be made downstream from these tables if needed.
"""

from __future__ import annotations

import csv
import hashlib
import os
import re
import tempfile
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any


class SchemaError(ValueError):
    """Raised when a tabular input cannot be interpreted safely."""


def ensure_parent(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def sha256_file(path: str | Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: str | Path, required: Sequence[str] = ()) -> list[dict[str, str]]:
    """Read a headered tab-separated file and validate required columns."""
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SchemaError(f"{path}: missing header row")
        fields = [field.strip() for field in reader.fieldnames]
        missing = [field for field in required if field not in fields]
        if missing:
            raise SchemaError(
                f"{path}: required column(s) missing: {', '.join(missing)}; "
                f"found: {', '.join(fields)}"
            )
        return [
            {key.strip(): "" if value is None else value.strip() for key, value in row.items()}
            for row in reader
            if any(value not in (None, "") for value in row.values())
        ]


def write_tsv(
    path: str | Path,
    rows: Iterable[Mapping[str, Any]],
    fieldnames: Sequence[str],
) -> Path:
    """Atomically write a headered UTF-8 TSV with a fixed column order."""
    output = ensure_parent(path)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", delete=False, dir=output.parent,
        prefix=f".{output.name}.", suffix=".tmp"
    ) as temporary:
        writer = csv.DictWriter(
            temporary,
            fieldnames=list(fieldnames),
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in fieldnames})
        temporary_name = temporary.name
    os.replace(temporary_name, output)
    return output


def concatenate_tsv(paths: Sequence[str | Path], output_path: str | Path) -> Path:
    """Concatenate headered TSVs only when their schemas are exactly identical."""
    if not paths:
        raise SchemaError("at least one TSV input is required")
    expected_fields: list[str] | None = None
    combined: list[dict[str, str]] = []
    for path in paths:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            if not fields:
                raise SchemaError(f"{path}: missing TSV header")
            if expected_fields is None:
                expected_fields = fields
            elif fields != expected_fields:
                raise SchemaError(
                    f"{path}: schema differs from first input; expected {expected_fields}, found {fields}. "
                    "Do not mix workflow versions in one combined table."
                )
            combined.extend(
                {key: "" if value is None else value for key, value in row.items()}
                for row in reader if any(value not in (None, "") for value in row.values())
            )
    return write_tsv(output_path, combined, expected_fields or [])


def read_id_file(path: str | Path) -> set[str]:
    """Read first-column IDs from a one-column, TSV, CSV, or plain-text file.

    A header named ``reference_gene``, ``gene_id``, or ``query_id`` is skipped.
    Blank and comment lines are ignored.
    """
    identifiers: set[str] = set()
    header_values = {"reference_gene", "gene_id", "query_id", "lost_gene", "id"}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            first = re.split(r"[\t,\s]+", line, maxsplit=1)[0].strip()
            if first.lower() in header_values:
                continue
            identifiers.add(first)
    return identifiers


def natural_key(value: str) -> list[object]:
    """Sort chromosome labels naturally (Chr2 before Chr10)."""
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value)]


def parse_float(value: str, field: str, source: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise SchemaError(f"{source}: invalid {field!r} value {value!r}") from exc


def parse_int(value: str, field: str, source: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError) as exc:
        raise SchemaError(f"{source}: invalid {field!r} value {value!r}") from exc


def bh_adjust(p_values: Sequence[float | None]) -> list[float | None]:
    """Benjamini-Hochberg adjustment preserving ``None`` values and order."""
    indexed = [(index, value) for index, value in enumerate(p_values) if value is not None]
    if not indexed:
        return [None for _ in p_values]
    indexed.sort(key=lambda pair: pair[1])
    count = len(indexed)
    adjusted_sorted: list[float] = [0.0] * count
    running = 1.0
    for rank in range(count, 0, -1):
        _, p_value = indexed[rank - 1]
        candidate = min(1.0, p_value * count / rank)
        running = min(running, candidate)
        adjusted_sorted[rank - 1] = running
    adjusted: list[float | None] = [None for _ in p_values]
    for (original_index, _), value in zip(indexed, adjusted_sorted):
        adjusted[original_index] = value
    return adjusted


def format_number(value: float | int | None, digits: int = 12) -> str:
    if value is None:
        return ""
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}g}"
