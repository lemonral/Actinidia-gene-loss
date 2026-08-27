"""Write complete, privacy-safe publication figure bundles.

A bundle is assembled in a private staging directory and published with one
directory rename.  It contains both graphical formats, the exact table passed
to the plotting code, an English caption, machine-readable validation, and a
checksum manifest.  The manifest deliberately records basenames only: local
paths, environment details, and shell commands are not accepted as metadata.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import re
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any


_SIMPLE_BASENAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,99}\Z")
_PUBLIC_FILE_MODE = 0o644
_PUBLIC_DIRECTORY_MODE = 0o755


@dataclass(frozen=True)
class FigureBundle:
    """Paths in a successfully published figure bundle."""

    directory: Path
    png: Path
    pdf: Path
    plot_data: Path
    caption: Path
    validation: Path
    manifest: Path


def _validate_basename(value: str) -> str:
    if not isinstance(value, str) or not _SIMPLE_BASENAME.fullmatch(value):
        raise ValueError(
            "basename must contain only ASCII letters, digits, underscores, or "
            "hyphens; it must start with a letter or digit and be at most 100 characters"
        )
    return value


def _validate_columns(columns: Sequence[str], label: str) -> tuple[str, ...]:
    if isinstance(columns, (str, bytes)) or not columns:
        raise ValueError(f"{label} must be a nonempty sequence of column names")
    normalized = tuple(columns)
    if any(
        not isinstance(column, str)
        or not column
        or any(character in column for character in "\t\r\n")
        for column in normalized
    ):
        raise ValueError(f"{label} contains an empty or invalid column name")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{label} contains duplicate column names")
    return normalized


def _validate_scalar(value: Any, location: str) -> None:
    if value is None or isinstance(value, (str, int, bool)):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{location} contains a non-finite number")
        return
    raise TypeError(
        f"{location} contains unsupported value type {type(value).__name__}; "
        "use strings, finite numbers, booleans, or null"
    )


def _materialize_rows(
    rows: Sequence[Mapping[str, Any]],
    columns: Sequence[str],
    label: str,
) -> list[dict[str, Any]]:
    if isinstance(rows, (str, bytes, Mapping)) or not isinstance(rows, Sequence) or not rows:
        raise ValueError(f"{label} must be a nonempty sequence of row mappings")
    expected = set(columns)
    materialized: list[dict[str, Any]] = []
    for index, row in enumerate(rows, start=1):
        if not isinstance(row, Mapping):
            raise TypeError(f"{label} row {index} is not a mapping")
        actual = set(row)
        if actual != expected:
            missing = sorted(expected - actual)
            extra = sorted(actual - expected)
            raise ValueError(
                f"{label} row {index} does not exactly match its schema; "
                f"missing={missing}, extra={extra}"
            )
        copied: dict[str, Any] = {}
        for column in columns:
            value = row[column]
            _validate_scalar(value, f"{label} row {index}, column {column!r}")
            copied[column] = value
        materialized.append(copied)
    return materialized


def _tsv_bytes(rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=list(columns),
        delimiter="\t",
        extrasaction="raise",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def _validate_json(value: Any, location: str = "validation") -> None:
    if isinstance(value, Mapping):
        for key, nested in value.items():
            if not isinstance(key, str):
                raise TypeError(f"{location} contains a non-string object key")
            _validate_json(nested, f"{location}.{key}")
        return
    if isinstance(value, (list, tuple)):
        for index, nested in enumerate(value):
            _validate_json(nested, f"{location}[{index}]")
        return
    _validate_scalar(value, location)


def _atomic_write(path: Path, content: bytes) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary_name, _PUBLIC_FILE_MODE)
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


@lru_cache(maxsize=1)
def _register_arial_fonts() -> tuple[str, ...]:
    """Register Arial faces and fail rather than silently substituting a font.

    Set ``ARIAL_FONT_DIR`` when Arial is installed outside a standard system
    font directory. Arial itself is not redistributed by this repository.
    """

    from matplotlib import font_manager, rcParams

    faces = (
        ("normal", "normal"),
        ("italic", "normal"),
        ("normal", "bold"),
        ("italic", "bold"),
    )

    def resolve_faces() -> tuple[str, ...] | None:
        resolved: list[str] = []
        for style, weight in faces:
            properties = font_manager.FontProperties(
                family="Arial",
                style=style,
                weight=weight,
            )
            try:
                path = font_manager.findfont(properties, fallback_to_default=False)
            except ValueError:
                return None
            if font_manager.FontProperties(fname=path).get_name() != "Arial":
                return None
            resolved.append(path)
        return tuple(resolved)

    resolved = resolve_faces()
    if resolved is None:
        configured = os.environ.get("ARIAL_FONT_DIR", "").strip()
        directories = [
            Path(configured) if configured else None,
            Path("/usr/share/fonts/truetype/msttcorefonts"),
            Path("/usr/local/share/fonts"),
            Path("/System/Library/Fonts/Supplemental"),
            Path("/Library/Fonts"),
        ]
        filenames = (
            "Arial.ttf",
            "Arial Bold.ttf",
            "Arial Italic.ttf",
            "Arial Bold Italic.ttf",
            "arial.ttf",
            "arialbd.ttf",
            "ariali.ttf",
            "arialbi.ttf",
            "Arial_Bold.ttf",
            "Arial_Italic.ttf",
            "Arial_Bold_Italic.ttf",
        )
        for directory in directories:
            if directory is None or not directory.is_dir():
                continue
            for filename in filenames:
                candidate = directory / filename
                if candidate.is_file():
                    font_manager.fontManager.addfont(str(candidate))
        resolved = resolve_faces()

    if resolved is None:
        raise RuntimeError(
            "Arial regular, italic, bold, and bold-italic faces are required. "
            "Install Arial or set ARIAL_FONT_DIR to the directory containing them."
        )

    rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "mathtext.cal": "Arial:italic",
            "mathtext.sf": "Arial",
            "mathtext.tt": "Arial",
            "mathtext.fallback": None,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    return resolved


def _atomic_save_figure(figure: Any, path: Path, image_format: str, dpi: int) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=f".{image_format}.tmp"
    )
    os.close(descriptor)
    try:
        figure.savefig(
            temporary_name,
            format=image_format,
            dpi=dpi,
            bbox_inches="tight",
        )
        temporary = Path(temporary_name)
        if not temporary.is_file() or temporary.stat().st_size == 0:
            raise RuntimeError(f"matplotlib produced an empty {image_format.upper()} file")
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        os.chmod(temporary, _PUBLIC_FILE_MODE)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _checksum_entry(path: Path, role: str) -> dict[str, str | int]:
    return {
        "role": role,
        "basename": path.name,
        "sha256": _sha256(path),
        "bytes": path.stat().st_size,
    }


def _validated_inputs(input_paths: Sequence[str | Path]) -> list[dict[str, str | int]]:
    if isinstance(input_paths, (str, bytes, Path)):
        raise TypeError("input_paths must be a sequence of file paths, not one path")
    entries: list[dict[str, str | int]] = []
    seen: set[str] = set()
    for raw_path in input_paths:
        path = Path(raw_path)
        if not path.is_file():
            raise FileNotFoundError(f"figure input is not a regular file: {path}")
        name = path.name
        if name in seen:
            raise ValueError(f"input basenames must be unique in a privacy-safe manifest: {name}")
        seen.add(name)
        entries.append(
            {
                "basename": name,
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
        )
    return sorted(entries, key=lambda entry: str(entry["basename"]))


def _prepare_destination(output: Path) -> bool:
    """Validate destination and return whether an empty directory already exists."""
    if output.is_symlink():
        raise ValueError(f"refusing symlink output directory: {output}")
    if not output.exists():
        return False
    if not output.is_dir():
        raise ValueError(f"output path exists and is not a directory: {output}")
    if any(output.iterdir()):
        raise FileExistsError(f"refusing nonempty output directory: {output}")
    return True


def write_figure_bundle(
    *,
    figure: Any,
    output_dir: str | Path,
    basename: str,
    plot_rows: Sequence[Mapping[str, Any]],
    plot_columns: Sequence[str],
    caption: str,
    validation: Mapping[str, Any] | Sequence[Mapping[str, Any]],
    validation_format: str = "json",
    validation_columns: Sequence[str] | None = None,
    input_paths: Sequence[str | Path] = (),
    dpi: int = 300,
) -> FigureBundle:
    """Atomically publish one complete PNG/PDF figure bundle.

    ``plot_rows`` must be the exact rows plotted, in their plotting order.  A
    JSON validation object is used by default.  For TSV validation, provide a
    nonempty sequence of exact row mappings and ``validation_columns``.

    The destination may be absent or an existing empty directory.  A nonempty
    directory, symlink, schema mismatch, missing input, non-finite value, or
    partial figure save aborts publication and removes all staged files.
    """
    name = _validate_basename(basename)
    if not callable(getattr(figure, "savefig", None)):
        raise TypeError("figure must be a matplotlib Figure-like object with savefig()")
    if not isinstance(dpi, int) or isinstance(dpi, bool) or not 72 <= dpi <= 2400:
        raise ValueError("dpi must be an integer from 72 through 2400")
    if not isinstance(caption, str) or not caption.strip():
        raise ValueError("caption must be nonempty English text")
    if "\x00" in caption:
        raise ValueError("caption contains a NUL character")

    data_columns = _validate_columns(plot_columns, "plot_columns")
    data_rows = _materialize_rows(plot_rows, data_columns, "plot_rows")
    input_entries = _validated_inputs(input_paths)

    normalized_validation_format = validation_format.lower()
    if normalized_validation_format == "json":
        if validation_columns is not None:
            raise ValueError("validation_columns is only valid for TSV validation")
        _validate_json(validation)
        validation_content = (
            json.dumps(
                validation,
                ensure_ascii=False,
                allow_nan=False,
                indent=2,
                sort_keys=True,
            )
            + "\n"
        ).encode("utf-8")
    elif normalized_validation_format == "tsv":
        if validation_columns is None:
            raise ValueError("validation_columns is required for TSV validation")
        columns = _validate_columns(validation_columns, "validation_columns")
        rows = _materialize_rows(validation, columns, "validation")  # type: ignore[arg-type]
        validation_content = _tsv_bytes(rows, columns)
    else:
        raise ValueError("validation_format must be 'json' or 'tsv'")

    output = Path(output_dir)
    existed_empty = _prepare_destination(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.{name}.", suffix=".tmp", dir=output.parent)
    )
    os.chmod(staging, _PUBLIC_DIRECTORY_MODE)

    png_name = f"{name}.png"
    pdf_name = f"{name}.pdf"
    data_name = f"{name}.plot_data.tsv"
    caption_name = f"{name}.caption.txt"
    validation_name = f"{name}.validation.{normalized_validation_format}"
    manifest_name = f"{name}.manifest.json"
    try:
        png = staging / png_name
        pdf = staging / pdf_name
        plot_data = staging / data_name
        caption_path = staging / caption_name
        validation_path = staging / validation_name
        manifest_path = staging / manifest_name

        _atomic_save_figure(figure, png, "png", dpi)
        _atomic_save_figure(figure, pdf, "pdf", dpi)
        _atomic_write(plot_data, _tsv_bytes(data_rows, data_columns))
        _atomic_write(caption_path, (caption.rstrip() + "\n").encode("utf-8"))
        _atomic_write(validation_path, validation_content)

        outputs = [
            _checksum_entry(png, "figure_png"),
            _checksum_entry(pdf, "figure_pdf"),
            _checksum_entry(plot_data, "plot_data_tsv"),
            _checksum_entry(caption_path, "caption_text"),
            _checksum_entry(validation_path, f"validation_{normalized_validation_format}"),
        ]
        manifest = {
            "schema_version": "1.0",
            "bundle_basename": name,
            "inputs": input_entries,
            "outputs": outputs,
        }
        _atomic_write(
            manifest_path,
            (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"),
        )

        if existed_empty:
            output.rmdir()
        os.replace(staging, output)
        os.chmod(output, _PUBLIC_DIRECTORY_MODE)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    return FigureBundle(
        directory=output,
        png=output / png_name,
        pdf=output / pdf_name,
        plot_data=output / data_name,
        caption=output / caption_name,
        validation=output / validation_name,
        manifest=output / manifest_name,
    )
