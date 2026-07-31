#!/usr/bin/env python3
"""Fail-closed merge of complete, disjoint assembly-unit loss matrices."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path


class MergeError(RuntimeError):
    pass


SOURCE_COLUMNS = ("source_id", "matrix_dir", "expected_unit_count")
MATRIX_COLUMNS = (
    "reference_gene_id",
    "assembly_unit_id",
    "classification",
    "callable",
    "evidence_source",
    "primary_search_state",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    source = path.resolve()
    if not source.is_file() or (not allow_empty and source.stat().st_size == 0):
        raise MergeError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def resolve(root: Path, value: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise MergeError(f"unsafe data-root-relative path: {value!r}")
    path = (root / relative).absolute()
    if not path.is_relative_to(root) or not path.is_dir():
        raise MergeError(f"missing or unsafe matrix directory: {value!r}")
    return path


def read_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise MergeError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise MergeError(f"{path}: JSON root is not an object")
    return value


def read_sources(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SOURCE_COLUMNS:
            raise MergeError("source manifest columns differ from exact schema")
        rows = list(reader)
    ids = [row["source_id"] for row in rows]
    if not rows or any(not value for row in rows for value in row.values()) or len(ids) != len(set(ids)):
        raise MergeError("source manifest is empty or has empty/duplicate values")
    for row in rows:
        try:
            count = int(row["expected_unit_count"])
        except ValueError as error:
            raise MergeError("expected_unit_count is not an integer") from error
        if count < 1:
            raise MergeError("expected_unit_count must be positive")
    return rows


def validate_checksums(root: Path) -> dict[str, dict[str, object]]:
    path = root / "checksums.tsv"
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("file", "bytes", "sha256"):
            raise MergeError(f"{root.name}: checksum columns differ from exact schema")
        rows = list(reader)
    result: dict[str, dict[str, object]] = {}
    for row in rows:
        name = row["file"]
        if not name or Path(name).name != name or name in result:
            raise MergeError(f"{root.name}: unsafe/duplicate checksum filename")
        file_path = root / name
        observed = binding(file_path, allow_empty=True)
        try:
            expected_size = int(row["bytes"])
        except ValueError as error:
            raise MergeError(f"{root.name}: invalid checksum byte count") from error
        if observed["bytes"] != expected_size or observed["sha256"] != row["sha256"]:
            raise MergeError(f"{root.name}: checksum mismatch for {name}")
        result[name] = observed
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--expected-total-units", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    staging: Path | None = None
    try:
        if args.expected_total_units < 1:
            raise MergeError("expected-total-units must be positive")
        root = args.data_root.resolve()
        sources = read_sources(args.sources)
        if sum(int(row["expected_unit_count"]) for row in sources) != args.expected_total_units:
            raise MergeError("source expected-unit counts do not close to expected-total-units")
        output = args.output_dir.absolute()
        if output.exists():
            raise MergeError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        matrix_output = staging / "complete_unit_loss_matrix.tsv"
        all_pairs: set[tuple[str, str]] = set()
        all_units: set[str] = set()
        reference_ids: set[str] | None = None
        source_audits: list[dict[str, object]] = []
        total_rows = 0
        with matrix_output.open("w", encoding="utf-8", newline="") as output_handle:
            writer = csv.DictWriter(
                output_handle, fieldnames=MATRIX_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for source in sources:
                source_id = source["source_id"]
                matrix_root = resolve(root, source["matrix_dir"])
                checksums = validate_checksums(matrix_root)
                report_path = matrix_root / "run_manifest.json"
                matrix_path = matrix_root / "complete_unit_loss_matrix.tsv"
                if "run_manifest.json" not in checksums or "complete_unit_loss_matrix.tsv" not in checksums:
                    raise MergeError(f"{source_id}: required files absent from checksum table")
                report = read_json(report_path)
                expected_units = int(source["expected_unit_count"])
                if report.get("status") != "PASS" or report.get("assembly_unit_count") != expected_units:
                    raise MergeError(f"{source_id}: source matrix report is not exact PASS")
                output_binding = report.get("outputs")
                if not isinstance(output_binding, dict) or output_binding.get(
                    "complete_unit_loss_matrix"
                ) != binding(matrix_path):
                    raise MergeError(f"{source_id}: source report/output binding failed")
                source_units: set[str] = set()
                source_genes: set[str] = set()
                source_rows = 0
                with matrix_path.open("r", encoding="utf-8", newline="") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    if tuple(reader.fieldnames or ()) != MATRIX_COLUMNS:
                        raise MergeError(f"{source_id}: source matrix schema differs")
                    for line_number, row in enumerate(reader, 2):
                        gene = row["reference_gene_id"]
                        unit = row["assembly_unit_id"]
                        if not gene or not unit:
                            raise MergeError(f"{source_id}:{line_number}: empty gene/unit")
                        pair = (gene, unit)
                        if pair in all_pairs:
                            raise MergeError(f"{source_id}:{line_number}: duplicate unit/gene pair")
                        all_pairs.add(pair)
                        source_units.add(unit)
                        source_genes.add(gene)
                        writer.writerow(row)
                        source_rows += 1
                        total_rows += 1
                if len(source_units) != expected_units or all_units.intersection(source_units):
                    raise MergeError(f"{source_id}: unit count or cross-source disjointness failed")
                if source_rows != len(source_units) * len(source_genes):
                    raise MergeError(f"{source_id}: source unit-by-gene grid is incomplete")
                if reference_ids is None:
                    reference_ids = source_genes
                elif source_genes != reference_ids:
                    raise MergeError(f"{source_id}: reference-gene universe differs")
                all_units.update(source_units)
                source_audits.append(
                    {
                        "source_id": source_id,
                        "matrix_directory": source["matrix_dir"],
                        "assembly_unit_count": len(source_units),
                        "reference_gene_count": len(source_genes),
                        "matrix_rows": source_rows,
                        "run_manifest": binding(report_path),
                        "complete_unit_loss_matrix": binding(matrix_path),
                        "checksums": binding(matrix_root / "checksums.tsv"),
                    }
                )
        if reference_ids is None or len(all_units) != args.expected_total_units:
            raise MergeError("merged unit/reference closure failed")
        expected_rows = len(reference_ids) * len(all_units)
        if total_rows != expected_rows or len(all_pairs) != expected_rows:
            raise MergeError("merged unit-by-gene grid is incomplete")
        report = {
            "schema_version": 1,
            "workflow": "merged_primary_complete_loss_matrix",
            "status": "PASS",
            "assembly_unit_count": len(all_units),
            "reference_gene_count": len(reference_ids),
            "matrix_rows": total_rows,
            "expected_matrix_rows": expected_rows,
            "assembly_units": sorted(all_units),
            "source_manifest": binding(args.sources),
            "sources": source_audits,
            "outputs": {"complete_unit_loss_matrix": binding(matrix_output)},
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        with (staging / "checksums.tsv").open("w", encoding="utf-8", newline="") as handle:
            handle.write("file\tbytes\tsha256\n")
            for path in sorted(staging.iterdir()):
                if path.is_file() and path.name != "checksums.tsv":
                    item = binding(path, allow_empty=True)
                    handle.write(f"{path.name}\t{item['bytes']}\t{item['sha256']}\n")
        os.rename(staging, output)
        staging = None
        print(json.dumps({"status": "PASS", "units": len(all_units), "rows": total_rows}, sort_keys=True))
        return 0
    except (MergeError, OSError, UnicodeError, csv.Error, ValueError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    finally:
        if staging is not None and staging.exists():
            shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
