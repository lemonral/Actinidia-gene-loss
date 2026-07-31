#!/usr/bin/env python3
"""Merge disjoint legacy and new expected-deletion-locus spatial inputs."""

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


SOURCE_COLUMNS = ("source_id", "input_dir", "expected_unit_count")
CALL_COLUMNS = ("assembly_unit_id", "reference_gene_id", "classification", "callable")
COORDINATE_COLUMNS = (
    "assembly_unit_id", "reference_gene_id", "classification", "chromosome",
    "expected_locus_start_1based", "expected_locus_end_1based", "coordinate_semantics",
)
ASSEMBLY_COLUMNS = (
    "assembly_unit_id", "biological_species", "haplotype_or_subgenome", "assembly_scope",
    "genome", "gff", "genome_local_sha256", "gff_local_sha256",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.resolve().open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    source = path.resolve()
    if not source.is_file() or (not allow_empty and source.stat().st_size == 0):
        raise MergeError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def read_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise MergeError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise MergeError(f"{path}: JSON root is not an object")
    return value


def read_rows(path: Path, columns: tuple[str, ...], *, allow_empty: bool = False) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != columns:
            raise MergeError(f"{path.name}: columns differ from exact schema")
        rows = list(reader)
    if not allow_empty and not rows:
        raise MergeError(f"{path.name}: no data rows")
    return rows


def read_sources(path: Path) -> list[dict[str, str]]:
    rows = read_rows(path, SOURCE_COLUMNS)
    ids = [row["source_id"] for row in rows]
    if any(not value for row in rows for value in row.values()) or len(ids) != len(set(ids)):
        raise MergeError("source manifest has empty/duplicate values")
    for row in rows:
        if int(row["expected_unit_count"]) < 1:
            raise MergeError("expected unit counts must be positive")
    return rows


def validate_checksums(root: Path) -> None:
    rows = read_rows(root / "checksums.tsv", ("file", "bytes", "sha256"))
    seen: set[str] = set()
    for row in rows:
        name = row["file"]
        if Path(name).name != name or name in seen:
            raise MergeError(f"{root.name}: unsafe/duplicate checksum row")
        seen.add(name)
        observed = binding(root / name, allow_empty=True)
        if observed["bytes"] != int(row["bytes"]) or observed["sha256"] != row["sha256"]:
            raise MergeError(f"{root.name}: checksum mismatch for {name}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--expected-total-units", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    staging: Path | None = None
    try:
        root = args.data_root.resolve()
        sources = read_sources(args.sources)
        if sum(int(row["expected_unit_count"]) for row in sources) != args.expected_total_units:
            raise MergeError("source unit counts do not close to expected total")
        output = args.output_dir.resolve()
        if output.exists():
            raise MergeError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        calls: list[dict[str, str]] = []
        coordinates: list[dict[str, str]] = []
        assemblies: list[dict[str, str]] = []
        source_audits: list[dict[str, object]] = []
        all_units: set[str] = set()
        all_call_keys: set[tuple[str, str]] = set()
        for source in sources:
            relative = Path(source["input_dir"])
            if relative.is_absolute() or ".." in relative.parts:
                raise MergeError("unsafe source input directory")
            source_root = (root / relative).absolute()
            if not source_root.is_relative_to(root) or not source_root.is_dir():
                raise MergeError(f"missing source input directory: {relative}")
            validate_checksums(source_root)
            report_path = source_root / "run_manifest.json"
            report = read_json(report_path)
            allowed = {
                "callable_positive_deleted_expected_locus_spatial_inputs",
                "legacy_conservative_deleted_expected_locus_spatial_inputs",
            }
            expected_units = int(source["expected_unit_count"])
            if (
                report.get("status") != "PASS"
                or report.get("workflow") not in allowed
                or report.get("unit_count") != expected_units
            ):
                raise MergeError(f"{source['source_id']}: source report is not exact PASS")
            outputs = report.get("outputs")
            if not isinstance(outputs, dict):
                raise MergeError(f"{source['source_id']}: missing report outputs")
            call_path = source_root / "positive_deleted_calls.tsv"
            coordinate_path = source_root / "expected_deleted_locus_coordinates.tsv"
            assembly_path = source_root / "assembly_manifest.tsv"
            for key, path in (
                ("positive_calls", call_path),
                ("expected_locus_coordinates", coordinate_path),
                ("assembly_manifest", assembly_path),
            ):
                if outputs.get(key) != binding(path, allow_empty=key != "assembly_manifest"):
                    raise MergeError(f"{source['source_id']}: output binding failed for {key}")
            source_calls = read_rows(call_path, CALL_COLUMNS, allow_empty=True)
            source_coordinates = read_rows(coordinate_path, COORDINATE_COLUMNS, allow_empty=True)
            source_assemblies = read_rows(assembly_path, ASSEMBLY_COLUMNS)
            source_units = {row["assembly_unit_id"] for row in source_assemblies}
            if (
                len(source_units) != expected_units
                or len(source_units) != len(source_assemblies)
                or source_units.intersection(all_units)
            ):
                raise MergeError(f"{source['source_id']}: unit scope/disjointness failed")
            source_call_keys = {(row["assembly_unit_id"], row["reference_gene_id"]) for row in source_calls}
            source_coordinate_keys = {
                (row["assembly_unit_id"], row["reference_gene_id"]) for row in source_coordinates
            }
            if (
                len(source_call_keys) != len(source_calls)
                or len(source_coordinate_keys) != len(source_coordinates)
                or source_call_keys != source_coordinate_keys
                or all_call_keys.intersection(source_call_keys)
            ):
                raise MergeError(f"{source['source_id']}: call/coordinate closure failed")
            if any(
                row["assembly_unit_id"] not in source_units
                or row["classification"] != "positive_deleted"
                or row["callable"] != "true"
                for row in source_calls
            ):
                raise MergeError(f"{source['source_id']}: invalid positive call row")
            for assembly in source_assemblies:
                rewritten = dict(assembly)
                for field, checksum_field in (
                    ("genome", "genome_local_sha256"), ("gff", "gff_local_sha256")
                ):
                    asset = Path(assembly[field])
                    if not asset.is_absolute():
                        asset = (assembly_path.parent / asset).resolve()
                    else:
                        asset = asset.resolve()
                    if not asset.is_file() or sha256(asset) != assembly[checksum_field]:
                        raise MergeError(f"{source['source_id']}: assembly asset binding failed")
                    rewritten[field] = os.path.relpath(asset, staging)
                assemblies.append(rewritten)
            calls.extend(source_calls)
            coordinates.extend(source_coordinates)
            all_units.update(source_units)
            all_call_keys.update(source_call_keys)
            source_audits.append(
                {
                    "source_id": source["source_id"], "input_dir": source["input_dir"],
                    "unit_count": len(source_units), "positive_deleted_count": len(source_calls),
                    "run_manifest": binding(report_path), "checksums": binding(source_root / "checksums.tsv"),
                }
            )
        if len(all_units) != args.expected_total_units:
            raise MergeError("merged unit count does not close")
        calls.sort(key=lambda row: (row["assembly_unit_id"], row["reference_gene_id"]))
        coordinates.sort(key=lambda row: (row["assembly_unit_id"], row["reference_gene_id"]))
        assemblies.sort(key=lambda row: row["assembly_unit_id"])
        call_out = staging / "positive_deleted_calls.tsv"
        coordinate_out = staging / "expected_deleted_locus_coordinates.tsv"
        assembly_out = staging / "assembly_manifest.tsv"
        for path, rows, columns in (
            (call_out, calls, CALL_COLUMNS),
            (coordinate_out, coordinates, COORDINATE_COLUMNS),
            (assembly_out, assemblies, ASSEMBLY_COLUMNS),
        ):
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        report = {
            "schema_version": 1, "workflow": "merged_primary_expected_deleted_locus_spatial_inputs",
            "status": "PASS", "unit_count": len(all_units), "positive_deleted_count": len(calls),
            "coordinate_semantics": "expected bilateral-SynOrths-bounded interval midpoint; not an observed remnant",
            "source_manifest": binding(args.sources), "sources": source_audits,
            "outputs": {
                "positive_calls": binding(call_out, allow_empty=True),
                "expected_locus_coordinates": binding(coordinate_out, allow_empty=True),
                "assembly_manifest": binding(assembly_out),
            },
        }
        (staging / "run_manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        with (staging / "checksums.tsv").open("w", encoding="utf-8") as handle:
            handle.write("file\tbytes\tsha256\n")
            for path in sorted(staging.iterdir()):
                if path.is_file() and path.name != "checksums.tsv":
                    item = binding(path, allow_empty=True)
                    handle.write(f"{path.name}\t{item['bytes']}\t{item['sha256']}\n")
        os.rename(staging, output)
        staging = None
        print(json.dumps({"status": "PASS", "units": len(all_units), "calls": len(calls)}))
        return 0
    except (MergeError, OSError, UnicodeError, csv.Error, ValueError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    finally:
        if staging is not None and staging.exists():
            shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
