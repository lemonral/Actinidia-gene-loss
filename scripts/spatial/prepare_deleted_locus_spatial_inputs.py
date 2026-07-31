#!/usr/bin/env python3
"""Build exact spatial inputs for callable positive-deletion loci.

The coordinate of a deleted gene is not an observed fragment.  It is the
midpoint-bearing interval bounded by the two accepted SynOrths flanks used by
the translated-search workflow.  This adapter keeps that distinction explicit,
maps publisher chromosome IDs to the HY4A-compatible Chr01-Chr29 labels, and
binds every source and generated table by SHA-256.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path


class InputError(RuntimeError):
    pass


MANIFEST_COLUMNS = (
    "unit",
    "biological_species",
    "haplotype_or_subgenome",
    "assembly_scope",
    "search_dir",
    "label_map",
    "relabel_dir",
)
LABEL_COLUMNS = (
    "query_chromosome",
    "final_chromosome",
    "coordinate_reference",
    "assignment_method",
    "assigned_score",
    "reciprocal_coverage",
    "orientation_to_hy4a",
    "hy4p_and_jcvi_agree",
    "strict_homology_gates_pass",
    "confidence_flag",
)
STATE_COLUMNS = (
    "unit",
    "reference_gene",
    "callable",
    "callability_reason",
    "target_chromosome",
    "target_interval_start_1based",
    "target_interval_end_1based",
    "qualifying_genome_hit_count",
    "qualifying_local_hit_count",
    "best_hit_subject",
    "best_hit_percent_identity",
    "best_hit_alignment_length",
    "best_hit_evalue",
    "best_hit_bitscore",
    "primary_state",
    "positive_loss",
    "historical_reproduction_state",
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
        raise InputError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def strict_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise InputError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise InputError(f"{path}: JSON root is not an object")
    return value


def resolve(root: Path, value: str, *, file: bool = False) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise InputError(f"unsafe data-root-relative path: {value!r}")
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise InputError(f"path escapes data root: {value!r}")
    if file and (not path.is_file() or path.stat().st_size == 0):
        raise InputError(f"missing input file: {path}")
    return path


def read_rows(path: Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != columns:
            raise InputError(f"{path.name}: columns differ from exact schema")
        rows = list(reader)
    if not rows:
        raise InputError(f"{path.name}: no data rows")
    return rows


def read_manifest(path: Path) -> list[dict[str, str]]:
    rows = read_rows(path, MANIFEST_COLUMNS)
    units = [row["unit"] for row in rows]
    if any(not value for row in rows for value in row.values()) or len(units) != len(set(units)):
        raise InputError("manifest contains empty values or duplicate units")
    return rows


def read_checksum_table(path: Path) -> dict[str, dict[str, object]]:
    rows = read_rows(path, ("file", "bytes", "sha256"))
    result: dict[str, dict[str, object]] = {}
    for row in rows:
        name = row["file"]
        if not name or Path(name).name != name or name in result:
            raise InputError(f"{path.name}: unsafe or duplicate filename {name!r}")
        try:
            byte_count = int(row["bytes"])
        except ValueError as error:
            raise InputError(f"{path.name}: non-integer byte count") from error
        if byte_count < 0 or len(row["sha256"]) != 64:
            raise InputError(f"{path.name}: invalid checksum row for {name!r}")
        result[name] = {"basename": name, "bytes": byte_count, "sha256": row["sha256"]}
    return result


def require_checksum(root: Path, table: dict[str, dict[str, object]], name: str) -> Path:
    path = root / name
    observed = binding(path, allow_empty=True)
    if table.get(name) != observed:
        raise InputError(f"{root.name}: checksum binding mismatch for {name}")
    return path


def read_label_map(path: Path) -> dict[str, str]:
    rows = read_rows(path, LABEL_COLUMNS)
    source = [row["query_chromosome"] for row in rows]
    final = [row["final_chromosome"] for row in rows]
    if (
        len(rows) != 29
        or len(set(source)) != 29
        or set(final) != {f"Chr{index:02d}" for index in range(1, 30)}
    ):
        raise InputError(f"{path.name}: label map is not a Chr01-Chr29 bijection")
    if any(
        row["coordinate_reference"] != "act_chinensis_hongyang_v4_hy4a"
        or row["assignment_method"] != "global_one_to_one_maximum_nucleotide_similarity"
        for row in rows
    ):
        raise InputError(f"{path.name}: label map does not use the frozen HY4A policy")
    return dict(zip(source, final, strict=True))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    staging: Path | None = None
    try:
        root = args.data_root.resolve()
        rows = read_manifest(args.manifest)
        output = args.output_dir.resolve()
        if output.exists():
            raise InputError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))

        calls: list[dict[str, object]] = []
        coordinates: list[dict[str, object]] = []
        assemblies: list[dict[str, object]] = []
        unit_audit: list[dict[str, object]] = []
        seen_calls: set[tuple[str, str]] = set()

        for row in rows:
            unit = row["unit"]
            search_dir = resolve(root, row["search_dir"])
            relabel_dir = resolve(root, row["relabel_dir"])
            label_path = resolve(root, row["label_map"], file=True)
            label_map = read_label_map(label_path)

            search_manifest_path = search_dir / "run_manifest.json"
            search_manifest = strict_json(search_manifest_path)
            if search_manifest.get("status") != "PASS" or search_manifest.get("unit") != unit:
                raise InputError(f"{unit}: translated-search manifest is not exact PASS")
            search_checksums = read_checksum_table(search_dir / "checksums.tsv")
            state_path = require_checksum(search_dir, search_checksums, "loss_states.tsv")
            require_checksum(search_dir, search_checksums, "run_manifest.json")
            state_rows = read_rows(state_path, STATE_COLUMNS)
            metrics = search_manifest.get("metrics")
            if not isinstance(metrics, dict) or metrics.get("candidate_rows") != len(state_rows):
                raise InputError(f"{unit}: translated-search row-count closure failed")

            relabel_manifest_path = relabel_dir / "run_manifest.json"
            relabel_manifest = strict_json(relabel_manifest_path)
            if relabel_manifest.get("status") != "PASS" or relabel_manifest.get("unit") != unit:
                raise InputError(f"{unit}: chromosome-relabel manifest is not exact PASS")
            relabel_inputs = relabel_manifest.get("inputs")
            if not isinstance(relabel_inputs, dict) or relabel_inputs.get("label_map") != binding(label_path):
                raise InputError(f"{unit}: relabel output is not bound to the selected label map")
            relabel_checksums = read_checksum_table(relabel_dir / "checksums.tsv")
            genome_name = f"{unit}.genome.fa.gz"
            gff_name = f"{unit}.primary.gff3"
            genome = require_checksum(relabel_dir, relabel_checksums, genome_name)
            gff = require_checksum(relabel_dir, relabel_checksums, gff_name)
            require_checksum(relabel_dir, relabel_checksums, "run_manifest.json")

            unit_seen: set[str] = set()
            positive_count = 0
            for state in state_rows:
                gene = state["reference_gene"]
                if state["unit"] != unit or not gene or gene in unit_seen:
                    raise InputError(f"{unit}: translated-search unit/gene closure failed")
                unit_seen.add(gene)
                if state["primary_state"] != "positive_deleted":
                    continue
                if state["callable"] != "true" or state["positive_loss"] != "true":
                    raise InputError(f"{unit}/{gene}: positive deletion is not callable/positive")
                source_chromosome = state["target_chromosome"]
                if source_chromosome not in label_map:
                    raise InputError(f"{unit}/{gene}: target chromosome is absent from label map")
                try:
                    start = int(state["target_interval_start_1based"])
                    end = int(state["target_interval_end_1based"])
                except ValueError as error:
                    raise InputError(f"{unit}/{gene}: invalid expected-locus interval") from error
                start, end = min(start, end), max(start, end)
                if start < 1:
                    raise InputError(f"{unit}/{gene}: expected-locus interval begins before 1")
                key = (unit, gene)
                if key in seen_calls:
                    raise InputError(f"duplicate positive call {unit}/{gene}")
                seen_calls.add(key)
                positive_count += 1
                calls.append(
                    {
                        "assembly_unit_id": unit,
                        "reference_gene_id": gene,
                        "classification": "positive_deleted",
                        "callable": "true",
                    }
                )
                coordinates.append(
                    {
                        "assembly_unit_id": unit,
                        "reference_gene_id": gene,
                        "classification": "positive_deleted",
                        "chromosome": label_map[source_chromosome],
                        "expected_locus_start_1based": start,
                        "expected_locus_end_1based": end,
                        "coordinate_semantics": "midpoint_of_bilateral_synorth_bounded_search_interval",
                    }
                )
            if metrics.get("positive_deleted") != positive_count:
                raise InputError(f"{unit}: positive-deletion count closure failed")

            assemblies.append(
                {
                    "assembly_unit_id": unit,
                    "biological_species": row["biological_species"],
                    "haplotype_or_subgenome": row["haplotype_or_subgenome"],
                    "assembly_scope": row["assembly_scope"],
                    "genome": os.path.relpath(genome, staging),
                    "gff": os.path.relpath(gff, staging),
                    "genome_local_sha256": sha256(genome),
                    "gff_local_sha256": sha256(gff),
                }
            )
            unit_audit.append(
                {
                    "unit": unit,
                    "translated_search_manifest": binding(search_manifest_path),
                    "loss_states": binding(state_path),
                    "label_map": binding(label_path),
                    "relabel_manifest": binding(relabel_manifest_path),
                    "genome": binding(genome),
                    "gff": binding(gff),
                    "state_rows": len(state_rows),
                    "positive_deleted": positive_count,
                }
            )

        calls.sort(key=lambda item: (str(item["assembly_unit_id"]), str(item["reference_gene_id"])))
        coordinates.sort(
            key=lambda item: (str(item["assembly_unit_id"]), str(item["reference_gene_id"]))
        )
        assemblies.sort(key=lambda item: str(item["assembly_unit_id"]))
        calls_path = staging / "positive_deleted_calls.tsv"
        coordinates_path = staging / "expected_deleted_locus_coordinates.tsv"
        assemblies_path = staging / "assembly_manifest.tsv"
        write_tsv(
            calls_path,
            calls,
            ("assembly_unit_id", "reference_gene_id", "classification", "callable"),
        )
        write_tsv(
            coordinates_path,
            coordinates,
            (
                "assembly_unit_id",
                "reference_gene_id",
                "classification",
                "chromosome",
                "expected_locus_start_1based",
                "expected_locus_end_1based",
                "coordinate_semantics",
            ),
        )
        write_tsv(
            assemblies_path,
            assemblies,
            (
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "assembly_scope",
                "genome",
                "gff",
                "genome_local_sha256",
                "gff_local_sha256",
            ),
        )
        report = {
            "schema_version": 1,
            "workflow": "callable_positive_deleted_expected_locus_spatial_inputs",
            "status": "PASS",
            "coordinate_semantics": (
                "expected locus interval bounded by bilateral accepted SynOrths anchors; "
                "the downstream position is its midpoint, not an observed remnant fragment"
            ),
            "chromosome_labels": "HY4A-compatible Chr01-Chr29; publisher direction preserved",
            "unit_count": len(rows),
            "positive_deleted_count": len(calls),
            "source_manifest": binding(args.manifest),
            "units": unit_audit,
            "outputs": {
                "positive_calls": binding(calls_path),
                "expected_locus_coordinates": binding(coordinates_path),
                "assembly_manifest": binding(assemblies_path),
            },
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        checksums = staging / "checksums.tsv"
        with checksums.open("w", encoding="utf-8", newline="") as handle:
            handle.write("file\tbytes\tsha256\n")
            for path in sorted(staging.iterdir()):
                if path.is_file() and path.name != checksums.name:
                    record = binding(path, allow_empty=True)
                    handle.write(f"{path.name}\t{record['bytes']}\t{record['sha256']}\n")
        os.rename(staging, output)
        staging = None
        print(f"PASS\t{len(rows)} units\t{len(calls)} positive deleted loci")
        return 0
    except (OSError, UnicodeError, csv.Error, json.JSONDecodeError, InputError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    finally:
        if staging is not None:
            import shutil

            shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
