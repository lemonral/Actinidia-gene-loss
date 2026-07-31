#!/usr/bin/env python3
"""Convert exact-bound legacy Actinidia evidence into complete call matrices.

The primary matrix is deliberately conservative.  An exact SynOrths anchor is
retained; a historical genome search with no qualifying translated hit is a
positive deletion; a historical ``decayed`` hit is uncertain because the old
output does not prove a disruptive mutation.  The manuscript-era
pseudogenized/deleted interpretation is written separately and never silently
mixed into the revised primary analysis.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from pathlib import Path


class LegacyMatrixError(RuntimeError):
    pass


MANIFEST_COLUMNS = (
    "assembly_unit_id",
    "legacy_sample",
    "synorth_audit",
    "synorth_pairs",
    "decayed_genes",
    "deleted_genes",
)
PRIMARY_COLUMNS = (
    "reference_gene_id",
    "assembly_unit_id",
    "classification",
    "callable",
    "evidence_source",
    "primary_search_state",
)
HISTORICAL_COLUMNS = (
    "reference_gene_id",
    "assembly_unit_id",
    "classification",
    "positive_loss",
    "evidence_source",
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
        raise LegacyMatrixError(f"missing or empty file: {source}")
    return {
        "basename": source.name,
        "bytes": source.stat().st_size,
        "sha256": sha256(source),
    }


def resolve(
    root: Path, value: str, *, file: bool = True, allow_empty: bool = False
) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise LegacyMatrixError(f"unsafe data-root-relative path: {value!r}")
    path = (root / relative).absolute()
    if not path.is_relative_to(root):
        raise LegacyMatrixError(f"path escapes data root: {value!r}")
    if file and (not path.is_file() or (not allow_empty and path.stat().st_size == 0)):
        raise LegacyMatrixError(f"missing input file: {path}")
    return path


def read_fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                identifier = raw[1:].strip().split()[0]
                if not identifier:
                    raise LegacyMatrixError(f"{path.name}:{line_number}: empty FASTA ID")
                identifiers.append(identifier)
    if not identifiers or len(identifiers) != len(set(identifiers)):
        raise LegacyMatrixError(f"{path.name}: empty or duplicate FASTA IDs")
    return identifiers


def read_id_list(path: Path) -> list[str]:
    values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(values) != len(set(values)):
        raise LegacyMatrixError(f"{path.name}: duplicate identifier rows")
    if any(len(value.split()) != 1 for value in values):
        raise LegacyMatrixError(f"{path.name}: identifiers must be one whitespace-free field")
    return values


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MANIFEST_COLUMNS:
            raise LegacyMatrixError("legacy manifest columns differ from exact schema")
        rows = list(reader)
    units = [row["assembly_unit_id"] for row in rows]
    samples = [row["legacy_sample"] for row in rows]
    if (
        not rows
        or any(not value for row in rows for value in row.values())
        or len(units) != len(set(units))
        or len(samples) != len(set(samples))
    ):
        raise LegacyMatrixError("legacy manifest is empty or has empty/duplicate unit identities")
    return rows


def read_audit(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise LegacyMatrixError(f"cannot read audit {path}: {error}") from error
    if not isinstance(value, dict):
        raise LegacyMatrixError(f"{path.name}: audit root is not an object")
    return value


def read_anchor_ids(path: Path, column_1_based: int, has_header: bool) -> set[str]:
    result: set[str] = set()
    seen_pairs: set[tuple[str, ...]] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if line_number == 1 and has_header:
                continue
            if len(fields) < column_1_based or not fields[column_1_based - 1]:
                raise LegacyMatrixError(f"{path.name}:{line_number}: invalid SynOrths row")
            pair = tuple(fields)
            if pair in seen_pairs:
                raise LegacyMatrixError(f"{path.name}:{line_number}: duplicate SynOrths pair row")
            seen_pairs.add(pair)
            result.add(fields[column_1_based - 1])
    if not result:
        raise LegacyMatrixError(f"{path.name}: no reference anchors")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--reference-cds", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    staging: Path | None = None
    try:
        root = args.data_root.resolve()
        rows = read_manifest(args.manifest)
        reference_cds = resolve(root, args.reference_cds)
        reference_ids = read_fasta_ids(reference_cds)
        reference_set = set(reference_ids)
        output = args.output_dir.absolute()
        if output.exists():
            raise LegacyMatrixError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        primary_path = staging / "complete_unit_loss_matrix.tsv"
        historical_path = staging / "historical_reproduction_loss_matrix.tsv"
        audits: list[dict[str, object]] = []
        with (
            primary_path.open("w", encoding="utf-8", newline="") as primary_handle,
            historical_path.open("w", encoding="utf-8", newline="") as historical_handle,
        ):
            primary_writer = csv.DictWriter(
                primary_handle, fieldnames=PRIMARY_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            historical_writer = csv.DictWriter(
                historical_handle, fieldnames=HISTORICAL_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            primary_writer.writeheader()
            historical_writer.writeheader()
            for row in rows:
                unit = row["assembly_unit_id"]
                sample = row["legacy_sample"]
                audit_path = resolve(root, row["synorth_audit"])
                pairs_path = resolve(root, row["synorth_pairs"])
                decayed_path = resolve(root, row["decayed_genes"], allow_empty=True)
                deleted_path = resolve(root, row["deleted_genes"], allow_empty=True)
                audit = read_audit(audit_path)
                inputs = audit.get("inputs")
                metrics = audit.get("metrics")
                if (
                    audit.get("schema_version") != 2
                    or audit.get("sample") != sample
                    or not isinstance(inputs, dict)
                    or not isinstance(metrics, dict)
                ):
                    raise LegacyMatrixError(f"{unit}: incompatible SynOrths audit")
                audited_pairs = Path(str(inputs.get("pairs", ""))).resolve()
                if audited_pairs != pairs_path.resolve():
                    raise LegacyMatrixError(f"{unit}: SynOrths audit/pair realpath mismatch")
                try:
                    reference_column = int(inputs["reference_column_1_based"])
                except (KeyError, TypeError, ValueError) as error:
                    raise LegacyMatrixError(f"{unit}: invalid audited reference column") from error
                has_header = inputs.get("pairs_has_header") is True
                anchors = read_anchor_ids(pairs_path, reference_column, has_header)
                if (
                    int(metrics.get("duplicate_pair_rows", -1)) != 0
                    or int(metrics.get("unique_reference_anchors", -1)) != len(anchors)
                    or int(metrics.get("reference_anchor_ids_absent_from_fasta_count", -1)) != 0
                    or int(metrics.get("reference_anchor_ids_absent_from_coordinates_count", -1)) != 0
                ):
                    raise LegacyMatrixError(f"{unit}: SynOrths audit metric closure failed")
                outside_anchors = anchors.difference(reference_set)
                if outside_anchors:
                    raise LegacyMatrixError(f"{unit}: {len(outside_anchors)} anchors outside reference CDS")
                decayed = set(read_id_list(decayed_path))
                deleted = set(read_id_list(deleted_path))
                if decayed.intersection(deleted):
                    raise LegacyMatrixError(f"{unit}: decayed/deleted lists overlap")
                positive = decayed.union(deleted)
                if positive.difference(reference_set):
                    raise LegacyMatrixError(f"{unit}: historical calls outside reference CDS")
                if positive.intersection(anchors):
                    raise LegacyMatrixError(f"{unit}: historical calls overlap SynOrths anchors")
                primary_counts: Counter[str] = Counter()
                historical_counts: Counter[str] = Counter()
                for gene in reference_ids:
                    if gene in anchors:
                        classification, callable_value, source, state = (
                            "retained", "true", "exact_bound_legacy_synorth_anchor", "not_searched_anchor"
                        )
                        historical_class, historical_positive = "retained", "false"
                    elif gene in deleted:
                        classification, callable_value, source, state = (
                            "deleted", "true", "historical_genomewide_translated_search", "historical_no_qualifying_hit"
                        )
                        historical_class, historical_positive = "deleted", "true"
                    elif gene in decayed:
                        classification, callable_value, source, state = (
                            "uncertain", "true", "historical_genomewide_translated_search", "historical_sequence_detected_no_mutation_proof"
                        )
                        historical_class, historical_positive = "pseudogenized", "true"
                    else:
                        classification, callable_value, source, state = (
                            "uncertain", "false", "outside_exact_bound_legacy_evidence", "not_called_loss"
                        )
                        historical_class, historical_positive = "not_called_loss", "false"
                    primary_counts[classification] += 1
                    historical_counts[historical_class] += 1
                    primary_writer.writerow(
                        {
                            "reference_gene_id": gene,
                            "assembly_unit_id": unit,
                            "classification": classification,
                            "callable": callable_value,
                            "evidence_source": source,
                            "primary_search_state": state,
                        }
                    )
                    historical_writer.writerow(
                        {
                            "reference_gene_id": gene,
                            "assembly_unit_id": unit,
                            "classification": historical_class,
                            "positive_loss": historical_positive,
                            "evidence_source": "manuscript_era_decayed_deleted_reproduction",
                        }
                    )
                if sum(primary_counts.values()) != len(reference_ids):
                    raise LegacyMatrixError(f"{unit}: primary count closure failed")
                audits.append(
                    {
                        "assembly_unit_id": unit,
                        "legacy_sample": sample,
                        "reference_gene_count": len(reference_ids),
                        "anchor_count": len(anchors),
                        "historical_decayed_count": len(decayed),
                        "historical_deleted_count": len(deleted),
                        "primary_counts": dict(sorted(primary_counts.items())),
                        "historical_counts": dict(sorted(historical_counts.items())),
                        "synorth_audit": binding(audit_path),
                        "synorth_pairs": binding(pairs_path),
                        "decayed_genes": binding(decayed_path, allow_empty=True),
                        "deleted_genes": binding(deleted_path, allow_empty=True),
                    }
                )
        expected_rows = len(reference_ids) * len(rows)
        report = {
            "schema_version": 1,
            "workflow": "complete_exact_bound_legacy_primary_loss_matrix",
            "status": "PASS",
            "definitions": {
                "primary_retained": "exact-bound legacy SynOrths reference anchor",
                "primary_deleted": "historical genome-wide translated search reported no qualifying hit",
                "primary_uncertain_sequence": "historical decayed hit lacks independently bound disruptive-mutation evidence",
                "primary_uncertain_non_callable": "gene lies outside exact-bound anchor and historical call evidence",
                "historical_reproduction": "preserves manuscript-era decayed=pseudogenized and deleted labels separately",
            },
            "reference_cds": binding(reference_cds),
            "reference_gene_count": len(reference_ids),
            "assembly_unit_count": len(rows),
            "matrix_rows": expected_rows,
            "source_manifest": binding(args.manifest),
            "units": audits,
            "outputs": {
                "complete_unit_loss_matrix": binding(primary_path),
                "historical_reproduction_loss_matrix": binding(historical_path),
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
                    item = binding(path, allow_empty=True)
                    handle.write(f"{path.name}\t{item['bytes']}\t{item['sha256']}\n")
        os.rename(staging, output)
        staging = None
        print(json.dumps({"status": "PASS", "units": len(rows), "rows": expected_rows}, sort_keys=True))
        return 0
    except (LegacyMatrixError, OSError, UnicodeError, csv.Error, ValueError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    finally:
        if staging is not None and staging.exists():
            shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
