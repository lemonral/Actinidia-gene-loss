#!/usr/bin/env python3
"""Build the complete callable-aware reference-gene matrix for new Actinidia units.

For each reference CDS gene and assembly unit:

* an exact SynOrths anchor is ``retained`` and callable;
* a searched callable locus with no qualifying local hit is ``deleted``;
* a searched callable locus containing genomic sequence is ``uncertain``;
* every other unanchored locus is ``uncertain`` and non-callable.

This prevents absence from a positive list from being interpreted as retained.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from pathlib import Path


class MatrixError(RuntimeError):
    pass


MANIFEST_COLUMNS = ("unit", "synorth_pairs", "candidate_dir", "search_dir")
STATE_COLUMNS = (
    "unit", "reference_gene", "callable", "callability_reason", "target_chromosome",
    "target_interval_start_1based", "target_interval_end_1based",
    "qualifying_genome_hit_count", "qualifying_local_hit_count", "best_hit_subject",
    "best_hit_percent_identity", "best_hit_alignment_length", "best_hit_evalue",
    "best_hit_bitscore", "primary_state", "positive_loss", "historical_reproduction_state",
)
OUTPUT_COLUMNS = (
    "reference_gene_id", "assembly_unit_id", "classification", "callable",
    "evidence_source", "primary_search_state",
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
        raise MatrixError(f"missing or empty file: {source}")
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def resolve(root: Path, value: str, *, file: bool = False) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise MatrixError(f"unsafe data-root-relative path: {value!r}")
    # Check lexical containment without dereferencing the final asset.  The
    # frozen data store intentionally exposes some immutable legacy inputs as
    # symlinks (for example the C. scandens reference CDS).  Rejecting those
    # merely because their registered target lives in the legacy data tree
    # would make the production manifest unusable.  ``absolute()`` removes
    # harmless ``.`` components but preserves the symlink for this check;
    # ``binding()`` still dereferences it and records the exact target bytes.
    path = (root / relative).absolute()
    if not path.is_relative_to(root):
        raise MatrixError(f"path escapes data root: {value!r}")
    if file and (not path.is_file() or path.stat().st_size == 0):
        raise MatrixError(f"missing input file: {path}")
    return path


def read_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise MatrixError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise MatrixError(f"{path}: JSON root is not an object")
    return value


def read_rows(path: Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != columns:
            raise MatrixError(f"{path.name}: columns differ from exact schema")
        rows = list(reader)
    if not rows:
        raise MatrixError(f"{path.name}: no data rows")
    return rows


def read_manifest(path: Path) -> list[dict[str, str]]:
    rows = read_rows(path, MANIFEST_COLUMNS)
    units = [row["unit"] for row in rows]
    if any(not value for row in rows for value in row.values()) or len(units) != len(set(units)):
        raise MatrixError("manifest has empty values or duplicate units")
    return rows


def read_fasta_ids(path: Path) -> set[str]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    identifiers: set[str] = set()
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            identifier = raw[1:].strip().split()[0]
            if not identifier or identifier in identifiers:
                raise MatrixError(f"{path.name}:{line_number}: duplicate/empty FASTA ID")
            identifiers.add(identifier)
    if not identifiers:
        raise MatrixError(f"{path.name}: no FASTA IDs")
    return identifiers


def read_synorth_reference_ids(path: Path) -> set[str]:
    identifiers: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if len(fields) < 8 or not fields[4]:
                raise MatrixError(f"{path.name}:{line_number}: invalid SynOrths row")
            identifiers.add(fields[4])
    if not identifiers:
        raise MatrixError(f"{path.name}: no reference anchors")
    return identifiers


def checksum_table(path: Path) -> dict[str, dict[str, object]]:
    rows = read_rows(path, ("file", "bytes", "sha256"))
    result: dict[str, dict[str, object]] = {}
    for row in rows:
        name = row["file"]
        if not name or Path(name).name != name or name in result:
            raise MatrixError(f"{path.name}: unsafe/duplicate checksum filename")
        try:
            size = int(row["bytes"])
        except ValueError as error:
            raise MatrixError(f"{path.name}: invalid byte count") from error
        if size < 0 or len(row["sha256"]) != 64:
            raise MatrixError(f"{path.name}: invalid checksum row")
        result[name] = {"basename": name, "bytes": size, "sha256": row["sha256"]}
    return result


def require_checksum(root: Path, table: dict[str, dict[str, object]], name: str) -> Path:
    path = root / name
    if table.get(name) != binding(path, allow_empty=True):
        raise MatrixError(f"{root.name}: checksum mismatch for {name}")
    return path


def candidate_ids(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "reference_gene" not in (reader.fieldnames or []):
            raise MatrixError(f"{path.name}: missing reference_gene")
        identifiers = [row["reference_gene"] for row in reader]
    if not identifiers or any(not item for item in identifiers) or len(identifiers) != len(set(identifiers)):
        raise MatrixError(f"{path.name}: empty/duplicate candidate IDs")
    return set(identifiers)


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
        reference_cds = resolve(root, args.reference_cds, file=True)
        reference_ids = read_fasta_ids(reference_cds)
        output = args.output_dir.resolve()
        if output.exists():
            raise MatrixError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        matrix_path = staging / "complete_unit_loss_matrix.tsv"
        audits: list[dict[str, object]] = []
        total_rows = 0
        with matrix_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                unit = row["unit"]
                synorth = resolve(root, row["synorth_pairs"], file=True)
                candidates_root = resolve(root, row["candidate_dir"])
                search_root = resolve(root, row["search_dir"])
                candidate_manifest = read_json(candidates_root / "run_manifest.json")
                candidate_inputs = candidate_manifest.get("inputs")
                candidate_outputs = candidate_manifest.get("outputs")
                if (
                    candidate_manifest.get("status") != "PASS"
                    or candidate_manifest.get("unit") != unit
                    or not isinstance(candidate_inputs, dict)
                    or not isinstance(candidate_outputs, dict)
                    or candidate_inputs.get("reference_cds") != binding(reference_cds)
                    or candidate_inputs.get("synorth_pairs") != binding(synorth)
                ):
                    raise MatrixError(f"{unit}: candidate manifest input binding failed")
                candidate_table = candidates_root / "candidates.tsv"
                if candidate_outputs.get("candidates") != binding(candidate_table):
                    raise MatrixError(f"{unit}: candidate-table output binding failed")
                candidates = candidate_ids(candidate_table)

                search_manifest = read_json(search_root / "run_manifest.json")
                if search_manifest.get("status") != "PASS" or search_manifest.get("unit") != unit:
                    raise MatrixError(f"{unit}: search manifest is not exact PASS")
                table = checksum_table(search_root / "checksums.tsv")
                states_path = require_checksum(search_root, table, "loss_states.tsv")
                require_checksum(search_root, table, "run_manifest.json")
                states = read_rows(states_path, STATE_COLUMNS)
                state_by_gene = {item["reference_gene"]: item for item in states}
                if len(state_by_gene) != len(states) or set(state_by_gene) != candidates:
                    raise MatrixError(f"{unit}: candidate/search gene-ID closure failed")
                search_metrics = search_manifest.get("metrics")
                if not isinstance(search_metrics, dict) or search_metrics.get("candidate_rows") != len(states):
                    raise MatrixError(f"{unit}: search row-count closure failed")

                anchored_all = read_synorth_reference_ids(synorth)
                overlap = anchored_all.intersection(candidates)
                if overlap:
                    raise MatrixError(f"{unit}: {len(overlap)} genes are both anchors and candidates")
                anchored = anchored_all.intersection(reference_ids)
                counts: Counter[str] = Counter()
                for gene in sorted(reference_ids):
                    if gene in anchored:
                        classification, callable_value, source, primary = (
                            "retained", "true", "exact_synorth_anchor", "not_searched_anchor",
                        )
                    elif gene in state_by_gene:
                        state = state_by_gene[gene]
                        primary = state["primary_state"]
                        if primary == "positive_deleted":
                            if state["callable"] != "true" or state["positive_loss"] != "true":
                                raise MatrixError(f"{unit}/{gene}: inconsistent positive state")
                            classification, callable_value = "deleted", "true"
                        elif primary == "uncertain_local_genomic_sequence_detected":
                            if state["callable"] != "true" or state["positive_loss"] != "false":
                                raise MatrixError(f"{unit}/{gene}: inconsistent local-sequence state")
                            classification, callable_value = "uncertain", "true"
                        elif primary == "uncertain_non_callable":
                            if state["callable"] != "false" or state["positive_loss"] != "false":
                                raise MatrixError(f"{unit}/{gene}: inconsistent non-callable state")
                            classification, callable_value = "uncertain", "false"
                        else:
                            raise MatrixError(f"{unit}/{gene}: unsupported primary state {primary!r}")
                        source = "callable_aware_translated_search"
                    else:
                        classification, callable_value, source, primary = (
                            "uncertain", "false", "outside_synorth_local_candidate_scope", "not_searched_non_callable",
                        )
                    counts[f"classification_{classification}"] += 1
                    counts[f"callable_{callable_value}"] += 1
                    writer.writerow(
                        {
                            "reference_gene_id": gene,
                            "assembly_unit_id": unit,
                            "classification": classification,
                            "callable": callable_value,
                            "evidence_source": source,
                            "primary_search_state": primary,
                        }
                    )
                    total_rows += 1
                if sum(counts[key] for key in counts if key.startswith("classification_")) != len(reference_ids):
                    raise MatrixError(f"{unit}: classification count closure failed")
                audits.append(
                    {
                        "unit": unit,
                        "reference_gene_count": len(reference_ids),
                        "synorth_reference_ids_all": len(anchored_all),
                        "synorth_reference_ids_in_callable_universe": len(anchored),
                        "candidate_rows": len(candidates),
                        "candidate_rows_outside_callable_universe": len(candidates.difference(reference_ids)),
                        "counts": dict(sorted(counts.items())),
                        "synorth_pairs": binding(synorth),
                        "candidate_manifest": binding(candidates_root / "run_manifest.json"),
                        "search_manifest": binding(search_root / "run_manifest.json"),
                        "loss_states": binding(states_path),
                    }
                )
        expected_rows = len(reference_ids) * len(rows)
        if total_rows != expected_rows:
            raise MatrixError(f"matrix row closure failed: {total_rows} != {expected_rows}")
        report = {
            "schema_version": 1,
            "workflow": "complete_callable_aware_new_unit_loss_matrix",
            "status": "PASS",
            "definitions": {
                "retained": "reference gene has at least one exact-bound SynOrths target anchor",
                "deleted": "callable bilateral local locus has no qualifying local translated hit",
                "uncertain_callable": "callable local locus contains translated genomic sequence but lacks disruptive-mutation evidence",
                "uncertain_non_callable": "unanchored reference gene lacks an accepted local callable interval",
                "absence_from_positive_list": "never treated as retained",
            },
            "reference_cds": binding(reference_cds),
            "reference_gene_count": len(reference_ids),
            "assembly_unit_count": len(rows),
            "matrix_rows": total_rows,
            "expected_matrix_rows": expected_rows,
            "source_manifest": binding(args.manifest),
            "units": audits,
            "outputs": {"complete_unit_loss_matrix": binding(matrix_path)},
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
        print(f"PASS\t{len(rows)} units\t{len(reference_ids)} genes\t{total_rows} rows")
        return 0
    except (OSError, UnicodeError, csv.Error, json.JSONDecodeError, MatrixError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    finally:
        if staging is not None:
            shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
