#!/usr/bin/env python3
"""Validate one completed NLR batch and prepare non-shared loss inputs.

The reference NLR universe is the set of unique reference-CDS sequence IDs in
the NLR-Annotator calls.  Genes positive-complete in every primary lineage are
excluded before downstream NLR loss rates are calculated.  For each assembly
unit the denominator is the catalog of remaining reference-NLR genes whose
complete primary loss-matrix row is resolved as retained, deleted, or strictly
supported pseudogenized.  The numerator is deleted plus strictly supported
pseudogenized.  Uncertain comparisons are excluded from both numerator and
denominator even when an interval itself was callable.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence


SHA256 = re.compile(r"^[0-9a-f]{64}$")
SELECTED_FIELDS = (
    "sample_id", "species", "ploidy", "analysis_role", "input_scope",
    "relative_fasta", "expected_fasta_records",
)
METADATA_FIELDS = (
    "assembly_unit_id", "biological_species", "haplotype_or_subgenome",
    "assembly_scope", "include", "analysis_cohort",
)
LOSS_FIELDS = (
    "reference_gene_id", "assembly_unit_id", "classification", "callable",
    "evidence_source", "primary_search_state",
)
SAMPLE_FILES = {
    "nlr_calls.txt", "nlr_loci.gff", "stdout.log", "stderr.log",
    "run_metadata.tsv", "output_checksums.tsv",
}
CHECKSUM_SAMPLE_FILES = {"nlr_calls.txt", "nlr_loci.gff", "stdout.log", "stderr.log"}
ROOT_FILES = {"selected_inputs.tsv", "batch_metadata.tsv"}
OPTIONAL_ROOT_FILES = {"resume_history.tsv"}
ALLOWED_CLASSES = {"deleted", "pseudogenized", "retained", "uncertain"}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nlr-root", type=Path, required=True)
    parser.add_argument("--input-bundle", type=Path, required=True)
    parser.add_argument("--unit-metadata", type=Path, required=True)
    parser.add_argument("--loss-matrix", type=Path, required=True)
    parser.add_argument("--shared-positive-genes", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-units", type=int, default=23)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    parser.add_argument("--expected-shared-positive", type=int, default=68)
    parser.add_argument("--expected-worker-threads", type=int, default=8)
    parser.add_argument("--expected-jar-sha256")
    parser.add_argument("--expected-motifs-sha256")
    parser.add_argument("--expected-store-sha256")
    return parser.parse_args(argv)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def regular_file(path: Path, *, allow_empty: bool = False) -> None:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected a regular non-symlink file: {path}")
    if not allow_empty and path.stat().st_size == 0:
        raise ValueError(f"Expected a non-empty file: {path}")


def read_tsv(
    path: Path,
    *,
    exact_fields: Sequence[str] | None = None,
    required_fields: Iterable[str] = (),
    allow_empty_file: bool = False,
) -> tuple[list[str], list[dict[str, str]]]:
    regular_file(path, allow_empty=allow_empty_file)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if exact_fields is not None and tuple(fields) != tuple(exact_fields):
            raise ValueError(f"Unexpected columns in {path.name}: {fields}")
        missing = set(required_fields) - set(fields)
        if missing:
            raise ValueError(f"{path.name} lacks columns: {sorted(missing)}")
        rows: list[dict[str, str]] = []
        for line_number, raw in enumerate(reader, 2):
            if None in raw:
                raise ValueError(f"Extra tab-delimited fields at {path}:{line_number}")
            rows.append({field: (raw.get(field) or "").strip() for field in fields})
    return fields, rows


def read_key_values(path: Path) -> dict[str, str]:
    _, rows = read_tsv(path, exact_fields=("key", "value"))
    values: dict[str, str] = {}
    for line_number, row in enumerate(rows, 2):
        key = row["key"]
        if not key or key in values:
            raise ValueError(f"Empty or duplicate metadata key at {path}:{line_number}")
        values[key] = row["value"]
    return values


def canonical_int(value: str, *, context: str, minimum: int = 0) -> int:
    if not re.fullmatch(r"0|[1-9][0-9]*", value):
        raise ValueError(f"{context} is not a canonical non-negative integer: {value!r}")
    parsed = int(value)
    if parsed < minimum:
        raise ValueError(f"{context} must be at least {minimum}: {parsed}")
    return parsed


def expected_digest(value: str | None, *, label: str) -> str | None:
    if value is None:
        return None
    normalized = value.lower()
    if SHA256.fullmatch(normalized) is None:
        raise ValueError(f"Invalid expected SHA-256 for {label}")
    return normalized


def read_fasta_ids(path: Path) -> list[str]:
    regular_file(path)
    identifiers: list[str] = []
    seen: set[str] = set()
    has_sequence = False
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(b">"):
                if identifiers and not has_sequence:
                    raise ValueError(f"Empty FASTA record before {path}:{line_number}")
                try:
                    identifier = raw[1:].strip().split()[0].decode("utf-8")
                except (IndexError, UnicodeDecodeError) as exc:
                    raise ValueError(f"Invalid FASTA header at {path}:{line_number}") from exc
                if not identifier or identifier in seen:
                    raise ValueError(f"Empty or duplicate FASTA ID at {path}:{line_number}")
                identifiers.append(identifier)
                seen.add(identifier)
                has_sequence = False
            elif raw.strip():
                if not identifiers:
                    raise ValueError(f"Sequence precedes first FASTA header in {path}")
                has_sequence = True
    if not identifiers or not has_sequence:
        raise ValueError(f"Empty FASTA or final empty record: {path}")
    return identifiers


def read_bundle(bundle: Path, expected_samples: int) -> tuple[dict[str, dict[str, object]], dict[str, str]]:
    if bundle.is_symlink() or not bundle.is_dir():
        raise ValueError(f"Input bundle is not a regular directory: {bundle}")
    manifest_path = bundle / "run_manifest.json"
    checksums_path = bundle / "checksums.tsv"
    regular_file(manifest_path)
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid input-bundle manifest JSON: {exc}") from exc
    if manifest.get("status") != "PASS" or manifest.get("workflow") != "plain_fasta_nlr_input_bundle":
        raise ValueError("NLR input bundle is not a PASS plain-FASTA bundle")
    if manifest.get("reference_count") != 1 or manifest.get("target_count") != expected_samples - 1:
        raise ValueError("NLR input-bundle reference/target counts are not exact")
    _, checksum_rows = read_tsv(checksums_path, exact_fields=("file", "bytes", "sha256"))
    checksums: dict[str, dict[str, object]] = {}
    for line_number, row in enumerate(checksum_rows, 2):
        name = row["file"]
        if not name or Path(name).name != name or name in checksums:
            raise ValueError(f"Unsafe or duplicate input-bundle checksum row {line_number}")
        digest = row["sha256"].lower()
        if SHA256.fullmatch(digest) is None:
            raise ValueError(f"Invalid input-bundle checksum at line {line_number}")
        checksums[name] = {
            "bytes": canonical_int(row["bytes"], context=f"checksums:{line_number}:bytes"),
            "sha256": digest,
        }
    outputs = manifest.get("outputs")
    if not isinstance(outputs, list) or len(outputs) != expected_samples + 1:
        raise ValueError("Input-bundle manifest output inventory is not exact")
    by_sample: dict[str, dict[str, object]] = {}
    seen_outputs: set[str] = set()
    for item in outputs:
        if not isinstance(item, dict):
            raise ValueError("Invalid item in input-bundle output inventory")
        sample = item.get("sample_id")
        name = item.get("basename")
        if not isinstance(sample, str) or not isinstance(name, str) or name in seen_outputs:
            raise ValueError("Invalid or duplicate input-bundle output inventory item")
        seen_outputs.add(name)
        recorded = checksums.get(name)
        if recorded is None or item.get("sha256") != recorded["sha256"] or item.get("bytes") != recorded["bytes"]:
            raise ValueError(f"Input-bundle manifest/checksum mismatch for {name}")
        if sample != "bundle":
            if sample in by_sample:
                raise ValueError(f"Duplicate sample in input-bundle outputs: {sample}")
            records = item.get("fasta_records")
            if not isinstance(records, int) or records < 1:
                raise ValueError(f"Invalid FASTA record count for {sample}")
            by_sample[sample] = {
                "basename": name,
                "bytes": item["bytes"],
                "sha256": item["sha256"],
                "fasta_records": records,
            }
    if len(by_sample) != expected_samples:
        raise ValueError("Input-bundle sample inventory is not exact")
    return by_sample, {
        "run_manifest_sha256": sha256_file(manifest_path),
        "checksums_sha256": sha256_file(checksums_path),
    }


def parse_nlr_calls(path: Path) -> tuple[int, set[str], set[str]]:
    regular_file(path, allow_empty=True)
    rows = 0
    sequence_ids: set[str] = set()
    locus_ids: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, values in enumerate(reader, 1):
            if not values:
                continue
            if len(values) != 7 or not all(value.strip() for value in values[:3]):
                raise ValueError(f"Invalid seven-column NLR call at {path}:{line_number}")
            sequence, locus = values[0].strip(), values[1].strip()
            if locus in locus_ids:
                raise ValueError(f"Duplicate NLR locus ID {locus!r} in {path}")
            rows += 1
            sequence_ids.add(sequence)
            locus_ids.add(locus)
    return rows, sequence_ids, locus_ids


def validate_nlr_root(
    root: Path,
    bundle_samples: Mapping[str, Mapping[str, object]],
    *,
    expected_units: int,
    expected_workers: int,
    expected_tools: Mapping[str, str | None],
) -> tuple[list[dict[str, str]], dict[str, dict[str, object]], dict[str, str]]:
    if root.is_symlink() or not root.is_dir():
        raise ValueError(f"NLR root is not a regular directory: {root}")
    _, selected = read_tsv(root / "selected_inputs.tsv", exact_fields=SELECTED_FIELDS)
    if len(selected) != expected_units + 1:
        raise ValueError(f"Observed {len(selected)} selected NLR inputs; expected {expected_units + 1}")
    if len({row["sample_id"] for row in selected}) != len(selected):
        raise ValueError("NLR selected-input snapshot has duplicate sample IDs")
    references = [row for row in selected if row["analysis_role"] == "reference_callable"]
    targets = [row for row in selected if row["analysis_role"] == "target_repertoire"]
    if len(references) != 1 or len(targets) != expected_units:
        raise ValueError("NLR selected-input roles are not one reference plus the exact target cohort")
    if references[0]["input_scope"] != "reference_transcript_cds" or references[0]["sample_id"] != "clem_scandens_reference":
        raise ValueError("Unexpected reference NLR input identity or scope")
    if {row["input_scope"] for row in targets} != {"whole_genome"}:
        raise ValueError("Every target NLR input must use whole_genome scope")
    if set(bundle_samples) != {row["sample_id"] for row in selected}:
        raise ValueError("NLR selected inputs and plain-input bundle sample sets differ")

    batch = read_key_values(root / "batch_metadata.tsv")
    if batch.get("completion_status") != "complete":
        raise ValueError("NLR batch is not marked complete")
    if canonical_int(batch.get("selected_inputs", ""), context="batch:selected_inputs", minimum=1) != len(selected):
        raise ValueError("NLR batch selected-input count mismatch")
    tool_keys = {
        "jar": "nlr_annotator_jar_sha256",
        "motifs": "motifs_sha256",
        "store": "store_sha256",
    }
    tool_hashes: dict[str, str] = {}
    for short, key in tool_keys.items():
        observed = batch.get(key, "").lower()
        if SHA256.fullmatch(observed) is None:
            raise ValueError(f"NLR batch has invalid {key}")
        if expected_tools[short] is not None and observed != expected_tools[short]:
            raise ValueError(f"NLR batch {key} differs from the declared expected checksum")
        tool_hashes[short] = observed

    allowed = ROOT_FILES | OPTIONAL_ROOT_FILES | {row["sample_id"] for row in selected}
    observed_entries = {entry.name for entry in root.iterdir()}
    if not ROOT_FILES.issubset(observed_entries) or observed_entries - allowed:
        raise ValueError(
            "NLR root contents are not atomic/exact: "
            f"missing={sorted(ROOT_FILES - observed_entries)}, extra={sorted(observed_entries - allowed)}"
        )

    sample_results: dict[str, dict[str, object]] = {}
    for row in selected:
        sample = row["sample_id"]
        sample_dir = root / sample
        if sample_dir.is_symlink() or not sample_dir.is_dir():
            raise ValueError(f"Missing regular sample directory: {sample_dir}")
        contents = {entry.name for entry in sample_dir.iterdir()}
        if contents != SAMPLE_FILES:
            raise ValueError(f"Unexpected completed-sample contents for {sample}: {sorted(contents)}")
        for name in SAMPLE_FILES:
            regular_file(sample_dir / name, allow_empty=name in {"nlr_calls.txt", "stdout.log", "stderr.log"})
        metadata = read_key_values(sample_dir / "run_metadata.tsv")
        bundle = bundle_samples[sample]
        exact = {
            "sample_id": sample,
            "species": row["species"],
            "ploidy": row["ploidy"],
            "analysis_role": row["analysis_role"],
            "input_scope": row["input_scope"],
            "input_fasta_sha256": str(bundle["sha256"]),
            "input_fasta_records": str(bundle["fasta_records"]),
            "configured_nlr_worker_threads": str(expected_workers),
            "jvm_processor_cap": str(expected_workers),
            "completion_status": "complete",
            "nlr_annotator_jar_sha256": tool_hashes["jar"],
            "motifs_sha256": tool_hashes["motifs"],
            "store_sha256": tool_hashes["store"],
        }
        for key, expected in exact.items():
            if metadata.get(key) != expected:
                raise ValueError(f"NLR sample {sample} metadata mismatch for {key}")
        if Path(metadata.get("input_fasta", "")).name != row["relative_fasta"]:
            raise ValueError(f"NLR sample {sample} input FASTA basename mismatch")
        if row["relative_fasta"] != bundle["basename"]:
            raise ValueError(f"NLR sample {sample} selected/bundle FASTA mismatch")

        _, checksum_rows = read_tsv(
            sample_dir / "output_checksums.tsv", exact_fields=("path", "sha256")
        )
        recorded: dict[str, str] = {}
        for checksum_row in checksum_rows:
            name, digest = checksum_row["path"], checksum_row["sha256"].lower()
            if name in recorded or name not in CHECKSUM_SAMPLE_FILES or SHA256.fullmatch(digest) is None:
                raise ValueError(f"Invalid NLR sample checksum inventory for {sample}")
            recorded[name] = digest
        if set(recorded) != CHECKSUM_SAMPLE_FILES:
            raise ValueError(f"Incomplete NLR sample checksum inventory for {sample}")
        for name, digest in recorded.items():
            if sha256_file(sample_dir / name) != digest:
                raise ValueError(f"NLR sample output checksum mismatch for {sample}/{name}")

        call_rows, sequence_ids, locus_ids = parse_nlr_calls(sample_dir / "nlr_calls.txt")
        recorded_counts = {
            "nlr_output_rows": call_rows,
            "nlr_output_sequence_ids": len(sequence_ids),
            "nlr_output_locus_ids": len(locus_ids),
        }
        for key, observed in recorded_counts.items():
            if canonical_int(metadata.get(key, ""), context=f"{sample}:{key}") != observed:
                raise ValueError(f"NLR sample {sample} recorded output count mismatch for {key}")
        sample_results[sample] = {
            "nlr_rows": call_rows,
            "sequence_ids": sequence_ids,
            "locus_ids": locus_ids,
            "nlr_calls_sha256": recorded["nlr_calls.txt"],
        }
    return selected, sample_results, tool_hashes


def read_metadata(path: Path, expected_units: int) -> tuple[list[dict[str, str]], dict[str, dict[str, str]], str]:
    _, rows = read_tsv(path, exact_fields=METADATA_FIELDS)
    if len(rows) != expected_units:
        raise ValueError(f"Observed {len(rows)} unit metadata rows; expected {expected_units}")
    mapping: dict[str, dict[str, str]] = {}
    cohorts: set[str] = set()
    for line_number, row in enumerate(rows, 2):
        unit = row["assembly_unit_id"]
        if not unit or unit in mapping:
            raise ValueError(f"Empty or duplicate unit metadata ID at line {line_number}")
        if row["include"] != "true" or not row["biological_species"] or not row["assembly_scope"] or not row["analysis_cohort"]:
            raise ValueError(f"Incomplete or excluded primary unit metadata at line {line_number}")
        mapping[unit] = row
        cohorts.add(row["analysis_cohort"])
    if len(cohorts) != 1:
        raise ValueError("Primary NLR unit metadata must declare exactly one cohort")
    return rows, mapping, next(iter(cohorts))


def read_shared(path: Path, reference: set[str], expected: int) -> set[str]:
    _, rows = read_tsv(path, required_fields={"reference_gene_id"})
    shared = [row["reference_gene_id"] for row in rows]
    if len(shared) != expected or len(set(shared)) != expected:
        raise ValueError(f"Shared-positive table must contain exactly {expected} unique genes")
    if not set(shared).issubset(reference):
        raise ValueError("Shared-positive genes are not a subset of the reference-CDS universe")
    return set(shared)


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def run(args: argparse.Namespace) -> dict[str, object]:
    if min(
        args.expected_units,
        args.expected_reference_genes,
        args.expected_worker_threads,
    ) < 1 or args.expected_shared_positive < 0:
        raise ValueError("Expected counts and worker threads must be positive")
    if args.output_dir.exists() or args.output_dir.is_symlink():
        raise ValueError(f"Output directory already exists: {args.output_dir}")
    for path in (args.unit_metadata, args.loss_matrix, args.shared_positive_genes):
        regular_file(path)

    expected_tools = {
        "jar": expected_digest(args.expected_jar_sha256, label="JAR"),
        "motifs": expected_digest(args.expected_motifs_sha256, label="motifs"),
        "store": expected_digest(args.expected_store_sha256, label="store"),
    }
    bundle_samples, bundle_hashes = read_bundle(args.input_bundle, args.expected_units + 1)
    selected, nlr_results, tool_hashes = validate_nlr_root(
        args.nlr_root,
        bundle_samples,
        expected_units=args.expected_units,
        expected_workers=args.expected_worker_threads,
        expected_tools=expected_tools,
    )
    metadata_rows, metadata, cohort = read_metadata(args.unit_metadata, args.expected_units)
    targets = [row for row in selected if row["analysis_role"] == "target_repertoire"]
    if set(metadata) != {row["sample_id"] for row in targets}:
        raise ValueError("Unit metadata and NLR target sets differ")
    for row in targets:
        if metadata[row["sample_id"]]["biological_species"] != row["species"]:
            raise ValueError(f"Species mismatch for NLR target {row['sample_id']}")

    reference_row = next(row for row in selected if row["analysis_role"] == "reference_callable")
    reference_fasta = args.input_bundle / reference_row["relative_fasta"]
    if sha256_file(reference_fasta) != bundle_samples[reference_row["sample_id"]]["sha256"]:
        raise ValueError("Reference CDS FASTA differs from its frozen plain-input bundle checksum")
    reference_order = read_fasta_ids(reference_fasta)
    if len(reference_order) != args.expected_reference_genes:
        raise ValueError(
            f"Observed {len(reference_order)} reference CDS IDs; expected {args.expected_reference_genes}"
        )
    reference = set(reference_order)
    reference_calls = nlr_results[reference_row["sample_id"]]
    reference_nlr = set(reference_calls["sequence_ids"])
    if not reference_nlr or not reference_nlr.issubset(reference):
        raise ValueError("Reference NLR call IDs are empty or outside the reference-CDS universe")
    shared = read_shared(args.shared_positive_genes, reference, args.expected_shared_positive)
    nonshared_nlr = reference_nlr - shared
    if not nonshared_nlr:
        raise ValueError("No non-shared reference NLR genes remain after shared-gene exclusion")
    universe_digest = hashlib.sha256(
        ("\n".join(sorted(nonshared_nlr)) + "\n").encode("utf-8")
    ).hexdigest()
    universe_id = f"clem_scandens_nonshared_nlr_v1_{universe_digest[:12]}"

    loss_sha = sha256_file(args.loss_matrix)
    unit_counts = defaultdict(int)
    seen_pairs: set[tuple[str, str]] = set()
    callable_ids: dict[str, list[str]] = {unit: [] for unit in metadata}
    positive_ids: dict[str, list[str]] = {unit: [] for unit in metadata}
    with args.loss_matrix.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(LOSS_FIELDS[:5]) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Primary loss matrix lacks required columns: {sorted(missing)}")
        for line_number, raw in enumerate(reader, 2):
            if None in raw:
                raise ValueError(f"Extra primary loss-matrix fields at line {line_number}")
            row = {key: (value or "").strip() for key, value in raw.items()}
            gene, unit = row["reference_gene_id"], row["assembly_unit_id"]
            if gene not in reference or unit not in metadata:
                raise ValueError(f"Primary loss row outside exact reference/unit universe at line {line_number}")
            pair = (gene, unit)
            if pair in seen_pairs:
                raise ValueError(f"Duplicate primary loss pair at line {line_number}: {pair}")
            seen_pairs.add(pair)
            unit_counts[unit] += 1
            classification, callable_value = row["classification"], row["callable"]
            if classification not in ALLOWED_CLASSES or callable_value not in {"true", "false"}:
                raise ValueError(f"Invalid primary loss state at line {line_number}")
            if classification in {"deleted", "pseudogenized", "retained"} and callable_value != "true":
                raise ValueError(f"Resolved primary loss state is non-callable at line {line_number}")
            if gene in nonshared_nlr and classification in {"deleted", "pseudogenized", "retained"}:
                callable_ids[unit].append(gene)
                if classification in {"deleted", "pseudogenized"}:
                    positive_ids[unit].append(gene)
    expected_rows = args.expected_units * args.expected_reference_genes
    if len(seen_pairs) != expected_rows or set(unit_counts) != set(metadata):
        raise ValueError(f"Primary loss matrix is not the exact {expected_rows}-row grid")
    if any(count != args.expected_reference_genes for count in unit_counts.values()):
        raise ValueError("At least one primary loss unit lacks exact reference-gene closure")

    for unit in metadata:
        if len(callable_ids[unit]) != len(set(callable_ids[unit])):
            raise ValueError(f"Duplicate callable NLR denominator IDs for {unit}")
        if not callable_ids[unit]:
            raise ValueError(
                f"Callable NLR catalog is empty for {unit}; catalog-mode summary cannot "
                "represent an absent assembly-unit row"
            )
        if not set(positive_ids[unit]).issubset(callable_ids[unit]):
            raise ValueError(f"Positive NLR loss calls exceed the callable catalog for {unit}")

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", suffix=".tmp", dir=args.output_dir.parent))
    try:
        write_tsv(staging / "assembly_units.tsv", METADATA_FIELDS, metadata_rows)
        write_tsv(
            staging / "reference_nlr_universe.tsv",
            (
                "reference_nlr_id", "included_in_nonshared_analysis", "exclusion_reason",
                "reference_nlr_universe_id", "reference_calls_basename", "reference_calls_sha256",
            ),
            (
                {
                    "reference_nlr_id": gene,
                    "included_in_nonshared_analysis": str(gene in nonshared_nlr).lower(),
                    "exclusion_reason": "" if gene in nonshared_nlr else "shared_positive_complete",
                    "reference_nlr_universe_id": universe_id,
                    "reference_calls_basename": "nlr_calls.txt",
                    "reference_calls_sha256": reference_calls["nlr_calls_sha256"],
                }
                for gene in sorted(reference_nlr)
            ),
        )
        write_tsv(
            staging / "repertoire_counts.tsv",
            (
                "assembly_unit_id", "assembly_scope", "total_nlr_count",
                "repertoire_source_basename", "repertoire_source_sha256",
            ),
            (
                {
                    "assembly_unit_id": unit,
                    "assembly_scope": metadata[unit]["assembly_scope"],
                    "total_nlr_count": nlr_results[unit]["nlr_rows"],
                    "repertoire_source_basename": "nlr_calls.txt",
                    "repertoire_source_sha256": nlr_results[unit]["nlr_calls_sha256"],
                }
                for unit in metadata
            ),
        )
        write_tsv(
            staging / "positive_reference_nlr_loss_calls.tsv",
            ("assembly_unit_id", "assembly_scope", "reference_nlr_id", "reference_nlr_universe_id"),
            (
                {
                    "assembly_unit_id": unit,
                    "assembly_scope": metadata[unit]["assembly_scope"],
                    "reference_nlr_id": gene,
                    "reference_nlr_universe_id": universe_id,
                }
                for unit in metadata for gene in sorted(positive_ids[unit])
            ),
        )
        write_tsv(
            staging / "callable_reference_nlr_denominators.tsv",
            (
                "assembly_unit_id", "assembly_scope", "reference_nlr_id",
                "reference_nlr_universe_id", "denominator_source_basename",
                "denominator_source_sha256",
            ),
            (
                {
                    "assembly_unit_id": unit,
                    "assembly_scope": metadata[unit]["assembly_scope"],
                    "reference_nlr_id": gene,
                    "reference_nlr_universe_id": universe_id,
                    "denominator_source_basename": args.loss_matrix.name,
                    "denominator_source_sha256": loss_sha,
                }
                for unit in metadata for gene in sorted(callable_ids[unit])
            ),
        )
        sample_audit_rows = []
        for index, row in enumerate(selected, 1):
            result = nlr_results[row["sample_id"]]
            sample_audit_rows.append(
                {
                    "manifest_order": index,
                    "sample_id": row["sample_id"],
                    "analysis_role": row["analysis_role"],
                    "input_scope": row["input_scope"],
                    "input_fasta_sha256": bundle_samples[row["sample_id"]]["sha256"],
                    "nlr_calls_sha256": result["nlr_calls_sha256"],
                    "nlr_output_rows": result["nlr_rows"],
                    "nlr_output_sequence_ids": len(result["sequence_ids"]),
                    "nlr_output_locus_ids": len(result["locus_ids"]),
                    "validation_status": "PASS",
                }
            )
        write_tsv(
            staging / "sample_validation.tsv",
            (
                "manifest_order", "sample_id", "analysis_role", "input_scope",
                "input_fasta_sha256", "nlr_calls_sha256", "nlr_output_rows",
                "nlr_output_sequence_ids", "nlr_output_locus_ids", "validation_status",
            ),
            sample_audit_rows,
        )
        validation = {
            "status": "PASS_PRIMARY_NLR_SUMMARY_INPUTS",
            "workflow": "primary_nonshared_nlr_summary_inputs",
            "analysis_cohort": cohort,
            "assembly_unit_count": args.expected_units,
            "reference_gene_count": len(reference),
            "shared_positive_complete_gene_count": len(shared),
            "reference_nlr_gene_count": len(reference_nlr),
            "shared_reference_nlr_gene_count_excluded": len(reference_nlr & shared),
            "nonshared_reference_nlr_gene_count": len(nonshared_nlr),
            "reference_nlr_universe_id": universe_id,
            "total_target_nlr_locus_count": sum(int(nlr_results[unit]["nlr_rows"]) for unit in metadata),
            "positive_reference_nlr_loss_call_count": sum(map(len, positive_ids.values())),
            "callable_reference_nlr_denominator_sum": sum(map(len, callable_ids.values())),
            "denominator_policy": (
                "Each non-shared reference-NLR gene with callable=true in the complete primary "
                "loss matrix contributes one resolved unit comparison; uncertain calls are "
                "excluded even when the local interval was callable."
            ),
            "tool_sha256": tool_hashes,
            "input_bundle_sha256": bundle_hashes,
            "checks": {
                "completed_atomic_nlr_batch": "PASS",
                "sample_input_tool_output_checksum_closure": "PASS",
                "exact_reference_and_23_unit_cohort": "PASS",
                "complete_23_by_35547_loss_grid": "PASS",
                "shared_gene_exclusion": "PASS",
                "positive_calls_subset_of_callable_catalog": "PASS",
            },
        }
        (staging / "validation.json").write_text(
            json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        input_paths = [
            args.nlr_root / "selected_inputs.tsv",
            args.nlr_root / "batch_metadata.tsv",
            args.input_bundle / "run_manifest.json",
            args.input_bundle / "checksums.tsv",
            args.unit_metadata,
            args.loss_matrix,
            args.shared_positive_genes,
            args.nlr_root / reference_row["sample_id"] / "nlr_calls.txt",
        ]
        write_tsv(
            staging / "input_checksums.tsv",
            ("role", "basename", "bytes", "sha256"),
            (
                {
                    "role": role,
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                for role, path in zip(
                    (
                        "nlr_selected_inputs", "nlr_batch_metadata", "input_bundle_manifest",
                        "input_bundle_checksums", "unit_metadata", "primary_loss_matrix",
                        "shared_positive_genes", "reference_nlr_calls",
                    ),
                    input_paths,
                    strict=True,
                )
            ),
        )
        output_files = sorted(path for path in staging.iterdir() if path.name != "output_checksums.tsv")
        write_tsv(
            staging / "output_checksums.tsv",
            ("file", "bytes", "sha256"),
            (
                {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
                for path in output_files
            ),
        )
        for path in staging.iterdir():
            os.chmod(path, 0o644)
        os.chmod(staging, 0o755)
        os.replace(staging, args.output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return validation


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validation = run(args)
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        f"{validation['status']}\tunits={validation['assembly_unit_count']}\t"
        f"nonshared_reference_nlr={validation['nonshared_reference_nlr_gene_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
