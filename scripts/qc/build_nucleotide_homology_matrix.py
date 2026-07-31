#!/usr/bin/env python3
"""Build one fail-closed 29x29 nucleotide chromosome-homology matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.chromosome_assignment import (  # noqa: E402
    NUCLEOTIDE_COLUMNS,
    REQUIRED_PAF_TAGS,
    load_homology_policy,
    read_score_matrix,
)
from geneloss_repro.chromosome_provenance import (  # noqa: E402
    NUCLEOTIDE_UPSTREAM_CHECKS,
    PROVENANCE_SCHEMA_VERSION,
    capture_snapshot,
    read_reference_registries,
    read_target_registry,
    validate_matrix_provenance,
)
from geneloss_repro.nucleotide_matrix import (  # noqa: E402
    NucleotideMatrixError,
    build_nucleotide_rows,
    fasta_lengths,
    read_role_normalized_paf,
)


ORIENTATION_COLUMNS = (
    "query_chromosome",
    "reference_chromosome",
    "plus_matching_bases",
    "minus_matching_bases",
    "total_matching_bases",
    "dominant_orientation",
    "dominant_fraction",
    "automatic_orientation_gate",
)


class BuildError(RuntimeError):
    pass


def stable_binding(path: Path) -> dict[str, str | int]:
    source = path.expanduser().resolve()
    if not source.is_file() or source.is_symlink():
        raise BuildError(f"Input must be a regular non-symlink file: {source.name}")
    before = source.stat()
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    after = source.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise BuildError(f"Input changed while hashing: {source.name}")
    if after.st_size <= 0:
        raise BuildError(f"Input is empty: {source.name}")
    return {"basename": source.name, "bytes": after.st_size, "sha256": digest.hexdigest()}


def require_binding(expected: object, observed: dict[str, str | int], label: str) -> None:
    if not isinstance(expected, dict) or {
        key: expected.get(key) for key in ("basename", "bytes", "sha256")
    } != observed:
        raise BuildError(f"{label} does not match its frozen binding")


def strict_json(path: Path) -> dict[str, object]:
    try:
        document = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=lambda value: (_ for _ in ()).throw(
                BuildError(f"Non-finite JSON constant: {value}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise BuildError(f"Cannot read JSON {path.name}: {error}") from error
    if not isinstance(document, dict):
        raise BuildError(f"{path.name}: top-level JSON must be an object")
    return document


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, document: dict[str, object]) -> None:
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def generation_parameters(policy) -> dict[str, object]:
    return {
        "minimap2_version": policy.minimap2_version,
        "minimap2_command_template": list(policy.minimap2_command_template),
        "required_paf_tags": list(REQUIRED_PAF_TAGS),
        "minimap2_primary_only": policy.minimap2_primary_only,
        "minimum_mapq": policy.minimum_mapq,
        "minimum_alignment_block_bp": policy.minimum_alignment_block_bp,
        "maximum_de": policy.maximum_de,
        "coverage_arithmetic": policy.coverage_arithmetic,
        "nucleotide_score_formula": policy.nucleotide_score_formula,
        "reciprocal_nucleotide_coverage_formula": (
            policy.reciprocal_nucleotide_coverage_formula
        ),
        "arithmetic_tolerance": policy.arithmetic_tolerance,
    }


def role_spec(matrix_role: str, policy) -> tuple[str, str, str, str, str]:
    if matrix_role == "nucleotide_hy4a":
        return (
            "HY4A",
            policy.coordinate_reference,
            "primary_coordinate_reference",
            "target_to_hy4a",
            "hy4a_to_target",
        )
    if matrix_role == "nucleotide_hy4p":
        return (
            "HY4P",
            policy.confirmation_reference,
            "independent_confirmation_reference",
            "target_to_hy4p",
            "hy4p_to_target",
        )
    raise BuildError("matrix role must be nucleotide_hy4a or nucleotide_hy4p")


def run(args: argparse.Namespace) -> Path:
    policy = load_homology_policy(args.parameters)
    slot, reference_id, reference_role, forward_role, reverse_role = role_spec(
        args.matrix_role, policy
    )
    target_registry = read_target_registry(capture_snapshot(args.target_asset_registry))
    references = read_reference_registries(
        capture_snapshot(args.reference_asset_registry),
        capture_snapshot(args.reference_chromosome_map_registry),
        coordinate_reference=policy.coordinate_reference,
        confirmation_reference=policy.confirmation_reference,
        expected_chromosomes=policy.expected_reference_chromosomes,
    )
    if references.asset_snapshot.sha256 != policy.reference_asset_registry_sha256:
        raise BuildError("Reference asset registry hash conflicts with frozen policy")
    if references.map_snapshot.sha256 != policy.reference_chromosome_map_registry_sha256:
        raise BuildError("Reference chromosome-map registry hash conflicts with frozen policy")
    target_key = (args.assembly_unit_id, args.target_scope_id)
    if target_key not in target_registry.assets:
        raise BuildError("Target registry lacks the invoked unit/scope")

    target_genome_binding = stable_binding(args.target_genome)
    reference_genome_binding = stable_binding(args.reference_genome)
    require_binding(
        target_registry.assets[target_key]["genome"].as_dict(),
        target_genome_binding,
        "target genome",
    )
    require_binding(
        references.assets[reference_id]["genome"].as_dict(),
        reference_genome_binding,
        "reference genome",
    )
    forward_binding = stable_binding(args.forward_paf)
    reverse_binding = stable_binding(args.reverse_paf)

    bundle = strict_json(args.bundle_validation)
    if (
        bundle.get("status") != "PASS"
        or bundle.get("workflow") != "bidirectional_chromosome_minimap_bundle_validation"
        or bundle.get("unit") != args.assembly_unit_id
    ):
        raise BuildError("Bundle validation is not the expected PASS for this unit")
    bundle_inputs = bundle.get("inputs")
    bundle_comparisons = bundle.get("comparisons")
    if not isinstance(bundle_inputs, dict) or not isinstance(bundle_comparisons, dict):
        raise BuildError("Bundle validation lacks exact inputs/comparisons")
    require_binding(bundle_inputs.get("target"), target_genome_binding, "bundle target genome")
    reference_bundle_role = "hy4a" if slot == "HY4A" else "hy4p"
    require_binding(
        bundle_inputs.get(reference_bundle_role),
        reference_genome_binding,
        f"bundle {reference_bundle_role} genome",
    )
    for role, binding in ((forward_role, forward_binding), (reverse_role, reverse_binding)):
        comparison = bundle_comparisons.get(role)
        if not isinstance(comparison, dict):
            raise BuildError(f"Bundle validation lacks {role}")
        require_binding(comparison.get("paf"), binding, f"bundle {role} PAF")

    target_lengths = fasta_lengths(args.target_genome)
    reference_lengths = fasta_lengths(args.reference_genome)
    if len(target_lengths) != policy.expected_query_chromosomes:
        raise BuildError("Target genome must contain exactly 29 chromosomes")
    if len(reference_lengths) != policy.expected_reference_chromosomes:
        raise BuildError("Reference genome must contain exactly 29 chromosomes")
    canonical = references.chromosome_maps[reference_id]
    forward_records, forward_audit = read_role_normalized_paf(
        args.forward_paf,
        query_role="target",
        target_lengths=target_lengths,
        reference_lengths=reference_lengths,
        minimum_mapq=policy.minimum_mapq,
        minimum_alignment_block_bp=policy.minimum_alignment_block_bp,
        maximum_de=Decimal(str(policy.maximum_de)),
    )
    reverse_records, reverse_audit = read_role_normalized_paf(
        args.reverse_paf,
        query_role="reference",
        target_lengths=target_lengths,
        reference_lengths=reference_lengths,
        minimum_mapq=policy.minimum_mapq,
        minimum_alignment_block_bp=policy.minimum_alignment_block_bp,
        maximum_de=Decimal(str(policy.maximum_de)),
    )
    rows, orientation_rows = build_nucleotide_rows(
        target_lengths=target_lengths,
        reference_lengths=reference_lengths,
        canonical_by_reference=canonical,
        forward_records=forward_records,
        reverse_records=reverse_records,
    )

    output = args.output_dir.expanduser().resolve()
    if output.exists():
        raise BuildError(f"Refusing to overwrite output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        matrix_path = staging / f"{args.matrix_role}.tsv"
        orientation_path = staging / f"{args.matrix_role}.orientation_support.tsv"
        report_path = staging / f"{args.matrix_role}.upstream_validation.json"
        sidecar_path = staging / f"{args.matrix_role}.provenance.json"
        audit_path = staging / f"{args.matrix_role}.input_audit.json"
        write_tsv(matrix_path, NUCLEOTIDE_COLUMNS, rows)
        write_tsv(orientation_path, ORIENTATION_COLUMNS, orientation_rows)
        matrix_binding = stable_binding(matrix_path)
        generation = generation_parameters(policy)
        report = {
            "schema_version": "1.0.0",
            "workflow": "chromosome_nucleotide_matrix_builder",
            "workflow_version": "1.0.0",
            "status": "PASS",
            "assembly_unit_id": args.assembly_unit_id,
            "target_scope_id": args.target_scope_id,
            "matrix_role": args.matrix_role,
            "matrix_kind": "nucleotide",
            "reference_slot": slot,
            "reference_id": reference_id,
            "reference_role": reference_role,
            "reference_map_id": references.map_ids[reference_id],
            "matrix_sha256": matrix_binding["sha256"],
            "target_asset_registry_sha256": target_registry.snapshot.sha256,
            "reference_asset_registry_sha256": references.asset_snapshot.sha256,
            "reference_chromosome_map_registry_sha256": references.map_snapshot.sha256,
            "generation_parameters": generation,
            "checks": {key: True for key in sorted(NUCLEOTIDE_UPSTREAM_CHECKS)},
        }
        write_json(report_path, report)
        report_binding = stable_binding(report_path)
        sidecar = {
            "schema_version": PROVENANCE_SCHEMA_VERSION,
            "status": "PASS",
            "assembly_unit_id": args.assembly_unit_id,
            "target_scope_id": args.target_scope_id,
            "matrix_role": args.matrix_role,
            "matrix_kind": "nucleotide",
            "reference_slot": slot,
            "reference_id": reference_id,
            "reference_role": reference_role,
            "reference_map_id": references.map_ids[reference_id],
            "matrix": matrix_binding,
            "target_asset_registry": target_registry.snapshot.public_binding(),
            "target_assets": {
                role: binding.as_dict()
                for role, binding in target_registry.assets[target_key].items()
            },
            "reference_asset_registry": references.asset_snapshot.public_binding(),
            "reference_assets": {
                role: binding.as_dict()
                for role, binding in references.assets[reference_id].items()
            },
            "reference_chromosome_map_registry": references.map_snapshot.public_binding(),
            "generation_parameters": generation,
            "upstream_validation_report": report_binding,
        }
        write_json(sidecar_path, sidecar)
        write_json(
            audit_path,
            {
                "schema_version": 1,
                "status": "PASS",
                "assembly_unit_id": args.assembly_unit_id,
                "target_scope_id": args.target_scope_id,
                "matrix_role": args.matrix_role,
                "bundle_validation": stable_binding(args.bundle_validation),
                "inputs": {
                    "target_genome": target_genome_binding,
                    "reference_genome": reference_genome_binding,
                    "forward_paf": forward_binding,
                    "reverse_paf": reverse_binding,
                },
                "paf_filter_counts": {
                    "forward": vars(forward_audit),
                    "reverse": vars(reverse_audit),
                },
                "matrix": matrix_binding,
                "orientation_support": stable_binding(orientation_path),
            },
        )

        parsed_matrix = read_score_matrix(
            capture_snapshot(matrix_path), role=args.matrix_role, kind="nucleotide", policy=policy
        )
        if len(parsed_matrix.normalized_rows) != 841:
            raise BuildError("Published matrix did not revalidate as exactly 841 rows")
        validate_matrix_provenance(
            sidecar=capture_snapshot(sidecar_path),
            matrix=capture_snapshot(matrix_path),
            assembly_unit_id=args.assembly_unit_id,
            target_scope_id=args.target_scope_id,
            matrix_role=args.matrix_role,
            matrix_kind="nucleotide",
            reference_slot=slot,
            reference_id=reference_id,
            reference_role=reference_role,
            reference_map_id=references.map_ids[reference_id],
            expected_generation_parameters=generation,
            target_registry=target_registry,
            reference_registries=references,
        )
        os.rename(staging, output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return output


def parser() -> argparse.ArgumentParser:
    argument_parser = argparse.ArgumentParser(description=__doc__)
    argument_parser.add_argument("--assembly-unit-id", required=True)
    argument_parser.add_argument("--target-scope-id", required=True)
    argument_parser.add_argument(
        "--matrix-role", required=True, choices=("nucleotide_hy4a", "nucleotide_hy4p")
    )
    argument_parser.add_argument("--bundle-validation", required=True, type=Path)
    argument_parser.add_argument("--forward-paf", required=True, type=Path)
    argument_parser.add_argument("--reverse-paf", required=True, type=Path)
    argument_parser.add_argument("--target-genome", required=True, type=Path)
    argument_parser.add_argument("--reference-genome", required=True, type=Path)
    argument_parser.add_argument("--parameters", required=True, type=Path)
    argument_parser.add_argument("--target-asset-registry", required=True, type=Path)
    argument_parser.add_argument("--reference-asset-registry", required=True, type=Path)
    argument_parser.add_argument(
        "--reference-chromosome-map-registry", required=True, type=Path
    )
    argument_parser.add_argument("--output-dir", required=True, type=Path)
    return argument_parser


def main() -> int:
    try:
        output = run(parser().parse_args())
        print(f"PASS\t{output}")
        return 0
    except (BuildError, NucleotideMatrixError, OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
