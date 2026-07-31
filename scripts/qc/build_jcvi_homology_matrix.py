#!/usr/bin/env python3
"""Build one fail-closed bidirectional 29x29 JCVI homology matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shlex
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.chromosome_assignment import (  # noqa: E402
    JCVI_COLUMNS,
    load_homology_policy,
    read_score_matrix,
)
from geneloss_repro.chromosome_provenance import (  # noqa: E402
    JCVI_UPSTREAM_CHECKS,
    PROVENANCE_SCHEMA_VERSION,
    capture_snapshot,
    read_reference_registries,
    read_target_registry,
    validate_matrix_provenance,
)
from geneloss_repro.jcvi_matrix import (  # noqa: E402
    JcviMatrixError,
    build_jcvi_rows,
    read_bed,
    read_normalized_anchor_pairs,
    relabel_reference_bed_from_canonical_truth,
    require_bed_gff_identity,
    require_bed_protein_identity,
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


def require_binding(expected: object, observed: dict[str, str | int], label: str) -> None:
    if not isinstance(expected, dict):
        raise BuildError(f"{label} lacks a binding object")
    normalized = {
        "basename": expected.get("basename", Path(str(expected.get("path", ""))).name),
        "bytes": expected.get("bytes", expected.get("size_bytes")),
        "sha256": expected.get("sha256"),
    }
    if normalized != observed:
        raise BuildError(f"{label} does not match its frozen binding")


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=JCVI_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, document: dict[str, object]) -> None:
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def role_spec(matrix_role: str, policy) -> tuple[str, str, str]:
    if matrix_role == "jcvi_hy4a":
        return "HY4A", policy.coordinate_reference, "primary_coordinate_reference"
    if matrix_role == "jcvi_hy4p":
        return "HY4P", policy.confirmation_reference, "independent_confirmation_reference"
    raise BuildError("matrix role must be jcvi_hy4a or jcvi_hy4p")


def validate_canonical_manifest(
    manifest_path: Path,
    *,
    reference_id: str,
    reference_assets,
    reference_protein: dict[str, str | int],
    reference_bed: dict[str, str | int],
) -> dict[str, str | int]:
    manifest_binding = stable_binding(manifest_path)
    manifest = strict_json(manifest_path)
    if manifest.get("status") != "PASS":
        raise BuildError("Reference canonical-input manifest is not PASS")
    inputs, outputs, checks = (
        manifest.get("inputs"),
        manifest.get("outputs"),
        manifest.get("checks"),
    )
    if not isinstance(inputs, dict) or not isinstance(outputs, dict) or not isinstance(checks, dict):
        raise BuildError("Reference canonical-input manifest is incomplete")
    require_binding(inputs.get("protein"), reference_assets[reference_id]["protein"].as_dict(),
                    "canonical source protein")
    require_binding(inputs.get("gff"), reference_assets[reference_id]["gff"].as_dict(),
                    "canonical source GFF3")
    require_binding(outputs.get("protein"), reference_protein, "canonical output protein")
    require_binding(outputs.get("bed"), reference_bed, "canonical output BED")
    required_checks = {
        "canonical_protein_ids_unique",
        "protein_gff_transcript_id_identity",
        "one_selected_transcript_per_gene",
        "complete_gene_closure",
    }
    if any(checks.get(key) is not True for key in required_checks):
        raise BuildError("Reference canonical-input checks did not all pass")
    return manifest_binding


def validate_run(
    run_dir: Path,
    *,
    expected_reference_id: str,
    expected_query_id: str,
    expected_inputs: dict[str, dict[str, str | int]],
) -> dict[str, object]:
    root = run_dir.expanduser().resolve()
    manifest_path = root / "run_manifest.json"
    summary_json_path = root / "jcvi_bidirectional_coverage.json"
    summary_tsv_path = root / "jcvi_bidirectional_coverage.tsv"
    manifest_binding = stable_binding(manifest_path)
    summary_json_binding = stable_binding(summary_json_path)
    summary_tsv_binding = stable_binding(summary_tsv_path)
    manifest = strict_json(manifest_path)
    if manifest.get("reference_id") != expected_reference_id:
        raise BuildError("JCVI run reference alias is wrong")
    if manifest.get("query_id") != expected_query_id:
        raise BuildError("JCVI run query alias is wrong")
    if manifest.get("threads") != 4:
        raise BuildError("JCVI run did not record exactly four threads")
    inputs = manifest.get("inputs")
    commands = manifest.get("commands")
    if not isinstance(inputs, dict) or not isinstance(commands, dict):
        raise BuildError("JCVI run manifest lacks inputs or commands")
    for role, expected in expected_inputs.items():
        require_binding(inputs.get(role), expected, f"JCVI run {role}")
    ortholog = commands.get("ortholog")
    if not isinstance(ortholog, str):
        raise BuildError("JCVI run manifest lacks ortholog command")
    tokens = shlex.split(ortholog)
    required_tokens = {
        "--dbtype=prot",
        "--align_soft=last",
        "--cpus=4",
        "--cscore=0.7",
        "--tandem_Nmax=10",
        "--dist=20",
        "--min_size=4",
        "--no_strip_names",
        "--no_dotplot",
    }
    if not required_tokens.issubset(tokens):
        raise BuildError("JCVI ortholog command differs from the frozen parameterization")
    anchors_path = root / "work" / f"{expected_reference_id}.{expected_query_id}.anchors"
    anchors_binding = stable_binding(anchors_path)
    summary = strict_json(summary_json_path)
    if summary.get("script_version") != "2.0.0":
        raise BuildError("JCVI coverage summary has the wrong script version")
    summary_inputs = summary.get("inputs")
    if not isinstance(summary_inputs, dict):
        raise BuildError("JCVI coverage summary lacks input bindings")
    for role, expected in expected_inputs.items():
        require_binding(summary_inputs.get(role), expected, f"JCVI summary {role}")
    require_binding(summary_inputs.get("anchors"), anchors_binding, "JCVI summary anchors")
    return {
        "run_manifest": manifest_binding,
        "anchors": anchors_binding,
        "summary_json": summary_json_binding,
        "summary_tsv": summary_tsv_binding,
        "anchors_path": anchors_path,
    }


def generation_parameters(policy) -> dict[str, object]:
    return {
        "aligner": policy.jcvi_aligner,
        "database_type": policy.jcvi_database_type,
        "cscore": policy.jcvi_cscore,
        "tandem_nmax": policy.jcvi_tandem_nmax,
        "maximum_gene_distance": policy.jcvi_maximum_gene_distance,
        "minimum_anchor_block_size": policy.jcvi_minimum_anchor_block_size,
        "coverage_anchor_source": policy.jcvi_coverage_anchor_source,
        "jcvi_anchor_counting": policy.jcvi_anchor_counting,
        "jcvi_score_formula": policy.jcvi_score_formula,
        "arithmetic_tolerance": policy.arithmetic_tolerance,
    }


def run(args: argparse.Namespace) -> Path:
    policy = load_homology_policy(args.parameters)
    slot, reference_id, reference_role = role_spec(args.matrix_role, policy)
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

    target_protein_binding = stable_binding(args.target_protein)
    target_gff_binding = stable_binding(args.target_gff)
    target_bed_binding = stable_binding(args.target_bed)
    reference_protein_binding = stable_binding(args.reference_protein)
    reference_bed_binding = stable_binding(args.reference_bed)
    require_binding(
        target_registry.assets[target_key]["protein"].as_dict(),
        target_protein_binding,
        "target protein",
    )
    require_binding(
        target_registry.assets[target_key]["gff"].as_dict(), target_gff_binding, "target GFF3"
    )
    canonical_manifest_binding = validate_canonical_manifest(
        args.reference_canonical_manifest,
        reference_id=reference_id,
        reference_assets=references.assets,
        reference_protein=reference_protein_binding,
        reference_bed=reference_bed_binding,
    )

    target_bed = read_bed(args.target_bed)
    reference_bed_publisher = read_bed(args.reference_bed)
    if len(target_bed.chromosomes) != policy.expected_query_chromosomes:
        raise BuildError("Target BED must contain exactly 29 chromosomes")
    if len(reference_bed_publisher.chromosomes) != policy.expected_reference_chromosomes:
        raise BuildError("Reference BED must contain exactly 29 chromosomes")
    require_bed_protein_identity(target_bed, args.target_protein, label="target")
    require_bed_protein_identity(
        reference_bed_publisher, args.reference_protein, label="reference"
    )
    require_bed_gff_identity(target_bed, args.target_gff)
    frozen_reference_map = references.chromosome_maps[reference_id]
    if set(reference_bed_publisher.chromosomes) == set(frozen_reference_map):
        reference_bed = reference_bed_publisher
        reference_bed_aliases = {chromosome: chromosome for chromosome in reference_bed.chromosomes}
    else:
        reference_bed, reference_bed_aliases = relabel_reference_bed_from_canonical_truth(
            reference_bed_publisher, frozen_reference_map
        )

    common_inputs = {
        "target": {
            "protein": target_protein_binding,
            "bed": target_bed_binding,
        },
        "reference": {
            "protein": reference_protein_binding,
            "bed": reference_bed_binding,
        },
    }
    forward = validate_run(
        args.forward_run_dir,
        expected_reference_id=slot,
        expected_query_id=args.target_alias,
        expected_inputs={
            "reference_protein": common_inputs["reference"]["protein"],
            "reference_bed": common_inputs["reference"]["bed"],
            "query_protein": common_inputs["target"]["protein"],
            "query_bed": common_inputs["target"]["bed"],
        },
    )
    reverse = validate_run(
        args.reverse_run_dir,
        expected_reference_id=args.target_alias,
        expected_query_id=slot,
        expected_inputs={
            "reference_protein": common_inputs["target"]["protein"],
            "reference_bed": common_inputs["target"]["bed"],
            "query_protein": common_inputs["reference"]["protein"],
            "query_bed": common_inputs["reference"]["bed"],
        },
    )
    forward_pairs, forward_audit = read_normalized_anchor_pairs(
        Path(forward["anchors_path"]),
        first_bed=reference_bed,
        second_bed=target_bed,
        first_role="reference",
    )
    reverse_pairs, reverse_audit = read_normalized_anchor_pairs(
        Path(reverse["anchors_path"]),
        first_bed=target_bed,
        second_bed=reference_bed,
        first_role="target",
    )
    rows, pair_audit = build_jcvi_rows(
        target_bed=target_bed,
        reference_bed=reference_bed,
        canonical_by_reference=references.chromosome_maps[reference_id],
        forward_pairs=forward_pairs,
        reverse_pairs=reverse_pairs,
    )

    output = args.output_dir.expanduser().resolve()
    if output.exists():
        raise BuildError(f"Refusing to overwrite output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        matrix_path = staging / f"{args.matrix_role}.tsv"
        report_path = staging / f"{args.matrix_role}.upstream_validation.json"
        sidecar_path = staging / f"{args.matrix_role}.provenance.json"
        audit_path = staging / f"{args.matrix_role}.input_audit.json"
        write_tsv(matrix_path, rows)
        matrix_binding = stable_binding(matrix_path)
        generation = generation_parameters(policy)
        report = {
            "schema_version": "1.0.0",
            "workflow": "chromosome_jcvi_matrix_builder",
            "workflow_version": "1.0.0",
            "status": "PASS",
            "assembly_unit_id": args.assembly_unit_id,
            "target_scope_id": args.target_scope_id,
            "matrix_role": args.matrix_role,
            "matrix_kind": "jcvi",
            "reference_slot": slot,
            "reference_id": reference_id,
            "reference_role": reference_role,
            "reference_map_id": references.map_ids[reference_id],
            "matrix_sha256": matrix_binding["sha256"],
            "target_asset_registry_sha256": target_registry.snapshot.sha256,
            "reference_asset_registry_sha256": references.asset_snapshot.sha256,
            "reference_chromosome_map_registry_sha256": references.map_snapshot.sha256,
            "generation_parameters": generation,
            "checks": {key: True for key in sorted(JCVI_UPSTREAM_CHECKS)},
        }
        write_json(report_path, report)
        report_binding = stable_binding(report_path)
        sidecar = {
            "schema_version": PROVENANCE_SCHEMA_VERSION,
            "status": "PASS",
            "assembly_unit_id": args.assembly_unit_id,
            "target_scope_id": args.target_scope_id,
            "matrix_role": args.matrix_role,
            "matrix_kind": "jcvi",
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
                "inputs": {
                    "target_protein": target_protein_binding,
                    "target_gff": target_gff_binding,
                    "target_bed": target_bed_binding,
                    "reference_protein": reference_protein_binding,
                    "reference_bed": reference_bed_binding,
                    "reference_canonical_manifest": canonical_manifest_binding,
                    "forward_run": {
                        key: value for key, value in forward.items() if key != "anchors_path"
                    },
                    "reverse_run": {
                        key: value for key, value in reverse.items() if key != "anchors_path"
                    },
                },
                "anchor_audit": {
                    "forward": vars(forward_audit),
                    "reverse": vars(reverse_audit),
                    **pair_audit,
                },
                "reference_bed_chromosome_aliases": reference_bed_aliases,
                "matrix": matrix_binding,
            },
        )
        parsed = read_score_matrix(
            capture_snapshot(matrix_path), role=args.matrix_role, kind="jcvi", policy=policy
        )
        if len(parsed.normalized_rows) != 841:
            raise BuildError("Published JCVI matrix did not revalidate as exactly 841 rows")
        validate_matrix_provenance(
            sidecar=capture_snapshot(sidecar_path),
            matrix=capture_snapshot(matrix_path),
            assembly_unit_id=args.assembly_unit_id,
            target_scope_id=args.target_scope_id,
            matrix_role=args.matrix_role,
            matrix_kind="jcvi",
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
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assembly-unit-id", required=True)
    p.add_argument("--target-scope-id", required=True)
    p.add_argument("--target-alias", required=True)
    p.add_argument("--matrix-role", required=True, choices=("jcvi_hy4a", "jcvi_hy4p"))
    p.add_argument("--forward-run-dir", required=True, type=Path)
    p.add_argument("--reverse-run-dir", required=True, type=Path)
    p.add_argument("--target-protein", required=True, type=Path)
    p.add_argument("--target-gff", required=True, type=Path)
    p.add_argument("--target-bed", required=True, type=Path)
    p.add_argument("--reference-protein", required=True, type=Path)
    p.add_argument("--reference-bed", required=True, type=Path)
    p.add_argument("--reference-canonical-manifest", required=True, type=Path)
    p.add_argument("--parameters", required=True, type=Path)
    p.add_argument("--target-asset-registry", required=True, type=Path)
    p.add_argument("--reference-asset-registry", required=True, type=Path)
    p.add_argument("--reference-chromosome-map-registry", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    return p


def main() -> int:
    try:
        output = run(parser().parse_args())
        print(f"PASS\t{output}")
        return 0
    except (BuildError, JcviMatrixError, OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
