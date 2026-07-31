#!/usr/bin/env python3
"""Atomically publish a chromosome-labelled and direction-harmonized bundle."""

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
from decimal import Decimal, InvalidOperation
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 server fallback
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.chromosome_assignment import (  # noqa: E402
    NUCLEOTIDE_COLUMNS,
    load_homology_policy,
    read_score_matrix,
)
from geneloss_repro.chromosome_harmonization import (  # noqa: E402
    HarmonizationError,
    build_actions,
    read_fasta,
    transform_genome,
    transform_gff,
    validate_cds_protein_closure,
    validate_sequence_closure,
    write_fasta,
)
from geneloss_repro.chromosome_provenance import (  # noqa: E402
    capture_snapshot,
    read_target_registry,
)


SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
COMMIT_RE = re.compile(r"^[0-9a-f]{40}(?:[0-9a-f]{24})?$")
FINAL_MAP_COLUMNS = (
    "query_chromosome",
    "final_chromosome",
    "coordinate_reference",
    "confirmation_reference",
    "status",
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
ACTION_COLUMNS = (
    "source_chromosome",
    "final_chromosome",
    "orientation",
    "source_length_bp",
    "dominant_fraction",
    "oriented_matching_bases",
    "status",
)


class PublicationError(RuntimeError):
    pass


def stable_binding(path: Path) -> dict[str, str | int]:
    source = path.expanduser().resolve()
    if not source.is_file() or source.is_symlink():
        raise PublicationError(f"Input must be a regular non-symlink file: {source.name}")
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
        raise PublicationError(f"Input changed while hashing: {source.name}")
    if after.st_size <= 0:
        raise PublicationError(f"Input is empty: {source.name}")
    return {"basename": source.name, "bytes": after.st_size, "sha256": digest.hexdigest()}


def strict_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=lambda item: (_ for _ in ()).throw(
                PublicationError(f"Non-finite JSON constant: {item}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise PublicationError(f"Cannot read JSON {path.name}: {error}") from error
    if not isinstance(value, dict):
        raise PublicationError(f"{path.name}: top-level JSON must be an object")
    return value


def normalize_binding(value: object) -> dict[str, str | int]:
    if not isinstance(value, dict):
        raise PublicationError("Expected a file-binding object")
    result = {
        "basename": value.get("basename", value.get("file_name")),
        "bytes": value.get("bytes", value.get("size_bytes")),
        "sha256": value.get("sha256"),
    }
    if (
        not isinstance(result["basename"], str)
        or not isinstance(result["bytes"], int)
        or not isinstance(result["sha256"], str)
        or not SHA256_RE.fullmatch(result["sha256"])
    ):
        raise PublicationError("Malformed file-binding object")
    return result


def require_binding(expected: object, observed: dict[str, str | int], label: str) -> None:
    if normalize_binding(expected) != observed:
        raise PublicationError(f"{label} does not match its frozen binding")


def read_exact_tsv(path: Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != columns:
            raise PublicationError(f"{path.name}: exact TSV columns differ")
        rows = list(reader)
    if not rows:
        raise PublicationError(f"{path.name}: TSV has no rows")
    return rows


def validate_checksum_bundle(root: Path) -> dict[str, dict[str, str | int]]:
    checksum_path = root / "checksums.tsv"
    rows = read_exact_tsv(checksum_path, ("file", "bytes", "sha256"))
    expected_files: set[str] = set()
    bindings: dict[str, dict[str, str | int]] = {}
    for line_number, row in enumerate(rows, 2):
        name = row["file"]
        if not name or name != Path(name).name or name in expected_files:
            raise PublicationError(
                f"{checksum_path.name}:{line_number}: invalid/duplicate file name"
            )
        try:
            size = int(row["bytes"])
        except ValueError as error:
            raise PublicationError(
                f"{checksum_path.name}:{line_number}: non-integer byte count"
            ) from error
        if size <= 0 or not SHA256_RE.fullmatch(row["sha256"]):
            raise PublicationError(f"{checksum_path.name}:{line_number}: invalid binding")
        observed = stable_binding(root / name)
        expected = {"basename": name, "bytes": size, "sha256": row["sha256"]}
        if observed != expected:
            raise PublicationError(f"{root.name}/{name}: checksum binding mismatch")
        expected_files.add(name)
        bindings[name] = observed
    actual_files = {path.name for path in root.iterdir() if path.is_file()}
    if actual_files != expected_files.union({"checksums.tsv"}):
        raise PublicationError(f"{root.name}: checksum manifest file closure failed")
    return bindings


def validate_assignment_bundle(
    root: Path,
    *,
    assembly_unit_id: str,
    target_scope_id: str,
    trusted_commit: str,
    nucleotide_matrix_binding: dict[str, str | int],
    nucleotide_provenance_binding: dict[str, str | int],
) -> tuple[Path, dict[str, str | int]]:
    if not COMMIT_RE.fullmatch(trusted_commit):
        raise PublicationError("Trusted repository commit is not a full Git commit ID")
    checksums = validate_checksum_bundle(root)
    manifest = strict_json(root / "run_manifest.json")
    validation = strict_json(root / "validation.json")
    for document, label in ((manifest, "assignment manifest"), (validation, "validation")):
        if (
            document.get("workflow") != "chromosome_homology_assignment"
            or document.get("assembly_unit_id") != assembly_unit_id
            or document.get("target_scope_id") != target_scope_id
            or document.get("trusted_repository_commit") != trusted_commit
            or document.get("status") != "PASS_AUTO"
            or document.get("publication_gate") != "PASS"
            or document.get("failure_states") != []
        ):
            raise PublicationError(f"{label} is not the exact requested PASS bundle")
    checks = validation.get("checks")
    if not isinstance(checks, dict) or not checks or any(value is not True for value in checks.values()):
        raise PublicationError("Assignment validation checks did not all pass")
    inputs = manifest.get("inputs")
    if not isinstance(inputs, dict):
        raise PublicationError("Assignment manifest lacks inputs")
    roles = inputs.get("matrix_roles")
    if not isinstance(roles, dict) or not isinstance(roles.get("nucleotide_hy4a"), dict):
        raise PublicationError("Assignment manifest lacks nucleotide_hy4a role binding")
    role = roles["nucleotide_hy4a"]
    require_binding(role.get("matrix"), nucleotide_matrix_binding, "assigned nucleotide matrix")
    require_binding(
        role.get("provenance_sidecar"),
        nucleotide_provenance_binding,
        "assigned nucleotide provenance",
    )
    map_name = f"{assembly_unit_id}.final_chromosome_map.tsv"
    if map_name not in checksums:
        raise PublicationError("PASS assignment bundle lacks its final chromosome map")
    return root / map_name, checksums[map_name]


def validate_primary_annotation_bundle(
    root: Path,
    *,
    gff_binding: dict[str, str | int],
    cds_binding: dict[str, str | int],
    protein_binding: dict[str, str | int],
) -> dict[str, str | int]:
    checksums = validate_checksum_bundle(root)
    manifest_path = root / "run_manifest.json"
    manifest = strict_json(manifest_path)
    if (
        manifest.get("status") != "PASS"
        or manifest.get("publication_gate") != "PASS"
        or manifest.get("workflow") != "primary_annotation_standardization"
        or manifest.get("workflow_version") != "1.2.0"
    ):
        raise PublicationError("Strict primary-annotation bundle is not PASS")
    comparison = manifest.get("gffread_comparison")
    if not isinstance(comparison, dict) or comparison.get("status") != "PASS":
        raise PublicationError("Strict primary-annotation gffread gate is not PASS")
    require_binding(
        {
            "file_name": comparison.get("comparison_GFF3_file_name"),
            "bytes": comparison.get("comparison_GFF3_bytes"),
            "sha256": comparison.get("comparison_GFF3_sha256"),
        },
        gff_binding,
        "strict primary GFF3",
    )
    for observed, label in ((gff_binding, "GFF3"), (cds_binding, "CDS"), (protein_binding, "protein")):
        name = str(observed["basename"])
        if name not in checksums or checksums[name] != observed:
            raise PublicationError(f"Strict primary {label} is absent from exact checksums")
    return stable_binding(manifest_path)


def parse_exact_decimal(value: str, label: str) -> Decimal:
    try:
        parsed = Decimal(value)
    except InvalidOperation as error:
        raise PublicationError(f"{label} is not numeric") from error
    if not parsed.is_finite():
        raise PublicationError(f"{label} must be finite")
    return parsed


def load_orientation_actions(
    *,
    final_map_path: Path,
    matrix_path: Path,
    orientation_path: Path,
    policy,
    minimum_fraction: Decimal,
    minimum_matching_bases: int,
) -> tuple[dict[str, str], dict[str, str], dict[str, tuple[str, int]]]:
    final_rows = read_exact_tsv(final_map_path, FINAL_MAP_COLUMNS)
    if len(final_rows) != policy.expected_query_chromosomes:
        raise PublicationError("Final chromosome map does not have exactly 29 rows")
    final_by_source: dict[str, str] = {}
    for row in final_rows:
        if row["status"] != "PASS_AUTO":
            raise PublicationError("Final map contains a non-PASS row")
        source, final = row["query_chromosome"], row["final_chromosome"]
        if source in final_by_source or final in final_by_source.values():
            raise PublicationError("Final chromosome map is not a bijection")
        final_by_source[source] = final
    matrix = read_score_matrix(matrix_path, role="nucleotide_hy4a", kind="nucleotide", policy=policy)
    orientation_rows = read_exact_tsv(orientation_path, ORIENTATION_COLUMNS)
    support_by_cell: dict[tuple[str, str], dict[str, str]] = {}
    for row in orientation_rows:
        key = (row["query_chromosome"], row["reference_chromosome"])
        if key in support_by_cell:
            raise PublicationError("Orientation support contains duplicate cells")
        support_by_cell[key] = row
    if set(final_by_source) != set(matrix.queries):
        raise PublicationError("Final map and nucleotide matrix query scopes differ")
    orientation_by_source: dict[str, str] = {}
    support_by_source: dict[str, tuple[str, int]] = {}
    for source, final in final_by_source.items():
        references = [
            reference
            for reference, canonical in matrix.canonical_by_reference.items()
            if canonical == final
        ]
        if len(references) != 1:
            raise PublicationError(f"{source}: final canonical label is not unique in matrix")
        reference = references[0]
        matrix_row = matrix.cells[(source, reference)].raw
        support = support_by_cell.get((source, reference))
        if support is None:
            raise PublicationError(f"{source}: assigned orientation-support cell is missing")
        orientation = matrix_row["orientation"]
        dominant = support["dominant_orientation"]
        if orientation not in {"+", "-"} or dominant != orientation:
            raise PublicationError(f"{source}: matrix and orientation audit disagree")
        try:
            plus = int(support["plus_matching_bases"])
            minus = int(support["minus_matching_bases"])
            total = int(support["total_matching_bases"])
        except ValueError as error:
            raise PublicationError(f"{source}: non-integer orientation support") from error
        if min(plus, minus, total) < 0 or plus + minus != total:
            raise PublicationError(f"{source}: orientation support arithmetic mismatch")
        fraction = parse_exact_decimal(support["dominant_fraction"], "dominant_fraction")
        expected_fraction = Decimal(max(plus, minus)) / Decimal(total) if total else Decimal(0)
        if abs(fraction - expected_fraction) > Decimal(str(policy.arithmetic_tolerance)):
            raise PublicationError(f"{source}: dominant fraction arithmetic mismatch")
        if (
            support["automatic_orientation_gate"] != "true"
            or fraction < minimum_fraction
            or total < minimum_matching_bases
        ):
            raise PublicationError(f"{source}: orientation support requires manual review")
        orientation_by_source[source] = orientation
        support_by_source[source] = (str(fraction), total)
    return final_by_source, orientation_by_source, support_by_source


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: dict[str, object]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )


def run(args: argparse.Namespace) -> Path:
    policy = load_homology_policy(args.parameters)
    with args.parameters.open("rb") as handle:
        configuration = tomllib.load(handle)
    orientation_policy = configuration.get("chromosome_orientation_harmonization")
    if not isinstance(orientation_policy, dict):
        raise PublicationError("Parameters lack chromosome_orientation_harmonization")
    minimum_fraction = Decimal(str(orientation_policy.get("minimum_dominant_orientation_fraction")))
    minimum_matching_bases = orientation_policy.get("minimum_oriented_matching_bases")
    if minimum_fraction != Decimal("0.80") or minimum_matching_bases != 1_000_000:
        raise PublicationError("Orientation thresholds differ from the frozen 0.80/1 Mb policy")
    if any(
        orientation_policy.get(key) is not True
        for key in (
            "transform_gff_coordinates_and_strands",
            "require_cds_and_protein_sequence_identity_after_transform",
            "retain_original_publisher_scope_assets",
        )
    ):
        raise PublicationError("Required orientation publication gates are disabled")

    target_registry = read_target_registry(capture_snapshot(args.target_asset_registry))
    target_key = (args.assembly_unit_id, args.target_scope_id)
    if target_key not in target_registry.assets:
        raise PublicationError("Target registry lacks the invoked unit/scope")
    genome_binding = stable_binding(args.source_genome)
    gff_binding = stable_binding(args.source_gff)
    cds_binding = stable_binding(args.source_cds)
    protein_binding = stable_binding(args.source_protein)
    for role, observed in (("genome", genome_binding), ("gff", gff_binding), ("protein", protein_binding)):
        require_binding(target_registry.assets[target_key][role].as_dict(), observed, f"target {role}")
    primary_manifest_binding = validate_primary_annotation_bundle(
        args.primary_annotation_dir,
        gff_binding=gff_binding,
        cds_binding=cds_binding,
        protein_binding=protein_binding,
    )

    matrix_binding = stable_binding(args.nucleotide_hy4a)
    matrix_provenance_binding = stable_binding(args.nucleotide_hy4a_provenance)
    orientation_binding = stable_binding(args.orientation_support)
    matrix_audit_binding = stable_binding(args.nucleotide_input_audit)
    matrix_audit = strict_json(args.nucleotide_input_audit)
    if matrix_audit.get("status") != "PASS" or matrix_audit.get("matrix_role") != "nucleotide_hy4a":
        raise PublicationError("Nucleotide input audit is not the requested PASS role")
    require_binding(matrix_audit.get("matrix"), matrix_binding, "nucleotide input-audit matrix")
    require_binding(
        matrix_audit.get("orientation_support"), orientation_binding, "orientation support audit"
    )
    final_map_path, final_map_binding = validate_assignment_bundle(
        args.assignment_dir,
        assembly_unit_id=args.assembly_unit_id,
        target_scope_id=args.target_scope_id,
        trusted_commit=args.trusted_repository_commit,
        nucleotide_matrix_binding=matrix_binding,
        nucleotide_provenance_binding=matrix_provenance_binding,
    )
    final_by_source, orientation_by_source, support_by_source = load_orientation_actions(
        final_map_path=final_map_path,
        matrix_path=args.nucleotide_hy4a,
        orientation_path=args.orientation_support,
        policy=policy,
        minimum_fraction=minimum_fraction,
        minimum_matching_bases=minimum_matching_bases,
    )
    source_genome = read_fasta(args.source_genome)
    actions = build_actions(source_genome, final_by_source, orientation_by_source)
    transformed_genome = transform_genome(source_genome, actions)

    output = args.output_dir.expanduser().resolve()
    if output.exists():
        raise PublicationError(f"Refusing to overwrite output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        genome_out = staging / f"{args.assembly_unit_id}.harmonized.genome.fa"
        gff_out = staging / f"{args.assembly_unit_id}.harmonized.primary.gff3"
        cds_out = staging / f"{args.assembly_unit_id}.harmonized.cds.fa"
        protein_out = staging / f"{args.assembly_unit_id}.harmonized.protein.faa"
        action_out = staging / f"{args.assembly_unit_id}.chromosome_actions.tsv"
        write_fasta(genome_out, transformed_genome)
        gff_audit = transform_gff(args.source_gff, gff_out, actions)
        shutil.copyfile(args.source_cds, cds_out)
        shutil.copyfile(args.source_protein, protein_out)
        reread_genome = read_fasta(genome_out)
        validate_sequence_closure(
            source_genome=source_genome, transformed_genome=reread_genome, actions=actions
        )
        closure = validate_cds_protein_closure(
            source_genome=source_genome,
            source_gff=args.source_gff,
            transformed_genome=reread_genome,
            transformed_gff=gff_out,
            expected_cds=read_fasta(args.source_cds),
            expected_proteins=read_fasta(args.source_protein),
        )
        rows = [
            {
                "source_chromosome": source,
                "final_chromosome": action.final_chromosome,
                "orientation": action.orientation,
                "source_length_bp": action.source_length,
                "dominant_fraction": support_by_source[source][0],
                "oriented_matching_bases": support_by_source[source][1],
                "status": "PASS_AUTO",
            }
            for source, action in sorted(actions.items())
        ]
        write_tsv(action_out, ACTION_COLUMNS, rows)
        output_bindings = {
            path.name: stable_binding(path)
            for path in (genome_out, gff_out, cds_out, protein_out, action_out)
        }
        validation = {
            "schema_version": 1,
            "workflow": "chromosome_bundle_harmonization",
            "workflow_version": "1.0.0",
            "status": "PASS",
            "publication_gate": "PASS",
            "assembly_unit_id": args.assembly_unit_id,
            "target_scope_id": args.target_scope_id,
            "trusted_repository_commit": args.trusted_repository_commit,
            "counts": {
                "chromosomes": len(actions),
                "positive_orientation": sum(x.orientation == "+" for x in actions.values()),
                "negative_orientation": sum(x.orientation == "-" for x in actions.values()),
                **gff_audit,
                **closure,
            },
            "checks": {
                "assignment_bundle_pass": True,
                "target_registry_reconciled": True,
                "strict_primary_annotation_reconciled": True,
                "orientation_support_reconciled": True,
                "complete_chromosome_bijection": True,
                "genome_sequence_transform_exact": True,
                "gff_coordinate_and_strand_transform_exact": True,
                "cds_sequence_identity": True,
                "protein_sequence_identity": True,
            },
        }
        validation_path = staging / "validation.json"
        write_json(validation_path, validation)
        output_bindings[validation_path.name] = stable_binding(validation_path)
        manifest = {
            "schema_version": 1,
            "workflow": "chromosome_bundle_harmonization",
            "workflow_version": "1.0.0",
            "status": "PASS",
            "publication_gate": "PASS",
            "assembly_unit_id": args.assembly_unit_id,
            "target_scope_id": args.target_scope_id,
            "trusted_repository_commit": args.trusted_repository_commit,
            "inputs": {
                "target_asset_registry": target_registry.snapshot.public_binding(),
                "source_genome": genome_binding,
                "source_gff": gff_binding,
                "source_cds": cds_binding,
                "source_protein": protein_binding,
                "primary_annotation_manifest": primary_manifest_binding,
                "assignment_final_map": final_map_binding,
                "nucleotide_hy4a": matrix_binding,
                "nucleotide_hy4a_provenance": matrix_provenance_binding,
                "nucleotide_input_audit": matrix_audit_binding,
                "orientation_support": orientation_binding,
            },
            "policy": {
                "minimum_dominant_orientation_fraction": str(minimum_fraction),
                "minimum_oriented_matching_bases": minimum_matching_bases,
                "negative_action": "reverse_complement",
                "gff_coordinate_formula": "L-end+1,L-start+1",
                "strand_action": "swap_plus_minus",
            },
            "outputs": output_bindings,
        }
        manifest_path = staging / "run_manifest.json"
        write_json(manifest_path, manifest)
        output_bindings[manifest_path.name] = stable_binding(manifest_path)
        checksum_rows = [
            {"file": name, "bytes": item["bytes"], "sha256": item["sha256"]}
            for name, item in sorted(output_bindings.items())
        ]
        write_tsv(staging / "checksums.tsv", ("file", "bytes", "sha256"), checksum_rows)
        os.rename(staging, output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return output


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assembly-unit-id", required=True)
    p.add_argument("--target-scope-id", required=True)
    p.add_argument("--trusted-repository-commit", required=True)
    p.add_argument("--assignment-dir", required=True, type=Path)
    p.add_argument("--nucleotide-hy4a", required=True, type=Path)
    p.add_argument("--nucleotide-hy4a-provenance", required=True, type=Path)
    p.add_argument("--nucleotide-input-audit", required=True, type=Path)
    p.add_argument("--orientation-support", required=True, type=Path)
    p.add_argument("--source-genome", required=True, type=Path)
    p.add_argument("--source-gff", required=True, type=Path)
    p.add_argument("--source-cds", required=True, type=Path)
    p.add_argument("--source-protein", required=True, type=Path)
    p.add_argument("--primary-annotation-dir", required=True, type=Path)
    p.add_argument("--parameters", required=True, type=Path)
    p.add_argument("--target-asset-registry", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    return p


def main() -> int:
    try:
        output = run(parser().parse_args())
        print(f"PASS\t{output}")
        return 0
    except (PublicationError, HarmonizationError, OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
