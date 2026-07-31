#!/usr/bin/env python3
"""Run non-publishing chromosome-assignment diagnostics before the final commit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from geneloss_repro.chromosome_assignment import (  # noqa: E402
    _generation_parameters,
    _matrix_pair_compatible,
    load_homology_policy,
    read_score_matrix,
    solve_score_matrix,
)
from geneloss_repro.chromosome_provenance import (  # noqa: E402
    capture_snapshot,
    read_reference_registries,
    read_target_registry,
    validate_matrix_provenance,
)


class DiagnosticError(RuntimeError):
    pass


ROLES = {
    "nucleotide_hy4a": ("nucleotide", "HY4A", "primary_coordinate_reference"),
    "jcvi_hy4a": ("jcvi", "HY4A", "primary_coordinate_reference"),
    "nucleotide_hy4p": ("nucleotide", "HY4P", "independent_confirmation_reference"),
    "jcvi_hy4p": ("jcvi", "HY4P", "independent_confirmation_reference"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink() or resolved.stat().st_size == 0:
        raise DiagnosticError(f"missing, empty, or symlink file: {resolved}")
    return resolved


def under(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute() or ".." in Path(relative).parts:
        raise DiagnosticError(f"unsafe data-root-relative path: {relative!r}")
    resolved = (root / relative).resolve()
    if not resolved.is_relative_to(root):
        raise DiagnosticError(f"path escapes data root: {relative!r}")
    return resolved


def binding(path: Path) -> dict[str, object]:
    source = regular(path)
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def write_unit(
    *,
    row: dict[str, str],
    root: Path,
    output: Path,
    policy,
    parameters_snapshot,
    reference_asset_snapshot,
    reference_map_snapshot,
    references,
) -> dict[str, object]:
    unit = row["assembly_unit_id"]
    scope = row["target_scope_id"]
    matrix_root = under(root, row["matrix_root"])
    target_registry_snapshot = capture_snapshot(under(root, row["target_asset_registry"]))
    target_registry = read_target_registry(target_registry_snapshot)
    if (unit, scope) not in target_registry.assets:
        raise DiagnosticError(f"target registry lacks {unit}/{scope}")

    matrices = {}
    inputs = {}
    for role, (kind, slot, reference_role) in ROLES.items():
        reference_id = (
            policy.coordinate_reference if slot == "HY4A" else policy.confirmation_reference
        )
        matrix_path = regular(matrix_root / role / f"{role}.tsv")
        sidecar_path = regular(matrix_root / role / f"{role}.provenance.json")
        matrix_snapshot = capture_snapshot(matrix_path)
        sidecar_snapshot = capture_snapshot(sidecar_path)
        matrix = read_score_matrix(matrix_snapshot, role=role, kind=kind, policy=policy)
        if dict(matrix.canonical_by_reference) != dict(references.chromosome_maps[reference_id]):
            raise DiagnosticError(f"{unit}/{role}: canonical map conflicts with frozen truth")
        provenance = validate_matrix_provenance(
            sidecar=sidecar_snapshot,
            matrix=matrix_snapshot,
            assembly_unit_id=unit,
            target_scope_id=scope,
            matrix_role=role,
            matrix_kind=kind,
            reference_slot=slot,
            reference_id=reference_id,
            reference_role=reference_role,
            reference_map_id=references.map_ids[reference_id],
            expected_generation_parameters=_generation_parameters(policy, kind),
            target_registry=target_registry,
            reference_registries=references,
        )
        matrices[role] = matrix
        inputs[role] = {
            "matrix": matrix_snapshot.public_binding(),
            "provenance": sidecar_snapshot.public_binding(),
            "upstream_validation": provenance.upstream_report_snapshot.public_binding(),
        }

    _matrix_pair_compatible(matrices["nucleotide_hy4a"], matrices["jcvi_hy4a"], "HY4A")
    _matrix_pair_compatible(matrices["nucleotide_hy4p"], matrices["jcvi_hy4p"], "HY4P")
    query_sets = {matrix.queries for matrix in matrices.values()}
    if len(query_sets) != 1:
        raise DiagnosticError(f"{unit}: four matrices disagree on query IDs")
    assignments = {role: solve_score_matrix(matrix, policy) for role, matrix in matrices.items()}
    queries = matrices["nucleotide_hy4a"].queries

    rows = []
    label_rows = []
    agreement_a = agreement_p = cross_reference = all_gates = 0
    for query in queries:
        na, ja = assignments["nucleotide_hy4a"][query], assignments["jcvi_hy4a"][query]
        np, jp = assignments["nucleotide_hy4p"][query], assignments["jcvi_hy4p"][query]
        agree_a = na.canonical == ja.canonical
        agree_p = np.canonical == jp.canonical
        cross = agree_a and agree_p and na.canonical == np.canonical
        gates = all(item.matrix_gate for item in (na, ja, np, jp))
        agreement_a += int(agree_a)
        agreement_p += int(agree_p)
        cross_reference += int(cross)
        all_gates += int(gates)
        failures = sorted({reason for item in (na, ja, np, jp) for reason in item.failure_reasons})
        if not agree_a or not agree_p:
            failures.append("CONFLICT_NUCLEOTIDE_JCVI")
        elif not cross:
            failures.append("HY4A_HY4P_DISAGREEMENT")
        rows.append(
            {
                "query_chromosome": query,
                "diagnostic_candidate": na.canonical if cross else "",
                "nucleotide_hy4a": na.canonical,
                "jcvi_hy4a": ja.canonical,
                "nucleotide_hy4p": np.canonical,
                "jcvi_hy4p": jp.canonical,
                "nucleotide_hy4a_reciprocal_coverage": "" if na.reciprocal_coverage is None else f"{na.reciprocal_coverage:.12g}",
                "jcvi_hy4a_score": f"{ja.score:.12g}",
                "orientation_hy4a": na.orientation,
                "all_four_matrix_gates": "true" if gates else "false",
                "diagnostic_status": "PASS_DIAGNOSTIC" if cross and gates else "REVIEW_REQUIRED",
                "failure_reasons": ";".join(sorted(set(failures))),
            }
        )
        confidence = "HIGH" if cross and gates else ("SUPPORTED" if cross else "LOW")
        label_rows.append(
            {
                "query_chromosome": query,
                "final_chromosome": na.canonical,
                "coordinate_reference": policy.coordinate_reference,
                "assignment_method": "global_one_to_one_maximum_nucleotide_similarity",
                "assigned_score": f"{na.score:.12g}",
                "reciprocal_coverage": (
                    "" if na.reciprocal_coverage is None else f"{na.reciprocal_coverage:.12g}"
                ),
                "orientation_to_hy4a": na.orientation,
                "hy4p_and_jcvi_agree": "true" if cross else "false",
                "strict_homology_gates_pass": "true" if gates else "false",
                "confidence_flag": confidence,
            }
        )

    complete = len(queries) == policy.expected_query_chromosomes
    candidates = [item["diagnostic_candidate"] for item in rows]
    bijective = complete and all(candidates) and len(set(candidates)) == len(candidates)
    diagnostic_pass = (
        agreement_a == len(queries)
        and agreement_p == len(queries)
        and cross_reference == len(queries)
        and all_gates == len(queries)
        and bijective
    )
    label_candidates = [item["final_chromosome"] for item in label_rows]
    naming_pass = (
        complete
        and len(label_candidates) == len(set(label_candidates))
        and not any(not item for item in label_candidates)
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        table = staging / "prefinal_assignment_diagnostic.tsv"
        with table.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        label_table = staging / "similarity_label_map.tsv"
        with label_table.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(label_rows[0]),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(label_rows)
        report = {
            "schema_version": 1,
            "workflow": "prefinal_nonpublishing_chromosome_assignment_diagnostic",
            "status": "PASS_DIAGNOSTIC" if diagnostic_pass else "REVIEW_REQUIRED",
            "chromosome_naming_status": "PASS_LABELS" if naming_pass else "FAIL_LABELS",
            "chromosome_naming_policy": (
                "HY4A global one-to-one maximum nucleotide similarity; absolute support is QC only"
            ),
            "publication_allowed": False,
            "publication_block": "final assignment must be regenerated after the one reviewed final repository commit",
            "assembly_unit_id": unit,
            "target_scope_id": scope,
            "query_count": len(queries),
            "nucleotide_jcvi_agreement_hy4a": agreement_a,
            "nucleotide_jcvi_agreement_hy4p": agreement_p,
            "hy4a_hy4p_label_agreement": cross_reference,
            "rows_passing_all_four_matrix_gates": all_gates,
            "diagnostic_candidates_form_bijection": bool(bijective),
            "parameters": parameters_snapshot.public_binding(),
            "target_asset_registry": target_registry_snapshot.public_binding(),
            "reference_asset_registry": reference_asset_snapshot.public_binding(),
            "reference_chromosome_map_registry": reference_map_snapshot.public_binding(),
            "matrix_inputs": inputs,
            "diagnostic_table": binding(table),
            "similarity_label_map": binding(label_table),
        }
        report_path = staging / "diagnostic.json"
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
        (staging / "checksums.tsv").write_text(
            "file\tsha256\n" + "".join(
                f"{path.name}\t{sha256(path)}\n" for path in sorted(staging.iterdir()) if path.is_file()
            ),
            encoding="utf-8",
        )
        os.rename(staging, output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--parameters", required=True, type=Path)
    parser.add_argument("--reference-assets", required=True, type=Path)
    parser.add_argument("--reference-maps", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    batch_staging: Path | None = None
    try:
        root = args.data_root.expanduser().resolve()
        output_root = args.output_root.expanduser().resolve()
        if output_root.exists():
            raise DiagnosticError(f"refusing to overwrite output root: {output_root}")
        manifest = regular(args.manifest)
        with manifest.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        required = {"assembly_unit_id", "target_scope_id", "target_asset_registry", "matrix_root"}
        if not rows or not required.issubset(rows[0]):
            raise DiagnosticError("manifest is empty or lacks required columns")
        if len({row["assembly_unit_id"] for row in rows}) != len(rows):
            raise DiagnosticError("duplicate units in manifest")

        parameters_snapshot = capture_snapshot(args.parameters)
        reference_asset_snapshot = capture_snapshot(args.reference_assets)
        reference_map_snapshot = capture_snapshot(args.reference_maps)
        policy = load_homology_policy(parameters_snapshot)
        references = read_reference_registries(
            reference_asset_snapshot,
            reference_map_snapshot,
            coordinate_reference=policy.coordinate_reference,
            confirmation_reference=policy.confirmation_reference,
            expected_chromosomes=policy.expected_reference_chromosomes,
        )
        output_root.parent.mkdir(parents=True, exist_ok=True)
        batch_staging = Path(
            tempfile.mkdtemp(prefix=f".{output_root.name}.staging.", dir=output_root.parent)
        )
        summaries = []
        for row in rows:
            summaries.append(
                write_unit(
                    row=row,
                    root=root,
                    output=batch_staging / row["assembly_unit_id"],
                    policy=policy,
                    parameters_snapshot=parameters_snapshot,
                    reference_asset_snapshot=reference_asset_snapshot,
                    reference_map_snapshot=reference_map_snapshot,
                    references=references,
                )
            )
        overall = {
            "schema_version": 1,
            "workflow": "prefinal_nonpublishing_chromosome_assignment_diagnostic_batch",
            "status": "DIAGNOSTIC_COMPLETE",
            "publication_allowed": False,
            "manifest": binding(manifest),
            "unit_count": len(summaries),
            "pass_diagnostic_count": sum(item["status"] == "PASS_DIAGNOSTIC" for item in summaries),
            "review_required_count": sum(item["status"] == "REVIEW_REQUIRED" for item in summaries),
            "pass_label_count": sum(
                item["chromosome_naming_status"] == "PASS_LABELS" for item in summaries
            ),
            "units": [
                {
                    "assembly_unit_id": item["assembly_unit_id"],
                    "status": item["status"],
                    "agreement": item["hy4a_hy4p_label_agreement"],
                    "all_gates": item["rows_passing_all_four_matrix_gates"],
                    "chromosome_naming_status": item["chromosome_naming_status"],
                }
                for item in summaries
            ],
        }
        (batch_staging / "batch_diagnostic.json").write_text(
            json.dumps(overall, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        os.rename(batch_staging, output_root)
        batch_staging = None
        print(f"DIAGNOSTIC_COMPLETE\t{len(summaries)} units")
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, DiagnosticError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if batch_staging is not None:
            shutil.rmtree(batch_staging, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
