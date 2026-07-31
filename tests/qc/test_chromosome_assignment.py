"""Synthetic regression tests for fail-closed chromosome assignment."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

import geneloss_repro.chromosome_assignment as chromosome_assignment_module
from geneloss_repro.chromosome_assignment import (
    JCVI_COLUMNS,
    MatrixCell,
    NUCLEOTIDE_COLUMNS,
    ScoreMatrix,
    ChromosomeAssignmentError,
    _generation_parameters,
    _hungarian_maximize,
    assign_chromosome_homology,
    load_homology_policy,
    solve_score_matrix,
)
from geneloss_repro.chromosome_provenance import (
    JCVI_UPSTREAM_CHECKS,
    NUCLEOTIDE_UPSTREAM_CHECKS,
    PROVENANCE_SCHEMA_VERSION,
    ChromosomeProvenanceError,
    capture_snapshot,
    parse_exact_int,
    read_reference_registries,
    read_target_registry,
    verify_snapshot,
)


ROOT = Path(__file__).resolve().parents[2]
PARAMETERS = ROOT / "config" / "analysis_parameters.toml"
REFERENCE_ASSET_REGISTRY = ROOT / "config" / "chromosome_coordinate_references.tsv"
REFERENCE_MAP_REGISTRY = ROOT / "config" / "chromosome_reference_maps.tsv"
COUNT = 29
ASSEMBLY_UNIT = "act_test_hap1"
TARGET_SCOPE = "act_test_hap1_chromosomes_v1"
TRUSTED_COMMIT = "0123456789abcdef0123456789abcdef01234567"


def harmonic(first: float, second: float) -> float:
    return 0.0 if first <= 0 or second <= 0 else 2 * first * second / (first + second)


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def target_indices() -> list[int]:
    return [(index + 7) % COUNT for index in range(COUNT)]


def make_nucleotide(
    path: Path,
    suffix: str,
    targets: list[int],
    *,
    low_coverage_query: int | None = None,
    tie: tuple[int, int] | None = None,
) -> None:
    rows: list[dict[str, str]] = []
    for query_index in range(COUNT):
        for reference_index in range(COUNT):
            selected = reference_index == targets[query_index] or tie == (
                query_index,
                reference_index,
            )
            coverage = 0.9 if selected else 0.02
            if selected and low_coverage_query == query_index:
                coverage = 0.049
            divergence = 0.01 if selected else 0.10
            score = harmonic(coverage, coverage) * (1.0 - divergence)
            covered = round(coverage * 1000)
            rows.append(
                {
                    "query_chromosome": f"PubChr{query_index + 1:02d}",
                    "reference_chromosome": f"Chr{reference_index + 1:02d}{suffix}",
                    "canonical_chromosome": f"Chr{reference_index + 1:02d}",
                    "score": f"{score:.12g}",
                    "query_covered_bp": str(covered),
                    "query_length_bp": "1000",
                    "query_coverage": f"{coverage:.12g}",
                    "reference_covered_bp": str(covered),
                    "reference_length_bp": "1000",
                    "reference_coverage": f"{coverage:.12g}",
                    "reciprocal_coverage": f"{coverage:.12g}",
                    "matching_bases": str(max(1, covered - 1)),
                    "weighted_divergence": f"{divergence:.12g}",
                    "orientation": "+",
                }
            )
    write_tsv(path, NUCLEOTIDE_COLUMNS, rows)


def make_jcvi(
    path: Path,
    suffix: str,
    targets: list[int],
    *,
    tie: tuple[int, int] | None = None,
) -> None:
    rows: list[dict[str, str]] = []
    for query_index in range(COUNT):
        for reference_index in range(COUNT):
            selected = reference_index == targets[query_index] or tie == (
                query_index,
                reference_index,
            )
            anchored = 90 if selected else 2
            coverage = anchored / 100
            score = harmonic(coverage, coverage)
            rows.append(
                {
                    "query_chromosome": f"PubChr{query_index + 1:02d}",
                    "reference_chromosome": f"Chr{reference_index + 1:02d}{suffix}",
                    "canonical_chromosome": f"Chr{reference_index + 1:02d}",
                    "score": f"{score:.12g}",
                    "query_anchored_genes": str(anchored),
                    "query_eligible_genes": "100",
                    "query_gene_coverage": f"{coverage:.12g}",
                    "reference_anchored_genes": str(anchored),
                    "reference_eligible_genes": "100",
                    "reference_gene_coverage": f"{coverage:.12g}",
                    "unique_anchor_pairs": str(anchored),
                }
            )
    write_tsv(path, JCVI_COLUMNS, rows)


def make_inputs(
    root: Path,
    *,
    jcvi_a_targets: list[int] | None = None,
    nucleotide_p_targets: list[int] | None = None,
    jcvi_p_targets: list[int] | None = None,
    low_coverage_query: int | None = None,
    nucleotide_a_tie: tuple[int, int] | None = None,
) -> dict[str, Path]:
    targets = target_indices()
    paths = {
        "nucleotide_hy4a": root / "nucleotide_hy4a.tsv",
        "jcvi_hy4a": root / "jcvi_hy4a.tsv",
        "nucleotide_hy4p": root / "nucleotide_hy4p.tsv",
        "jcvi_hy4p": root / "jcvi_hy4p.tsv",
    }
    make_nucleotide(
        paths["nucleotide_hy4a"],
        "A",
        targets,
        low_coverage_query=low_coverage_query,
        tie=nucleotide_a_tie,
    )
    make_jcvi(paths["jcvi_hy4a"], "A", jcvi_a_targets or targets)
    make_nucleotide(
        paths["nucleotide_hy4p"],
        "P",
        nucleotide_p_targets or targets,
        low_coverage_query=low_coverage_query,
    )
    make_jcvi(paths["jcvi_hy4p"], "P", jcvi_p_targets or targets)
    target_registry = root / "target_asset_registry.tsv"
    target_rows = []
    for role, file_name in (
        ("genome", "act_test_hap1.chromosomes.fasta"),
        ("gff", "act_test_hap1.chromosomes.gff3"),
        ("protein", "act_test_hap1.primary.protein.fasta"),
    ):
        payload = f"{ASSEMBLY_UNIT}:{TARGET_SCOPE}:{role}".encode("utf-8")
        target_rows.append(
            {
                "assembly_unit_id": ASSEMBLY_UNIT,
                "target_scope_id": TARGET_SCOPE,
                "asset_role": role,
                "file_name": file_name,
                "bytes": str(len(payload)),
                "sha256": hashlib.sha256(payload).hexdigest(),
                "status": "verified",
            }
        )
    write_tsv(
        target_registry,
        (
            "assembly_unit_id",
            "target_scope_id",
            "asset_role",
            "file_name",
            "bytes",
            "sha256",
            "status",
        ),
        target_rows,
    )
    paths["target_asset_registry"] = target_registry
    refresh_provenance(root, paths)
    return paths


def binding(path: Path) -> dict[str, str | int]:
    payload = path.read_bytes()
    return {
        "basename": path.name,
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def refresh_provenance(root: Path, paths: dict[str, Path]) -> None:
    """Create exact role-specific PASS reports and checksum-bound sidecars."""

    policy = load_homology_policy(PARAMETERS)
    target_registry = read_target_registry(capture_snapshot(paths["target_asset_registry"]))
    references = read_reference_registries(
        capture_snapshot(REFERENCE_ASSET_REGISTRY),
        capture_snapshot(REFERENCE_MAP_REGISTRY),
        coordinate_reference=policy.coordinate_reference,
        confirmation_reference=policy.confirmation_reference,
        expected_chromosomes=COUNT,
    )
    role_specs = {
        "nucleotide_hy4a": (
            "nucleotide",
            "HY4A",
            policy.coordinate_reference,
            "primary_coordinate_reference",
        ),
        "jcvi_hy4a": (
            "jcvi",
            "HY4A",
            policy.coordinate_reference,
            "primary_coordinate_reference",
        ),
        "nucleotide_hy4p": (
            "nucleotide",
            "HY4P",
            policy.confirmation_reference,
            "independent_confirmation_reference",
        ),
        "jcvi_hy4p": (
            "jcvi",
            "HY4P",
            policy.confirmation_reference,
            "independent_confirmation_reference",
        ),
    }
    target_assets = {
        role: asset.as_dict()
        for role, asset in target_registry.assets[(ASSEMBLY_UNIT, TARGET_SCOPE)].items()
    }
    for matrix_role, (kind, slot, reference_id, reference_role) in role_specs.items():
        matrix = paths[matrix_role]
        generation = _generation_parameters(policy, kind)
        report_path = root / f"{matrix_role}.upstream_validation.json"
        report = {
            "schema_version": "1.0.0",
            "workflow": (
                "chromosome_nucleotide_matrix_builder"
                if kind == "nucleotide"
                else "chromosome_jcvi_matrix_builder"
            ),
            "workflow_version": "1.0.0",
            "status": "PASS",
            "assembly_unit_id": ASSEMBLY_UNIT,
            "target_scope_id": TARGET_SCOPE,
            "matrix_role": matrix_role,
            "matrix_kind": kind,
            "reference_slot": slot,
            "reference_id": reference_id,
            "reference_role": reference_role,
            "reference_map_id": references.map_ids[reference_id],
            "matrix_sha256": binding(matrix)["sha256"],
            "target_asset_registry_sha256": target_registry.snapshot.sha256,
            "reference_asset_registry_sha256": references.asset_snapshot.sha256,
            "reference_chromosome_map_registry_sha256": references.map_snapshot.sha256,
            "generation_parameters": generation,
            "checks": {
                key: True
                for key in sorted(
                    NUCLEOTIDE_UPSTREAM_CHECKS
                    if kind == "nucleotide"
                    else JCVI_UPSTREAM_CHECKS
                )
            },
        }
        report_path.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        sidecar_path = root / f"{matrix_role}.provenance.json"
        sidecar = {
            "schema_version": PROVENANCE_SCHEMA_VERSION,
            "status": "PASS",
            "assembly_unit_id": ASSEMBLY_UNIT,
            "target_scope_id": TARGET_SCOPE,
            "matrix_role": matrix_role,
            "matrix_kind": kind,
            "reference_slot": slot,
            "reference_id": reference_id,
            "reference_role": reference_role,
            "reference_map_id": references.map_ids[reference_id],
            "matrix": binding(matrix),
            "target_asset_registry": target_registry.snapshot.public_binding(),
            "target_assets": target_assets,
            "reference_asset_registry": references.asset_snapshot.public_binding(),
            "reference_assets": {
                role: asset.as_dict()
                for role, asset in references.assets[reference_id].items()
            },
            "reference_chromosome_map_registry": references.map_snapshot.public_binding(),
            "generation_parameters": generation,
            "upstream_validation_report": binding(report_path),
        }
        sidecar_path.write_text(
            json.dumps(sidecar, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        paths[f"{matrix_role}_provenance"] = sidecar_path


def run_assignment(root: Path, paths: dict[str, Path], **overrides):
    arguments = dict(paths)
    assembly_unit_id = overrides.pop("assembly_unit_id", ASSEMBLY_UNIT)
    target_scope_id = overrides.pop("target_scope_id", TARGET_SCOPE)
    trusted_repository_commit = overrides.pop(
        "trusted_repository_commit", TRUSTED_COMMIT
    )
    output_dir = overrides.pop("output_dir", root / "assignment")
    parameters = overrides.pop("parameters", PARAMETERS)
    reference_asset_registry = overrides.pop(
        "reference_asset_registry", REFERENCE_ASSET_REGISTRY
    )
    reference_map_registry = overrides.pop(
        "reference_chromosome_map_registry", REFERENCE_MAP_REGISTRY
    )
    arguments.update(overrides)
    return assign_chromosome_homology(
        **arguments,
        parameters=parameters,
        reference_asset_registry=reference_asset_registry,
        reference_chromosome_map_registry=reference_map_registry,
        assembly_unit_id=assembly_unit_id,
        target_scope_id=target_scope_id,
        trusted_repository_commit=trusted_repository_commit,
        output_dir=output_dir,
    )


def rewrite_json(path: Path, mutate) -> None:
    document = json.loads(path.read_text(encoding="utf-8"))
    mutate(document)
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def rebind_report(sidecar_path: Path, report_path: Path) -> None:
    rewrite_json(
        sidecar_path,
        lambda document: document.__setitem__(
            "upstream_validation_report", binding(report_path)
        ),
    )


class ChromosomeAssignmentTests(unittest.TestCase):
    def test_hungarian_finds_global_optimum_that_row_greedy_misses(self) -> None:
        scores = ((10.0, 9.0, 0.0), (9.0, 0.0, 0.0), (0.0, 8.0, 7.0))
        self.assertEqual(_hungarian_maximize(scores), (1, 0, 2))

    def test_cli_help_works_without_pythonpath(self) -> None:
        script = ROOT / "scripts" / "qc" / "assign_chromosome_homology.py"
        environment = dict(os.environ)
        environment.pop("PYTHONPATH", None)
        with tempfile.TemporaryDirectory() as temporary:
            completed = subprocess.run(
                [sys.executable, str(script), "--help"],
                cwd=temporary,
                env=environment,
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--nucleotide-hy4a", completed.stdout)

    def test_exact_four_matrix_agreement_publishes_complete_path_free_map(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            result = run_assignment(root, paths)
            self.assertEqual(result.status, "PASS_AUTO")
            self.assertEqual(result.publication_gate, "PASS")
            self.assertEqual(result.final_map_row_count, 29)
            output = root / "assignment"
            final_rows = read_tsv(output / "act_test_hap1.final_chromosome_map.tsv")
            self.assertEqual(len(final_rows), 29)
            self.assertEqual(len({row["final_chromosome"] for row in final_rows}), 29)
            first = {row["query_chromosome"]: row for row in final_rows}["PubChr01"]
            self.assertEqual(first["final_chromosome"], "Chr08")
            self.assertEqual(first["status"], "PASS_AUTO")

            validation = json.loads((output / "validation.json").read_text(encoding="utf-8"))
            self.assertEqual(validation["publication_gate"], "PASS")
            self.assertTrue(validation["checks"]["final_map_is_complete_bijection"])
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["matrix_schema_version"], "1.0.0")
            self.assertEqual(
                manifest["inputs"]["matrix_roles"]["nucleotide_hy4a"]["matrix"]["basename"],
                "nucleotide_hy4a.tsv",
            )
            self.assertEqual(manifest["target_scope_id"], TARGET_SCOPE)
            self.assertEqual(manifest["trusted_repository_commit"], TRUSTED_COMMIT)
            self.assertEqual(manifest["policy"]["minimum_unique_anchor_pairs"], 30)
            self.assertEqual(manifest["policy"]["minimap2_version"], "2.28-r1209")
            self.assertEqual(
                manifest["policy"]["minimap2_command_template"],
                [
                    "minimap2",
                    "-x",
                    "asm5",
                    "--secondary=no",
                    "-c",
                    "--cs=long",
                    "{reference_fasta}",
                    "{query_fasta}",
                ],
            )

            checksum_rows = read_tsv(output / "checksums.tsv")
            self.assertNotIn("checksums.tsv", {row["file"] for row in checksum_rows})
            for row in checksum_rows:
                payload = (output / row["file"]).read_bytes()
                self.assertEqual(len(payload), int(row["bytes"]))
                self.assertEqual(hashlib.sha256(payload).hexdigest(), row["sha256"])
            all_text = "".join(
                path.read_text(encoding="utf-8") for path in output.iterdir() if path.is_file()
            )
            self.assertNotIn(str(root), all_text)

    def test_nucleotide_jcvi_conflict_publishes_failure_audit_without_final_map(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            altered = target_indices()
            altered[0], altered[1] = altered[1], altered[0]
            paths = make_inputs(root, jcvi_a_targets=altered)
            result = run_assignment(root, paths)
            self.assertEqual(result.status, "CONFLICT_NUCLEOTIDE_JCVI")
            self.assertEqual(result.publication_gate, "FAIL")
            self.assertFalse(
                (
                    root
                    / "assignment"
                    / "act_test_hap1.final_chromosome_map.tsv"
                ).exists()
            )
            summary = read_tsv(
                root / "assignment" / "act_test_hap1.chromosome_assignment.summary.tsv"
            )[0]
            self.assertEqual(summary["nucleotide_jcvi_agreement_count_hy4a"], "27")
            self.assertIn("CONFLICT_NUCLEOTIDE_JCVI", summary["failure_states"])
            combined = read_tsv(
                root / "assignment" / "act_test_hap1.chromosome_assignment.tsv"
            )
            self.assertIn(
                "NOT_PUBLISHED_UNIT_FAILURE", {row["status"] for row in combined}
            )

    def test_cli_returns_nonzero_but_keeps_valid_failure_audit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            altered = target_indices()
            altered[0], altered[1] = altered[1], altered[0]
            paths = make_inputs(root, jcvi_a_targets=altered)
            script = ROOT / "scripts" / "qc" / "assign_chromosome_homology.py"
            output = root / "cli_failure"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--nucleotide-hy4a",
                    str(paths["nucleotide_hy4a"]),
                    "--jcvi-hy4a",
                    str(paths["jcvi_hy4a"]),
                    "--nucleotide-hy4p",
                    str(paths["nucleotide_hy4p"]),
                    "--jcvi-hy4p",
                    str(paths["jcvi_hy4p"]),
                    "--nucleotide-hy4a-provenance",
                    str(paths["nucleotide_hy4a_provenance"]),
                    "--jcvi-hy4a-provenance",
                    str(paths["jcvi_hy4a_provenance"]),
                    "--nucleotide-hy4p-provenance",
                    str(paths["nucleotide_hy4p_provenance"]),
                    "--jcvi-hy4p-provenance",
                    str(paths["jcvi_hy4p_provenance"]),
                    "--parameters",
                    str(PARAMETERS),
                    "--target-asset-registry",
                    str(paths["target_asset_registry"]),
                    "--reference-asset-registry",
                    str(REFERENCE_ASSET_REGISTRY),
                    "--reference-chromosome-map-registry",
                    str(REFERENCE_MAP_REGISTRY),
                    "--assembly-unit-id",
                    ASSEMBLY_UNIT,
                    "--target-scope-id",
                    TARGET_SCOPE,
                    "--trusted-repository-commit",
                    TRUSTED_COMMIT,
                    "--output-dir",
                    str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 2, completed.stderr)
            self.assertIn("status\tCONFLICT_NUCLEOTIDE_JCVI", completed.stdout)
            self.assertTrue((output / "validation.json").is_file())
            self.assertFalse((output / "act_test_hap1.final_chromosome_map.tsv").exists())

    def test_confirmation_reference_disagreement_is_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            altered = target_indices()
            altered[0], altered[1] = altered[1], altered[0]
            paths = make_inputs(
                root,
                nucleotide_p_targets=altered,
                jcvi_p_targets=altered,
            )
            result = run_assignment(root, paths)
            self.assertEqual(result.status, "HY4A_HY4P_DISAGREEMENT")
            self.assertIn("HY4A_HY4P_DISAGREEMENT", result.failure_states)
            self.assertEqual(result.final_map_row_count, 0)

    def test_low_reciprocal_nucleotide_coverage_stops_entire_unit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = run_assignment(root, make_inputs(root, low_coverage_query=0))
            self.assertEqual(result.status, "LOW_RECIPROCAL_COVERAGE")
            self.assertEqual(result.final_map_row_count, 0)
            evidence = {
                row["query_chromosome"]: row
                for row in read_tsv(root / "assignment" / "nucleotide_hy4a.assignment_evidence.tsv")
            }
            self.assertEqual(evidence["PubChr01"]["coverage_gate"], "false")
            self.assertEqual(evidence["PubChr01"]["assigned_reciprocal_coverage"], "0.049")

    def test_tied_evidence_is_not_resolved_by_hungarian_tie_break(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            targets = target_indices()
            extra = (0, targets[1])
            result = run_assignment(root, make_inputs(root, nucleotide_a_tie=extra))
            self.assertIn("AMBIGUOUS_RECIPROCAL_BEST", result.failure_states)
            self.assertIn("AMBIGUOUS_SEPARATION", result.failure_states)
            self.assertEqual(result.publication_gate, "FAIL")
            self.assertEqual(result.final_map_row_count, 0)

    def test_incomplete_matrix_is_rejected_before_any_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["jcvi_hy4p"])
            write_tsv(paths["jcvi_hy4p"], JCVI_COLUMNS, rows[:-1])
            refresh_provenance(root, paths)
            with self.assertRaisesRegex(ChromosomeAssignmentError, "requires 841 rows"):
                run_assignment(root, paths)
            self.assertFalse((root / "assignment").exists())
            self.assertFalse(list(root.glob(".assignment.staging.*")))

    def test_exact_reference_id_sets_must_match_between_evidence_types(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["jcvi_hy4a"])
            for row in rows:
                if row["reference_chromosome"] == "Chr29A":
                    row["reference_chromosome"] = "renamed_Chr29A"
            write_tsv(paths["jcvi_hy4a"], JCVI_COLUMNS, rows)
            refresh_provenance(root, paths)
            with self.assertRaisesRegex(
                ChromosomeAssignmentError,
                "matrix reference IDs do not equal the frozen HY4A map",
            ):
                run_assignment(root, paths)
            self.assertFalse((root / "assignment").exists())

    def test_arithmetic_mismatch_is_rejected_before_any_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["nucleotide_hy4a"])
            rows[0]["reciprocal_coverage"] = "0.5"
            write_tsv(paths["nucleotide_hy4a"], NUCLEOTIDE_COLUMNS, rows)
            refresh_provenance(root, paths)
            with self.assertRaisesRegex(ChromosomeAssignmentError, "arithmetic mismatch"):
                run_assignment(root, paths)
            self.assertFalse((root / "assignment").exists())

    def test_existing_output_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            output = root / "assignment"
            output.mkdir()
            sentinel = output / "keep.txt"
            sentinel.write_text("keep\n", encoding="utf-8")
            with self.assertRaisesRegex(ChromosomeAssignmentError, "already exists"):
                run_assignment(root, paths)
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep\n")

    def test_copied_or_hardlinked_matrix_role_substitution_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            shutil.copyfile(paths["nucleotide_hy4a"], paths["nucleotide_hy4p"])
            with self.assertRaisesRegex(ChromosomeProvenanceError, "identical bytes"):
                run_assignment(root, paths)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            hardlink = root / "hardlinked_as_hy4p.tsv"
            os.link(paths["nucleotide_hy4a"], hardlink)
            with self.assertRaisesRegex(ChromosomeProvenanceError, "one inode"):
                run_assignment(root, paths, nucleotide_hy4p=hardlink)

    def test_a_p_matrix_and_sidecar_swaps_are_rejected_by_exact_role(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            with self.assertRaisesRegex(
                (ChromosomeAssignmentError, ChromosomeProvenanceError),
                "HY4A|matrix_role|reference",
            ):
                run_assignment(
                    root,
                    paths,
                    nucleotide_hy4a=paths["nucleotide_hy4p"],
                    jcvi_hy4a=paths["jcvi_hy4p"],
                    nucleotide_hy4p=paths["nucleotide_hy4a"],
                    jcvi_hy4p=paths["jcvi_hy4a"],
                    nucleotide_hy4a_provenance=paths["nucleotide_hy4p_provenance"],
                    jcvi_hy4a_provenance=paths["jcvi_hy4p_provenance"],
                    nucleotide_hy4p_provenance=paths["nucleotide_hy4a_provenance"],
                    jcvi_hy4p_provenance=paths["jcvi_hy4a_provenance"],
                )

    def test_wrong_assembly_unit_is_rejected_even_with_same_pubchr_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            with self.assertRaisesRegex(ChromosomeProvenanceError, "assembly_unit_id"):
                run_assignment(root, paths, assembly_unit_id="act_wrong_hap1")

    def test_target_registry_is_one_scope_with_distinct_role_assets(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["target_asset_registry"])
            rows.extend(
                {
                    **row,
                    "assembly_unit_id": "act_other_hap1",
                    "target_scope_id": "act_other_hap1_chromosomes_v1",
                }
                for row in list(rows)
            )
            write_tsv(paths["target_asset_registry"], tuple(rows[0]), rows)
            with self.assertRaisesRegex(
                ChromosomeProvenanceError, "exactly one assembly-unit/target-scope"
            ):
                read_target_registry(capture_snapshot(paths["target_asset_registry"]))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["target_asset_registry"])
            genome_sha256 = next(
                row["sha256"] for row in rows if row["asset_role"] == "genome"
            )
            next(row for row in rows if row["asset_role"] == "protein")[
                "sha256"
            ] = genome_sha256
            write_tsv(paths["target_asset_registry"], tuple(rows[0]), rows)
            with self.assertRaisesRegex(ChromosomeProvenanceError, "copy-identical"):
                read_target_registry(capture_snapshot(paths["target_asset_registry"]))

    def test_all_four_consistent_canonical_swap_conflicts_with_frozen_map(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            for role, columns in (
                ("nucleotide_hy4a", NUCLEOTIDE_COLUMNS),
                ("jcvi_hy4a", JCVI_COLUMNS),
                ("nucleotide_hy4p", NUCLEOTIDE_COLUMNS),
                ("jcvi_hy4p", JCVI_COLUMNS),
            ):
                rows = read_tsv(paths[role])
                for row in rows:
                    if row["canonical_chromosome"] == "Chr01":
                        row["canonical_chromosome"] = "Chr02"
                    elif row["canonical_chromosome"] == "Chr02":
                        row["canonical_chromosome"] = "Chr01"
                write_tsv(paths[role], columns, rows)
            refresh_provenance(root, paths)
            with self.assertRaisesRegex(
                ChromosomeAssignmentError, "canonical labels conflict with the frozen"
            ):
                run_assignment(root, paths)

    def test_reference_registry_checksum_substitution_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            altered_registry = root / "altered_reference_assets.tsv"
            rows = read_tsv(REFERENCE_ASSET_REGISTRY)
            rows[0]["sha256"] = hashlib.sha256(b"wrong-reference-genome").hexdigest()
            write_tsv(
                altered_registry,
                tuple(rows[0]),
                rows,
            )
            with self.assertRaisesRegex(
                ChromosomeAssignmentError, "checksum conflicts with the frozen policy"
            ):
                run_assignment(
                    root, paths, reference_asset_registry=altered_registry
                )

    def test_generation_parameter_tamper_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            sidecar = paths["nucleotide_hy4a_provenance"]
            rewrite_json(
                sidecar,
                lambda document: document["generation_parameters"].__setitem__(
                    "minimum_mapq", 19
                ),
            )
            with self.assertRaisesRegex(
                ChromosomeProvenanceError, "frozen generation parameter"
            ):
                run_assignment(root, paths)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            sidecar = paths["nucleotide_hy4a_provenance"]
            text = sidecar.read_text(encoding="utf-8").replace(
                '"maximum_de": 0.15',
                '"maximum_de": 0.1500000000000000000000000000000000001',
                1,
            )
            sidecar.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(
                ChromosomeProvenanceError, "frozen generation parameter"
            ):
                run_assignment(root, paths)

    def test_minimap2_version_and_exact_argv_are_bound_in_both_reports(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            sidecar = paths["nucleotide_hy4a_provenance"]
            rewrite_json(
                sidecar,
                lambda document: document["generation_parameters"].__setitem__(
                    "minimap2_version", "2.28-r1208"
                ),
            )
            with self.assertRaisesRegex(
                ChromosomeProvenanceError, "frozen generation parameter"
            ):
                run_assignment(root, paths)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            sidecar = paths["nucleotide_hy4a_provenance"]
            rewrite_json(
                sidecar,
                lambda document: document["generation_parameters"][
                    "required_paf_tags"
                ].__setitem__(0, "tp:A:S"),
            )
            with self.assertRaisesRegex(
                ChromosomeProvenanceError, "frozen generation parameter"
            ):
                run_assignment(root, paths)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            report = root / "nucleotide_hy4a.upstream_validation.json"
            rewrite_json(
                report,
                lambda document: document["generation_parameters"][
                    "minimap2_command_template"
                ].insert(1, "-t"),
            )
            rebind_report(paths["nucleotide_hy4a_provenance"], report)
            with self.assertRaisesRegex(
                ChromosomeProvenanceError, "frozen generation parameter"
            ):
                run_assignment(root, paths)

    def test_upstream_report_requires_exact_all_true_check_set(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            report = root / "nucleotide_hy4a.upstream_validation.json"
            rewrite_json(
                report,
                lambda document: document["checks"].__setitem__(
                    "divergence_filter_reconciled", False
                ),
            )
            rebind_report(paths["nucleotide_hy4a_provenance"], report)
            with self.assertRaisesRegex(
                ChromosomeProvenanceError, "every upstream validation check must be true"
            ):
                run_assignment(root, paths)

    def test_upstream_report_workflow_version_is_pinned(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            report = root / "nucleotide_hy4a.upstream_validation.json"
            rewrite_json(
                report,
                lambda document: document.__setitem__(
                    "workflow_version", "attacker_version"
                ),
            )
            rebind_report(paths["nucleotide_hy4a_provenance"], report)
            with self.assertRaisesRegex(ChromosomeProvenanceError, "workflow_version"):
                run_assignment(root, paths)

    def test_duplicate_json_keys_and_nonfinite_json_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            sidecar = paths["nucleotide_hy4a_provenance"]
            text = sidecar.read_text(encoding="utf-8")
            sidecar.write_text(
                text.replace(
                    '"status": "PASS",',
                    '"status": "PASS",\n  "status": "PASS",',
                    1,
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ChromosomeProvenanceError, "repeats key"):
                run_assignment(root, paths)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            sidecar = paths["nucleotide_hy4a_provenance"]
            text = sidecar.read_text(encoding="utf-8")
            sidecar.write_text(
                text.replace('"minimum_mapq": 20', '"minimum_mapq": NaN', 1),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ChromosomeProvenanceError, "non-finite"):
                run_assignment(root, paths)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            sidecar = paths["nucleotide_hy4a_provenance"]
            rewrite_json(
                sidecar,
                lambda document: document["upstream_validation_report"].__setitem__(
                    "basename", "bad\x00report.json"
                ),
            )
            with self.assertRaisesRegex(ChromosomeProvenanceError, "control character"):
                run_assignment(root, paths)

    def test_weighted_divergence_boundary_is_exact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["nucleotide_hy4a"])
            rows[0]["weighted_divergence"] = "0.15"
            coverage = float(rows[0]["query_coverage"])
            rows[0]["score"] = f"{coverage * 0.85:.12g}"
            write_tsv(paths["nucleotide_hy4a"], NUCLEOTIDE_COLUMNS, rows)
            refresh_provenance(root, paths)
            self.assertEqual(run_assignment(root, paths).publication_gate, "PASS")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["nucleotide_hy4a"])
            rows[0]["weighted_divergence"] = "0.1500000001"
            coverage = float(rows[0]["query_coverage"])
            rows[0]["score"] = f"{coverage * (1.0 - 0.1500000001):.12g}"
            write_tsv(paths["nucleotide_hy4a"], NUCLEOTIDE_COLUMNS, rows)
            refresh_provenance(root, paths)
            with self.assertRaisesRegex(
                ChromosomeAssignmentError, "exceeds frozen maximum_de"
            ):
                run_assignment(root, paths)

    def test_low_jcvi_absolute_support_fails_unit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["jcvi_hy4a"])
            for row in rows:
                if (
                    row["query_chromosome"] == "PubChr01"
                    and row["reference_chromosome"] == "Chr08A"
                ):
                    row["query_anchored_genes"] = "29"
                    row["query_gene_coverage"] = "0.29"
                    row["reference_anchored_genes"] = "29"
                    row["reference_gene_coverage"] = "0.29"
                    row["unique_anchor_pairs"] = "29"
                    row["score"] = "0.29"
            write_tsv(paths["jcvi_hy4a"], JCVI_COLUMNS, rows)
            refresh_provenance(root, paths)
            result = run_assignment(root, paths)
            self.assertEqual(result.status, "LOW_JCVI_ABSOLUTE_SUPPORT")
            evidence = read_tsv(
                root / "assignment" / "jcvi_hy4a.assignment_evidence.tsv"
            )
            row = next(item for item in evidence if item["query_chromosome"] == "PubChr01")
            self.assertEqual(row["assigned_unique_anchor_pairs"], "29")
            self.assertEqual(row["absolute_support_gate"], "false")

    def test_jcvi_pair_count_cannot_exceed_cartesian_bound(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["jcvi_hy4a"])
            rows[0]["unique_anchor_pairs"] = "5"
            write_tsv(paths["jcvi_hy4a"], JCVI_COLUMNS, rows)
            refresh_provenance(root, paths)
            with self.assertRaisesRegex(
                ChromosomeAssignmentError, "Cartesian gene-pair bound"
            ):
                run_assignment(root, paths)

    def test_exact_integer_parser_has_no_float_round_trip(self) -> None:
        self.assertEqual(
            parse_exact_int("9007199254740993", label="large"),
            9007199254740993,
        )
        for invalid in ("1e3", "1.0", "+1", "01", "-1", "9" * 5000):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ChromosomeProvenanceError):
                    parse_exact_int(invalid, label="invalid")

    def test_decimal_bounds_use_original_lexemes_without_float_rounding(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["nucleotide_hy4a"])
            selected = next(
                row
                for row in rows
                if row["query_chromosome"] == "PubChr01"
                and row["reference_chromosome"] == "Chr08A"
            )
            selected["weighted_divergence"] = (
                "0.1500000000000000000000000000000000001"
            )
            write_tsv(paths["nucleotide_hy4a"], NUCLEOTIDE_COLUMNS, rows)
            refresh_provenance(root, paths)
            with self.assertRaisesRegex(
                ChromosomeAssignmentError, "exceeds frozen maximum_de"
            ):
                run_assignment(root, paths)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["nucleotide_hy4a"])
            selected = next(
                row
                for row in rows
                if row["query_chromosome"] == "PubChr01"
                and row["reference_chromosome"] == "Chr08A"
            )
            selected["weighted_divergence"] = "-1e-9999"
            selected["score"] = selected["query_coverage"]
            write_tsv(paths["nucleotide_hy4a"], NUCLEOTIDE_COLUMNS, rows)
            refresh_provenance(root, paths)
            with self.assertRaisesRegex(
                ChromosomeAssignmentError, "nonnegative decimal spelling"
            ):
                run_assignment(root, paths)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            rows = read_tsv(paths["jcvi_hy4a"])
            selected = next(
                row
                for row in rows
                if row["query_chromosome"] == "PubChr01"
                and row["reference_chromosome"] == "Chr08A"
            )
            selected.update(
                {
                    "query_anchored_genes": "100",
                    "query_gene_coverage": "1.0000000000000000000000000000000000001",
                    "reference_anchored_genes": "100",
                    "reference_gene_coverage": "1",
                    "unique_anchor_pairs": "100",
                    "score": "1",
                }
            )
            write_tsv(paths["jcvi_hy4a"], JCVI_COLUMNS, rows)
            refresh_provenance(root, paths)
            with self.assertRaisesRegex(
                ChromosomeAssignmentError, r"must lie in \[0,1\]"
            ):
                run_assignment(root, paths)

    def test_ratio_gate_is_inclusive_and_tolerance_bounded(self) -> None:
        policy = load_homology_policy(PARAMETERS)

        def matrix(second: float) -> ScoreMatrix:
            scores = {
                ("q1", "r1"): 0.9,
                ("q1", "r2"): second,
                ("q2", "r1"): second,
                ("q2", "r2"): 0.9,
            }
            cells = {
                key: MatrixCell(
                    query=key[0],
                    reference=key[1],
                    canonical="Chr01" if key[1] == "r1" else "Chr02",
                    score=value,
                    raw={"reciprocal_coverage": "0.9", "orientation": "+"},
                )
                for key, value in scores.items()
            }
            return ScoreMatrix(
                role="synthetic",
                kind="nucleotide",
                source=Path("synthetic.tsv"),
                columns=(),
                queries=("q1", "q2"),
                references=("r1", "r2"),
                canonical_by_reference={"r1": "Chr01", "r2": "Chr02"},
                cells=cells,
                normalized_rows=(),
            )

        self.assertTrue(solve_score_matrix(matrix(0.6), policy)["q1"].separation_gate)
        self.assertTrue(
            solve_score_matrix(matrix(0.6000000001), policy)["q1"].separation_gate
        )
        self.assertFalse(
            solve_score_matrix(matrix(0.60000001), policy)["q1"].separation_gate
        )

    def test_total_order_is_hash_seed_independent_under_natural_key_collision(self) -> None:
        code = (
            "from geneloss_repro.chromosome_assignment import _total_key;"
            "print('|'.join(sorted({'PubChr1','PubChr01','PubChr001'},key=_total_key)))"
        )
        observed = set()
        for seed in ("1", "2", "99", "random"):
            environment = dict(os.environ)
            environment["PYTHONHASHSEED"] = seed
            environment["PYTHONPATH"] = str(ROOT / "src")
            completed = subprocess.run(
                [sys.executable, "-c", code],
                cwd=ROOT,
                env=environment,
                check=True,
                capture_output=True,
                text=True,
            )
            observed.add(completed.stdout.strip())
        self.assertEqual(len(observed), 1)

    def test_policy_rejects_weakened_filter_and_failure_contracts(self) -> None:
        source = PARAMETERS.read_text(encoding="utf-8")
        replacements = (
            ("minimap2_primary_only = true", "minimap2_primary_only = false"),
            ("minimum_mapq = 20", "minimum_mapq = 19"),
            ("minimum_alignment_block_bp = 10000", "minimum_alignment_block_bp = 9999"),
            ("maximum_de = 0.15", "maximum_de = 0.16"),
            (
                "maximum_de = 0.15",
                "maximum_de = 0.149999999999999999999999999999999999",
            ),
            ("maximum_de = 0.15", "maximum_de = " + "9" * 500),
            (
                'minimap2_version = "2.28-r1209"',
                'minimap2_version = "2.28-r1208"',
            ),
            (
                '"--secondary=no", "-c"',
                '"--secondary=yes", "-c"',
            ),
            (
                "minimum_assigned_reciprocal_nucleotide_coverage = 0.05",
                "minimum_assigned_reciprocal_nucleotide_coverage = 0.0",
            ),
            (
                "minimum_unique_anchor_pairs = 30",
                "minimum_unique_anchor_pairs = " + "9" * 5000,
            ),
            ("arithmetic_tolerance = 1e-9", "arithmetic_tolerance = 1e-6"),
            ('aligner = "LAST"', 'aligner = "UNVETTED"'),
            ('database_type = "protein"', 'database_type = "nucleotide"'),
            ("cscore = 0.7", "cscore = 0.01"),
            (
                "cscore = 0.7",
                "cscore = 0.7000000000000000000000000000000000001",
            ),
            ("tandem_nmax = 10", "tandem_nmax = 999"),
            ("maximum_gene_distance = 20", "maximum_gene_distance = 999"),
            ("minimum_anchor_block_size = 4", "minimum_anchor_block_size = 1"),
            (
                'coverage_anchor_source = "raw JCVI anchors"',
                'coverage_anchor_source = "screened anchors"',
            ),
            (
                'coverage_arithmetic = "interval union in both query and reference directions; never raw PAF-row sums"',
                'coverage_arithmetic = "raw row sums"',
            ),
            (
                'failure_policy = "retain provisional PubChr labels and stop; never force a final Chr label"',
                'failure_policy = "force a label"',
            ),
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for index, (old, new) in enumerate(replacements):
                weakened = root / f"weakened_{index}.toml"
                weakened.write_text(source.replace(old, new, 1), encoding="utf-8")
                with self.subTest(replacement=new):
                    with self.assertRaises(ChromosomeAssignmentError):
                        load_homology_policy(weakened)
            nonfinite = root / "nonfinite.toml"
            nonfinite.write_text(
                source.replace("maximum_de = 0.15", "maximum_de = nan", 1),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ChromosomeAssignmentError, "must be finite"):
                load_homology_policy(nonfinite)

    def test_trusted_repository_commit_is_full_and_recorded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            with self.assertRaisesRegex(
                ChromosomeAssignmentError, "trusted_repository_commit"
            ):
                run_assignment(root, paths, trusted_repository_commit="deadbeef")
            self.assertFalse((root / "assignment").exists())

    def test_verify_snapshot_detects_path_change_when_explicitly_requested(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input.tsv"
            path.write_text("first\n", encoding="utf-8")
            snapshot = capture_snapshot(path)
            path.write_text("second\n", encoding="utf-8")
            with self.assertRaisesRegex(ChromosomeProvenanceError, "changed after snapshot"):
                verify_snapshot(snapshot)

    def test_snapshot_capture_rejects_fifo_without_blocking(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fifo = Path(temporary) / "input.fifo"
            os.mkfifo(fifo)
            code = (
                "from geneloss_repro.chromosome_provenance import "
                "ChromosomeProvenanceError,capture_snapshot;"
                f"p={str(fifo)!r};"
                "\ntry: capture_snapshot(p)"
                "\nexcept ChromosomeProvenanceError as e: print(e)"
                "\nelse: raise SystemExit(3)"
            )
            environment = dict(os.environ)
            environment["PYTHONPATH"] = str(ROOT / "src")
            completed = subprocess.run(
                [sys.executable, "-c", code],
                cwd=ROOT,
                env=environment,
                check=False,
                capture_output=True,
                text=True,
                timeout=2,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("regular file", completed.stdout)

    def test_later_path_mutation_cannot_change_captured_hash_bound_result(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            matrix = paths["nucleotide_hy4a"]
            captured_sha256 = hashlib.sha256(matrix.read_bytes()).hexdigest()
            original = chromosome_assignment_module._rename_no_replace

            def mutate_after_capture(source: Path, destination: Path) -> None:
                matrix.write_text("changed after capture\n", encoding="utf-8")
                original(source, destination)

            with mock.patch.object(
                chromosome_assignment_module,
                "_rename_no_replace",
                side_effect=mutate_after_capture,
            ):
                result = run_assignment(root, paths)
            self.assertEqual(result.publication_gate, "PASS")
            manifest = json.loads(
                (root / "assignment" / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                manifest["inputs"]["matrix_roles"]["nucleotide_hy4a"]["matrix"][
                    "sha256"
                ],
                captured_sha256,
            )
            self.assertNotEqual(hashlib.sha256(matrix.read_bytes()).hexdigest(), captured_sha256)

    def test_existing_lock_and_dangling_output_symlink_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            lock = root / ".assignment.chromosome_assignment.lock"
            lock.write_text("owned\n", encoding="utf-8")
            with self.assertRaisesRegex(ChromosomeAssignmentError, "owns the output lock"):
                run_assignment(root, paths)
            self.assertEqual(lock.read_text(encoding="utf-8"), "owned\n")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            output = root / "assignment"
            output.symlink_to(root / "missing-directory", target_is_directory=True)
            with self.assertRaisesRegex(ChromosomeAssignmentError, "already exists"):
                run_assignment(root, paths)
            self.assertTrue(output.is_symlink())

    def test_atomic_no_replace_preserves_racing_destination(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = make_inputs(root)
            original = chromosome_assignment_module._rename_no_replace

            def inject_destination(source: Path, destination: Path) -> None:
                destination.mkdir()
                (destination / "sentinel.txt").write_text("keep\n", encoding="utf-8")
                original(source, destination)

            with mock.patch.object(
                chromosome_assignment_module,
                "_rename_no_replace",
                side_effect=inject_destination,
            ):
                with self.assertRaisesRegex(ChromosomeAssignmentError, "refusing overwrite"):
                    run_assignment(root, paths)
            self.assertEqual(
                (root / "assignment" / "sentinel.txt").read_text(encoding="utf-8"),
                "keep\n",
            )
            self.assertFalse(list(root.glob(".assignment.staging.*")))


if __name__ == "__main__":
    unittest.main()
