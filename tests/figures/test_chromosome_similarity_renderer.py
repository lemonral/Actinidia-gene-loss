"""Tests for the author-approved chromosome similarity naming figure."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

try:
    import matplotlib  # noqa: F401
except ImportError:
    matplotlib = None


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "figures" / "render_chromosome_similarity_qc.py"
SPEC = importlib.util.spec_from_file_location("render_chromosome_similarity_qc", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fixture(root: Path) -> tuple[Path, Path]:
    metadata = root / "metadata.tsv"
    write_tsv(metadata, MODULE.METADATA_REQUIRED, [{
        "assembly_unit_id":"act_eriantha_hap1_2026", "biological_species":"Actinidia eriantha",
        "haplotype_or_subgenome":"HAP1", "accession":"GCA_TEST", "decision_status":"current",
    }])
    assignment_root = root / "assignments"
    unit = assignment_root / "act_eriantha_hap1_2026"
    unit.mkdir(parents=True)
    label_rows = []
    diagnostic_rows = []
    for index in range(1, 4):
        query = f"PubChr{index:02d}"
        final = f"Chr{index:02d}"
        strict = "true" if index != 1 else "false"
        label_rows.append({
            "query_chromosome":query, "final_chromosome":final,
            "coordinate_reference":"act_chinensis_hongyang_v4_hy4a",
            "assignment_method":"global_one_to_one_maximum_nucleotide_similarity",
            "assigned_score":"0.1", "reciprocal_coverage":"0.08", "orientation_to_hy4a":"+",
            "hy4p_and_jcvi_agree":"true", "strict_homology_gates_pass":strict,
            "confidence_flag":"HIGH" if strict == "true" else "SUPPORTED",
        })
        diagnostic_rows.append({
            "query_chromosome":query, "diagnostic_candidate":final, "nucleotide_hy4a":final,
            "jcvi_hy4a":final, "nucleotide_hy4p":final, "jcvi_hy4p":final,
            "nucleotide_hy4a_reciprocal_coverage":"0.08", "jcvi_hy4a_score":"0.65",
            "orientation_hy4a":"+", "all_four_matrix_gates":strict,
            "diagnostic_status":"PASS_DIAGNOSTIC" if strict == "true" else "REVIEW_REQUIRED",
            "failure_reasons":"" if strict == "true" else "LOW_RECIPROCAL_COVERAGE",
        })
    write_tsv(unit / "similarity_label_map.tsv", MODULE.LABEL_COLUMNS, label_rows)
    write_tsv(unit / "prefinal_assignment_diagnostic.tsv", MODULE.DIAGNOSTIC_COLUMNS, diagnostic_rows)
    (unit / "diagnostic.json").write_text(json.dumps({
        "assembly_unit_id":"act_eriantha_hap1_2026", "chromosome_naming_status":"PASS_LABELS",
        "chromosome_naming_policy":"HY4A global one-to-one maximum nucleotide similarity; absolute support is QC only",
    }))
    return metadata, assignment_root


class ChromosomeSimilarityRendererTest(unittest.TestCase):
    def test_support_does_not_block_complete_labels(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            metadata, assignment_root = fixture(Path(temporary))
            rows, validation = MODULE.prepare(
                metadata_path=metadata, assignment_roots=[assignment_root],
                expected_unit_count=1, expected_chromosome_count=3,
            )
        self.assertEqual(rows[0]["named_chromosome_count"], 3)
        self.assertEqual(rows[0]["strict_all_four_gate_count"], 2)
        self.assertEqual(rows[0]["qc_support_only_count"], 1)
        self.assertEqual(rows[0]["four_matrix_agreement_count"], 3)
        self.assertEqual(validation["status"], "PASS_CHROMOSOME_SIMILARITY_NAMING_PUBLICATION")
        self.assertIn(r"$\mathit{A.\ eriantha}$", rows[0]["display_label"])
        self.assertIn(r"$\mathrm{HAP1}$", rows[0]["display_label"])

    def test_nonbijective_map_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            metadata, assignment_root = fixture(Path(temporary))
            path = assignment_root / "act_eriantha_hap1_2026" / "similarity_label_map.tsv"
            rows = MODULE.read_tsv(path, MODULE.LABEL_COLUMNS, "labels")
            rows[2]["final_chromosome"] = "Chr02"
            write_tsv(path, MODULE.LABEL_COLUMNS, rows)
            with self.assertRaisesRegex(MODULE.ChromosomeSimilarityPlotError, "bijection"):
                MODULE.prepare(
                    metadata_path=metadata, assignment_roots=[assignment_root],
                    expected_unit_count=1, expected_chromosome_count=3,
                )

    @unittest.skipUnless(matplotlib is not None, "optional matplotlib is not installed")
    def test_render_complete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata, assignment_root = fixture(root)
            bundle = MODULE.render_bundle(
                metadata_path=metadata, assignment_roots=[assignment_root], output_dir=root / "out",
                basename="chromosome_naming", expected_unit_count=1, expected_chromosome_count=3, dpi=90,
            )
            self.assertTrue(bundle.png.read_bytes().startswith(b"\x89PNG"))
            self.assertIn("QC only", bundle.caption.read_text())


if __name__ == "__main__":
    unittest.main()
