"""Focused tests for the generic bidirectional JCVI renderer."""

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


PROJECT_ROOT = Path(__file__).parents[2]
SCRIPT = PROJECT_ROOT / "scripts" / "figures" / "make_jcvi_figure.py"
SPEC = importlib.util.spec_from_file_location("make_jcvi_figure", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
JCVI = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(JCVI)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def metadata_rows() -> list[dict[str, str]]:
    return [
        {
            "assembly_unit_id": "clem_scandens_reference",
            "biological_species": "Clematoclethra scandens",
            "haplotype_or_subgenome": "reference",
            "accession": "REF_1",
            "assembly_scope": "chromosome_anchored",
            "decision_status": "current",
            "comparison_role": "reference",
        },
        {
            "assembly_unit_id": "act_deliciosa_hap_a",
            "biological_species": "Actinidia deliciosa",
            "haplotype_or_subgenome": "A",
            "accession": "GWH_TEST_A",
            "assembly_scope": "chromosome_anchored",
            "decision_status": "current",
            "comparison_role": "target",
        },
        {
            "assembly_unit_id": "act_eriantha_hap1",
            "biological_species": "Actinidia eriantha",
            "haplotype_or_subgenome": "HAP1",
            "accession": "GCA_000001.1",
            "assembly_scope": "anchored_pseudochromosomes",
            "decision_status": "candidate",
            "comparison_role": "target",
        },
        {
            "assembly_unit_id": "act_rufa_subset",
            "biological_species": "Actinidia rufa",
            "haplotype_or_subgenome": "unphased",
            "accession": "GCA_000002.1",
            "assembly_scope": "chromosome_only_subset",
            "decision_status": "excluded",
            "comparison_role": "target",
        },
    ]


def coverage_row(
    assembly_unit_id: str,
    decision_status: str,
    *,
    reference_gene_covered: int,
    target_gene_covered: int,
    reference_sequence_covered: int,
    target_sequence_covered: int,
) -> dict[str, str]:
    reference_gene_total = 100
    target_gene_total = 200
    reference_sequence_total = 1000
    target_sequence_total = 2000
    return {
        "assembly_unit_id": assembly_unit_id,
        "reference_assembly_unit_id": "clem_scandens_reference",
        "decision_status": decision_status,
        "reference_gene_covered": str(reference_gene_covered),
        "reference_gene_total": str(reference_gene_total),
        "reference_gene_coverage_percent": f"{reference_gene_covered * 100 / reference_gene_total:.6f}",
        "target_gene_covered": str(target_gene_covered),
        "target_gene_total": str(target_gene_total),
        "target_gene_coverage_percent": f"{target_gene_covered * 100 / target_gene_total:.6f}",
        "reference_sequence_covered_bp": str(reference_sequence_covered),
        "reference_sequence_total_bp": str(reference_sequence_total),
        "reference_sequence_coverage_percent": f"{reference_sequence_covered * 100 / reference_sequence_total:.6f}",
        "target_sequence_covered_bp": str(target_sequence_covered),
        "target_sequence_total_bp": str(target_sequence_total),
        "target_sequence_coverage_percent": f"{target_sequence_covered * 100 / target_sequence_total:.6f}",
    }


def coverage_rows() -> list[dict[str, str]]:
    return [
        coverage_row(
            "act_deliciosa_hap_a",
            "current",
            reference_gene_covered=80,
            target_gene_covered=150,
            reference_sequence_covered=700,
            target_sequence_covered=1300,
        ),
        coverage_row(
            "act_eriantha_hap1",
            "candidate",
            reference_gene_covered=60,
            target_gene_covered=100,
            reference_sequence_covered=500,
            target_sequence_covered=900,
        ),
        coverage_row(
            "act_rufa_subset",
            "excluded",
            reference_gene_covered=40,
            target_gene_covered=60,
            reference_sequence_covered=300,
            target_sequence_covered=500,
        ),
    ]


def make_inputs(root: Path) -> tuple[Path, Path]:
    metadata = root / "comparison_metadata.tsv"
    coverage = root / "jcvi_coverage.tsv"
    write_tsv(metadata, metadata_rows())
    write_tsv(coverage, coverage_rows())
    return metadata, coverage


class JCVIPublicationRendererTest(unittest.TestCase):
    def test_prepare_requires_and_retains_all_four_directions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = make_inputs(Path(temporary_directory))
            rows, validation = JCVI.prepare_plot_rows(*paths)

        self.assertEqual(
            [row["assembly_unit_id"] for row in rows],
            ["act_deliciosa_hap_a", "act_eriantha_hap1", "act_rufa_subset"],
        )
        self.assertEqual(rows[0]["reference_gene_coverage_percent"], 80.0)
        self.assertEqual(rows[0]["target_gene_coverage_percent"], 75.0)
        self.assertEqual(rows[0]["reference_sequence_coverage_percent"], 70.0)
        self.assertEqual(rows[0]["target_sequence_coverage_percent"], 65.0)
        self.assertEqual(rows[2]["decision_status"], "excluded")
        self.assertEqual(rows[2]["assembly_scope"], "chromosome_only_subset")
        self.assertIn(r"$\mathit{A.\ rufa}$", rows[2]["display_label"])
        self.assertIn(
            r"$\mathit{C.\ scandens}$", rows[0]["reference_display_label"]
        )
        self.assertEqual(rows[0]["reference_decision_status"], "current")
        self.assertEqual(validation["excluded_target_count"], 1)
        self.assertEqual(
            validation["checks"][
                "all_four_bidirectional_gene_and_sequence_metrics_present"
            ],
            "pass",
        )

    def test_missing_direction_and_bad_denominator_reconciliation_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            metadata, coverage = make_inputs(root)
            incomplete_rows = coverage_rows()
            for row in incomplete_rows:
                row.pop("target_sequence_coverage_percent")
            write_tsv(coverage, incomplete_rows)
            with self.assertRaisesRegex(JCVI.JCVIPlotError, "missing JCVI coverage"):
                JCVI.prepare_plot_rows(metadata, coverage)

            inconsistent = coverage_rows()
            inconsistent[0]["reference_gene_coverage_percent"] = "79.0"
            write_tsv(coverage, inconsistent)
            with self.assertRaisesRegex(JCVI.JCVIPlotError, "does not reconcile"):
                JCVI.prepare_plot_rows(metadata, coverage)

    def test_target_id_and_decision_status_reconciliation_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            metadata, coverage = make_inputs(root)
            reduced = coverage_rows()[:-1]
            write_tsv(coverage, reduced)
            with self.assertRaisesRegex(JCVI.JCVIPlotError, "differs from target metadata"):
                JCVI.prepare_plot_rows(metadata, coverage)

            mismatched = coverage_rows()
            mismatched[1]["decision_status"] = "excluded"
            write_tsv(coverage, mismatched)
            with self.assertRaisesRegex(JCVI.JCVIPlotError, "differs from metadata"):
                JCVI.prepare_plot_rows(metadata, coverage)

    @unittest.skipUnless(matplotlib is not None, "optional matplotlib is not installed")
    def test_render_publishes_complete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            metadata, coverage = make_inputs(root)
            bundle = JCVI.render_bundle(
                metadata_path=metadata,
                coverage_path=coverage,
                output_dir=root / "publication_jcvi",
                basename="jcvi_bidirectional_coverage",
                dpi=100,
            )
            expected = {
                "jcvi_bidirectional_coverage.png",
                "jcvi_bidirectional_coverage.pdf",
                "jcvi_bidirectional_coverage.plot_data.tsv",
                "jcvi_bidirectional_coverage.caption.txt",
                "jcvi_bidirectional_coverage.validation.json",
                "jcvi_bidirectional_coverage.manifest.json",
            }
            self.assertEqual({path.name for path in bundle.directory.iterdir()}, expected)
            self.assertTrue(bundle.png.read_bytes().startswith(b"\x89PNG"))
            self.assertTrue(bundle.pdf.read_bytes().startswith(b"%PDF"))
            caption = bundle.caption.read_text()
            self.assertIn("No one-directional percentage", caption)
            validation = json.loads(bundle.validation.read_text())
            self.assertEqual(validation["target_assembly_unit_count"], 3)
            plot_header = bundle.plot_data.read_text().splitlines()[0].split("\t")
            for field in (
                "reference_gene_coverage_percent",
                "target_gene_coverage_percent",
                "reference_sequence_coverage_percent",
                "target_sequence_coverage_percent",
            ):
                self.assertIn(field, plot_header)


if __name__ == "__main__":
    unittest.main()
