"""Focused tests for the generic assembly/annotation QC renderer."""

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
SCRIPT = PROJECT_ROOT / "scripts" / "figures" / "make_qc_figure.py"
SPEC = importlib.util.spec_from_file_location("make_qc_figure", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
QC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(QC)


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
            "assembly_unit_id": "act_deliciosa_hap_a",
            "biological_species": "Actinidia deliciosa",
            "haplotype_or_subgenome": "A",
            "accession": "GWH_TEST_A",
            "assembly_scope": "chromosome_anchored",
            "decision_status": "current",
        },
        {
            "assembly_unit_id": "act_eriantha_hap1",
            "biological_species": "Actinidia eriantha",
            "haplotype_or_subgenome": "HAP1",
            "accession": "GCA_000001.1",
            "assembly_scope": "anchored_pseudochromosomes",
            "decision_status": "candidate",
        },
        {
            "assembly_unit_id": "act_rufa_subset",
            "biological_species": "Actinidia rufa",
            "haplotype_or_subgenome": "unphased",
            "accession": "GCA_000002.1",
            "assembly_scope": "chromosome_only_subset",
            "decision_status": "excluded",
        },
    ]


def basic_rows() -> list[dict[str, str]]:
    return [
        {
            "assembly_unit_id": "act_deliciosa_hap_a",
            "genome_total_bp": "650000000",
            "gff_gene_count": "42000",
        },
        {
            "assembly_unit_id": "act_eriantha_hap1",
            "genome_total_bp": "620000000",
            "gff_gene_count": "40500",
        },
        {
            "assembly_unit_id": "act_rufa_subset",
            "genome_total_bp": "500000000",
            "gff_gene_count": "30000",
        },
    ]


def busco_rows(mode: str, percentages: tuple[str, str, str]) -> list[dict[str, str]]:
    identifiers = (
        "act_deliciosa_hap_a",
        "act_eriantha_hap1",
        "act_rufa_subset",
    )
    return [
        {
            "assembly_unit_id": assembly_unit_id,
            "busco_version": "5.8.2",
            "dataset": "embryophyta_odb10",
            "dataset_creation_date": "2024-01-08",
            "mode": mode,
            "C_percent": percentage,
        }
        for assembly_unit_id, percentage in zip(identifiers, percentages)
    ]


def make_inputs(root: Path) -> tuple[Path, Path, Path, Path]:
    metadata = root / "assembly_metadata.tsv"
    basic = root / "basic_statistics.tsv"
    genome = root / "genome_busco.tsv"
    protein = root / "protein_busco.tsv"
    write_tsv(metadata, metadata_rows())
    write_tsv(basic, basic_rows())
    write_tsv(genome, busco_rows("euk_genome_min", ("98.0", "97.5", "86.9")))
    write_tsv(protein, busco_rows("proteins", ("96.0", "95.0", "58.0")))
    return metadata, basic, genome, protein


class QCPublicationRendererTest(unittest.TestCase):
    def test_prepare_retains_all_statuses_and_metadata_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            paths = make_inputs(Path(temporary_directory))
            rows, validation = QC.prepare_plot_rows(*paths)

        self.assertEqual(
            [row["assembly_unit_id"] for row in rows],
            ["act_deliciosa_hap_a", "act_eriantha_hap1", "act_rufa_subset"],
        )
        self.assertEqual(
            [row["decision_status"] for row in rows],
            ["current", "candidate", "excluded"],
        )
        self.assertEqual(rows[2]["assembly_scope"], "chromosome_only_subset")
        self.assertIn(r"$\mathit{A.\ deliciosa}$", rows[0]["display_label"])
        self.assertIn(r"$\mathrm{current}$", rows[0]["display_label"])
        self.assertEqual(validation["excluded_assembly_unit_count"], 1)
        self.assertEqual(
            validation["checks"]["all_metadata_rows_plotted_including_excluded"],
            "pass",
        )

    def test_identifier_and_busco_signature_mismatches_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            metadata, basic, genome, protein = make_inputs(root)
            incomplete = basic_rows()[:-1]
            write_tsv(basic, incomplete)
            with self.assertRaisesRegex(QC.QCPlotError, "differs from metadata"):
                QC.prepare_plot_rows(metadata, basic, genome, protein)

            write_tsv(basic, basic_rows())
            inconsistent = busco_rows("proteins", ("96.0", "95.0", "58.0"))
            inconsistent[2]["busco_version"] = "5.7.1"
            write_tsv(protein, inconsistent)
            with self.assertRaisesRegex(QC.QCPlotError, "one nonempty"):
                QC.prepare_plot_rows(metadata, basic, genome, protein)

    @unittest.skipUnless(matplotlib is not None, "optional matplotlib is not installed")
    def test_render_publishes_complete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            metadata, basic, genome, protein = make_inputs(root)
            bundle = QC.render_bundle(
                metadata_path=metadata,
                basic_stats_path=basic,
                genome_busco_path=genome,
                protein_busco_path=protein,
                output_dir=root / "publication_qc",
                basename="assembly_annotation_qc",
                dpi=100,
            )
            expected = {
                "assembly_annotation_qc.png",
                "assembly_annotation_qc.pdf",
                "assembly_annotation_qc.plot_data.tsv",
                "assembly_annotation_qc.caption.txt",
                "assembly_annotation_qc.validation.json",
                "assembly_annotation_qc.manifest.json",
            }
            self.assertEqual({path.name for path in bundle.directory.iterdir()}, expected)
            self.assertTrue(bundle.png.read_bytes().startswith(b"\x89PNG"))
            self.assertTrue(bundle.pdf.read_bytes().startswith(b"%PDF"))
            self.assertIn("including units later excluded", bundle.caption.read_text())
            validation = json.loads(bundle.validation.read_text())
            self.assertEqual(validation["assembly_unit_count"], 3)
            manifest_text = bundle.manifest.read_text()
            self.assertNotIn(str(root), manifest_text)


if __name__ == "__main__":
    unittest.main()
