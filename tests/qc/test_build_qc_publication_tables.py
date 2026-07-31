"""Focused tests for the fail-closed QC publication-table builder."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).parents[2]
SCRIPT = PROJECT_ROOT / "scripts" / "qc" / "build_qc_publication_tables.py"
SPEC = importlib.util.spec_from_file_location("build_qc_publication_tables", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
QC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(QC)

FIGURE_SCRIPT = PROJECT_ROOT / "scripts" / "figures" / "make_qc_figure.py"
FIGURE_SPEC = importlib.util.spec_from_file_location("make_qc_figure_for_builder", FIGURE_SCRIPT)
assert FIGURE_SPEC is not None and FIGURE_SPEC.loader is not None
QC_FIGURE = importlib.util.module_from_spec(FIGURE_SPEC)
FIGURE_SPEC.loader.exec_module(QC_FIGURE)


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metadata_rows() -> list[dict[str, str]]:
    specs = (
        ("act_deliciosa_a", "qc_a", "Actinidia deliciosa", "ADM", "A", "ACC_A", "current", ""),
        ("act_eriantha_hap1", "qc_e1", "Actinidia eriantha", "wild_1", "HAP1", "ACC_E1", "candidate", "paired-haplotype assessment pending"),
        ("act_rufa_old", "qc_r", "Actinidia rufa", "Fuchu", "unphased", "ACC_R", "excluded", "excluded from the primary cohort after scope review"),
    )
    return [
        {
            "assembly_unit_id": unit,
            "qc_sample": sample,
            "biological_species": species,
            "individual_id": individual,
            "haplotype_or_subgenome": haplotype,
            "accession": accession,
            "decision_status": status,
            "publisher_assembly_scope": "publisher_full_primary_assembly",
            "local_qc_scope": "matched_chromosome_and_annotation_scope",
            "publisher_assembly_provenance": f"https://example.org/assembly/{accession}",
            "publisher_annotation_provenance": f"https://example.org/gff/{accession}",
            "publisher_protein_provenance": f"https://example.org/protein/{accession}",
            "decision_reason": reason,
        }
        for unit, sample, species, individual, haplotype, accession, status, reason in specs
    ]


def basic_row(sample: str, accession: str, offset: int = 0) -> dict[str, str]:
    row = {column: "0" for column in QC.BASIC_INPUT_COLUMNS}
    row.update(
        {
            "sample": sample,
            "current_or_alternative": "current",
            "accession": accession,
            "source_url": f"https://example.org/source/{accession}",
            "genome_path": f"/private/genomes/{sample}.fa",
            "gff_path": f"/private/annotations/{sample}.gff3",
            "protein_path": f"/private/proteins/{sample}.faa",
            "genome_sequence_count": "2",
            "genome_total_bp": str(1000 + offset),
            "genome_ungapped_bp": "1000",
            "genome_n_bp": "10",
            "genome_n_percent": "1.000000",
            "genome_gc_bp": "400",
            "genome_gc_percent": "40.000000",
            "genome_longest_bp": str(600 + offset),
            "genome_n50_bp": str(600 + offset),
            "genome_l50": "1",
            "gff_feature_rows": "10",
            "gff_invalid_rows": "0",
            "gff_gene_count": "2",
            "gff_mrna_count": "2",
            "gff_transcript_count": "0",
            "gff_mrna_or_transcript_count": "2",
            "gff_cds_count": "3",
            "gff_exon_count": "3",
            "protein_sequence_count": "2",
            "protein_empty_sequence_count": "0",
            "protein_total_aa": "200",
            "protein_longest_aa": "120",
            "protein_n50_aa": "120",
            "protein_l50": "1",
            "protein_internal_stop_record_count": "0",
            "protein_terminal_stop_record_count": "2",
            "protein_internal_stop_character_count": "0",
            "protein_nonstandard_character_record_count": "0",
            "protein_nonstandard_character_count": "0",
        }
    )
    return row


def busco_row(sample: str, mode: str) -> dict[str, str]:
    return {
        "sample": sample,
        "busco_version": "5.8.2",
        "dataset": "embryophyta_odb10",
        "dataset_creation_date": "2024-01-08",
        "mode": mode,
        "input_path": f"/private/busco_inputs/{sample}.fa",
        "C_percent": "90.0",
        "S_percent": "70.0",
        "D_percent": "20.0",
        "F_percent": "5.0",
        "M_percent": "5.0",
        "n": "100",
        "C_count": "",
        "S_count": "70",
        "D_count": "20",
        "F_count": "5",
        "M_count": "5",
        "short_summary_path": f"/private/busco/{sample}/short_summary.txt",
    }


def make_inputs(root: Path) -> tuple[Path, Path, Path, Path]:
    metadata = root / "assembly_decisions.tsv"
    basic = root / "basic_stats.tsv"
    genome = root / "genome_busco.tsv"
    protein = root / "protein_busco.tsv"
    metadata_data = metadata_rows()
    write_tsv(metadata, QC.METADATA_INPUT_COLUMNS, metadata_data)
    write_tsv(
        basic,
        QC.BASIC_INPUT_COLUMNS,
        [
            basic_row(row["qc_sample"], row["accession"], index)
            for index, row in enumerate(metadata_data)
        ],
    )
    write_tsv(
        genome,
        QC.BUSCO_INPUT_COLUMNS,
        [busco_row(row["qc_sample"], "euk_genome_min") for row in metadata_data],
    )
    write_tsv(
        protein,
        QC.BUSCO_INPUT_COLUMNS,
        [busco_row(row["qc_sample"], "proteins") for row in metadata_data],
    )
    return metadata, basic, genome, protein


def namespace(paths: tuple[Path, Path, Path, Path], output: Path) -> argparse.Namespace:
    metadata, basic, genome, protein = paths
    return argparse.Namespace(
        metadata=metadata,
        basic_stats=basic,
        genome_busco=genome,
        protein_busco=protein,
        output_dir=output,
    )


class QCPublicationTableBuilderTest(unittest.TestCase):
    def test_complete_bundle_is_path_free_retains_statuses_and_feeds_renderer(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            paths = make_inputs(root)
            output = root / "publication_qc"
            summary = QC.run(namespace(paths, output))

            self.assertEqual(summary["assembly_unit_count"], 3)
            self.assertEqual(
                {path.name for path in output.iterdir()},
                set(QC.OUTPUT_FILENAMES.values()),
            )
            metadata = read_tsv(output / QC.OUTPUT_FILENAMES["metadata"])
            self.assertEqual(
                [row["decision_status"] for row in metadata],
                ["current", "candidate", "excluded"],
            )
            self.assertNotIn("qc_sample", metadata[0])
            genome = read_tsv(output / QC.OUTPUT_FILENAMES["genome_busco"])
            self.assertEqual(genome[0]["C_count"], "90")
            combined = read_tsv(output / QC.OUTPUT_FILENAMES["combined"])
            self.assertEqual(combined[2]["decision_status"], "excluded")
            self.assertEqual(combined[0]["genome_busco_C_percent"], "90.0")
            for path in output.iterdir():
                text = path.read_text(encoding="utf-8")
                self.assertNotIn("/private/", text)
                self.assertNotIn(str(root), text)

            plot_rows, validation = QC_FIGURE.prepare_plot_rows(
                output / QC.OUTPUT_FILENAMES["metadata"],
                output / QC.OUTPUT_FILENAMES["basic"],
                output / QC.OUTPUT_FILENAMES["genome_busco"],
                output / QC.OUTPUT_FILENAMES["protein_busco"],
            )
            self.assertEqual(len(plot_rows), 3)
            self.assertEqual(validation["excluded_assembly_unit_count"], 1)
            report = json.loads(
                (output / QC.OUTPUT_FILENAMES["validation"]).read_text()
            )
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["decision_status_counts"]["candidate"], 1)

    def test_identifier_mismatch_fails_before_output_publication(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            paths = make_inputs(root)
            protein_rows = read_tsv(paths[3])
            protein_rows[-1]["sample"] = "undeclared_sample"
            write_tsv(paths[3], QC.BUSCO_INPUT_COLUMNS, protein_rows)
            output = root / "must_not_exist"
            with self.assertRaisesRegex(QC.QCPublicationError, "differs from metadata"):
                QC.run(namespace(paths, output))
            self.assertFalse(output.exists())

    def test_busco_signature_and_arithmetic_mismatches_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            paths = make_inputs(root)
            protein_rows = read_tsv(paths[3])
            protein_rows[-1]["dataset_creation_date"] = "2023-01-01"
            write_tsv(paths[3], QC.BUSCO_INPUT_COLUMNS, protein_rows)
            with self.assertRaisesRegex(QC.QCPublicationError, "not uniform"):
                QC.run(namespace(paths, root / "bad_signature"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            paths = make_inputs(root)
            genome_rows = read_tsv(paths[2])
            genome_rows[0]["C_percent"] = "89.0"
            write_tsv(paths[2], QC.BUSCO_INPUT_COLUMNS, genome_rows)
            with self.assertRaisesRegex(QC.QCPublicationError, "inconsistent"):
                QC.run(namespace(paths, root / "bad_arithmetic"))

    def test_basic_stats_arithmetic_and_nonempty_destination_are_refused(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            paths = make_inputs(root)
            basic_rows = read_tsv(paths[1])
            basic_rows[0]["gff_mrna_or_transcript_count"] = "1"
            write_tsv(paths[1], QC.BASIC_INPUT_COLUMNS, basic_rows)
            with self.assertRaisesRegex(QC.QCPublicationError, "not mRNA plus transcript"):
                QC.run(namespace(paths, root / "bad_basic"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            paths = make_inputs(root)
            output = root / "occupied"
            output.mkdir()
            (output / "accepted.txt").write_text("keep\n", encoding="utf-8")
            with self.assertRaisesRegex(QC.QCPublicationError, "nonempty"):
                QC.run(namespace(paths, output))
            self.assertEqual((output / "accepted.txt").read_text(), "keep\n")


if __name__ == "__main__":
    unittest.main()
