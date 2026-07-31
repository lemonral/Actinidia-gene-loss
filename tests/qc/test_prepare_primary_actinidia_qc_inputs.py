"""Tests for primary Actinidia QC input preparation."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "qc" / "prepare_primary_actinidia_qc_inputs.py"


def write(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class PrimaryQCInputTests(unittest.TestCase):
    def test_rekeys_sources_and_binds_primary_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sys.path.insert(0, str(SCRIPT.parent))
            import prepare_primary_actinidia_qc_inputs as module

            basic = {column: "0" for column in module.BASIC_COLUMNS}
            basic.update({
                "sample": "source", "current_or_alternative": "current", "accession": "ACC",
                "source_url": "https://example.org", "genome_path": "/private/g.fa",
                "gff_path": "/private/g.gff3", "protein_path": "/private/p.faa",
                "genome_sequence_count": "29", "genome_total_bp": "1000",
                "genome_ungapped_bp": "1000", "genome_longest_bp": "100",
                "genome_n50_bp": "100", "genome_l50": "5", "gff_gene_count": "10",
                "gff_mrna_or_transcript_count": "12", "protein_sequence_count": "12",
            })
            busco = {column: "0" for column in module.BUSCO_COLUMNS}
            busco.update({
                "sample": "source", "busco_version": "5.8.2", "dataset": "embryophyta_odb10",
                "dataset_creation_date": "2024-01-08", "mode": "euk_genome_min",
                "input_path": "/private/g.fa", "C_percent": "90.0", "S_percent": "70.0",
                "D_percent": "20.0", "F_percent": "5.0", "M_percent": "5.0", "n": "100",
                "S_count": "70", "D_count": "20", "F_count": "5", "M_count": "5",
                "short_summary_path": "/private/summary.txt",
            })
            protein = dict(busco)
            protein["mode"] = "proteins"
            write(root / "basic.tsv", module.BASIC_COLUMNS, [basic])
            write(root / "genome.tsv", module.BUSCO_COLUMNS, [busco])
            write(root / "protein.tsv", module.BUSCO_COLUMNS, [protein])
            primary = {
                "status": "PASS", "publication_gate": "PASS",
                "counts": {"chromosome_sequences": 29, "selected_genes": 9,
                           "selected_transcripts": 9, "invalid_coding_genes": 1},
                "policy": {"fallback_selection": "longest_valid_spliced_CDS"},
            }
            (root / "primary.json").write_text(json.dumps(primary), encoding="utf-8")
            selection = {column: "value" for column in module.SELECTION_COLUMNS}
            selection.update({
                "assembly_unit_id": "unit", "source_sample": "source",
                "source_accession": "ACC", "accession": "CORRECTED_ACC",
                "basic_stats_table": "basic.tsv", "genome_busco_table": "genome.tsv",
                "protein_busco_table": "protein.tsv",
                "primary_annotation_manifest": "primary.json", "biological_species": "Actinidia test",
                "individual_id": "i", "haplotype_or_subgenome": "A",
                "publisher_assembly_scope": "full", "local_qc_scope": "29 chromosomes",
                "publisher_assembly_provenance": "https://example.org/g", "publisher_annotation_provenance": "https://example.org/a",
                "publisher_protein_provenance": "https://example.org/p", "decision_reason": "selected",
            })
            write(root / "selection.tsv", module.SELECTION_COLUMNS, [selection])
            output = root / "output"
            run = subprocess.run([
                sys.executable, str(SCRIPT), "--selection", str(root / "selection.tsv"),
                "--data-root", str(root), "--expected-units", "1", "--output-dir", str(output),
            ], text=True, capture_output=True, check=False)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            with (output / "basic_stats.tsv").open(newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["sample"], "unit")
            self.assertEqual(row["accession"], "CORRECTED_ACC")
            with (output / "analysis_scope_primary_annotation.tsv").open(newline="") as handle:
                scope = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(scope["analysis_primary_gene_count"], "9")
            self.assertEqual(scope["invalid_coding_gene_count"], "1")
            self.assertEqual(json.loads((output / "run_manifest.json").read_text())["status"], "PASS")

    def test_legacy_public_row_is_normalized_without_private_paths(self) -> None:
        sys.path.insert(0, str(SCRIPT.parent))
        import prepare_primary_actinidia_qc_inputs as module
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            row = {column: "0" for column in module.LEGACY_BASIC_COLUMNS}
            row.update({
                "assembly_unit_id": "legacy", "legacy_sample": "Old_name",
                "current_or_alternative": "current", "accession": "ACC",
                "source_url": "https://example.org", "source_basename": "old.tsv",
                "source_sha256": "a" * 64,
            })
            path = root / "legacy.tsv"
            write(path, module.LEGACY_BASIC_COLUMNS, [row])
            converted = module.normalized_source_row(
                path, module.BASIC_COLUMNS, "legacy", table_kind="basic"
            )
            self.assertEqual(converted["sample"], "legacy")
            self.assertEqual(converted["genome_path"], "legacy_exact_bound/legacy.genome")

    def test_missing_or_duplicate_source_fails_without_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selection = root / "selection.tsv"
            selection.write_text("assembly_unit_id\nunit\n", encoding="utf-8")
            output = root / "output"
            run = subprocess.run([
                sys.executable, str(SCRIPT), "--selection", str(selection),
                "--data-root", str(root), "--expected-units", "1", "--output-dir", str(output),
            ], text=True, capture_output=True, check=False)
            self.assertNotEqual(run.returncode, 0)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
