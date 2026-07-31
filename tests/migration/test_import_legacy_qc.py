"""Tests for the fail-closed legacy assembly-QC importer."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).parents[2]
SCRIPT = PROJECT_ROOT / "scripts" / "migration" / "import_legacy_qc.py"
PUBLIC_MAPPING = PROJECT_ROOT / "config" / "legacy_qc_sample_map.tsv"

SPEC = importlib.util.spec_from_file_location("import_legacy_qc", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
IMPORTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IMPORTER)


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def basic_row(sample: str) -> dict[str, str]:
    row = {column: "0" for column in IMPORTER.BASIC_STATS_COLUMNS}
    row.update(
        {
            "sample": sample,
            "current_or_alternative": "legacy_test",
            "accession": f"ACC_{sample}",
            "source_url": "https://example.org/public-record",
            "genome_path": f"/private/server/{sample}.fna",
            "gff_path": f"/private/server/{sample}.gff3",
            "protein_path": f"/private/server/{sample}.faa",
        }
    )
    return row


def busco_row(sample: str, mode: str, version: str = "5.8.2") -> dict[str, str]:
    row = {column: "0" for column in IMPORTER.BUSCO_COLUMNS}
    row.update(
        {
            "sample": sample,
            "busco_version": version,
            "dataset": "embryophyta_odb10",
            "dataset_creation_date": "2024-01-08",
            "mode": mode,
            "input_path": f"/private/server/{sample}.input",
            "C_percent": "98.0",
            "S_percent": "70.0",
            "D_percent": "28.0",
            "F_percent": "1.0",
            "M_percent": "1.0",
            "n": "1614",
            "short_summary_path": f"/private/server/{sample}.short_summary.txt",
        }
    )
    return row


def write_mapping(path: Path, samples: list[tuple[str, str]]) -> None:
    write_tsv(
        path,
        IMPORTER.MAPPING_COLUMNS,
        [
            {"legacy_sample": legacy_sample, "assembly_unit_id": assembly_unit_id}
            for legacy_sample, assembly_unit_id in samples
        ],
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class LegacyQCImportTest(unittest.TestCase):
    def run_import(
        self,
        root: Path,
        *,
        mapping: Path,
        basic_paths: list[Path],
        genome_paths: list[Path],
        protein_paths: list[Path],
    ) -> subprocess.CompletedProcess[str]:
        command = [sys.executable, str(SCRIPT), "--mapping", str(mapping)]
        for path in basic_paths:
            command.extend(("--basic-stats", str(path)))
        for path in genome_paths:
            command.extend(("--genome-busco", str(path)))
        for path in protein_paths:
            command.extend(("--protein-busco", str(path)))
        command.extend(("--output-dir", str(root / "public")))
        return subprocess.run(command, text=True, capture_output=True, check=False)

    def make_single_sample_inputs(
        self,
        root: Path,
        sample: str = "Old_A",
    ) -> tuple[Path, Path, Path]:
        basic = root / "basic.tsv"
        genome = root / "genome.tsv"
        protein = root / "protein.tsv"
        write_tsv(basic, IMPORTER.BASIC_STATS_COLUMNS, [basic_row(sample)])
        write_tsv(genome, IMPORTER.BUSCO_COLUMNS, [busco_row(sample, "euk_genome_min")])
        write_tsv(protein, IMPORTER.BUSCO_COLUMNS, [busco_row(sample, "proteins")])
        return basic, genome, protein

    def test_multiple_inputs_emit_path_free_traceable_tables(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            mapping = root / "mapping.tsv"
            # Reverse mapping order and IDs to confirm deterministic output sorting.
            write_mapping(mapping, [("Old_A", "unit_z"), ("Old_B", "unit_a")])

            basic_a = root / "current_basic.tsv"
            basic_b = root / "alternative_basic.tsv"
            genome_a = root / "current_genome.tsv"
            genome_b = root / "alternative_genome.tsv"
            protein_a = root / "current_protein.tsv"
            protein_b = root / "alternative_protein.tsv"
            write_tsv(basic_a, IMPORTER.BASIC_STATS_COLUMNS, [basic_row("Old_A")])
            write_tsv(basic_b, IMPORTER.BASIC_STATS_COLUMNS, [basic_row("Old_B")])
            write_tsv(genome_a, IMPORTER.BUSCO_COLUMNS, [busco_row("Old_A", "euk_genome_min")])
            write_tsv(genome_b, IMPORTER.BUSCO_COLUMNS, [busco_row("Old_B", "euk_genome_min")])
            write_tsv(protein_a, IMPORTER.BUSCO_COLUMNS, [busco_row("Old_A", "proteins")])
            write_tsv(protein_b, IMPORTER.BUSCO_COLUMNS, [busco_row("Old_B", "proteins")])

            completed = self.run_import(
                root,
                mapping=mapping,
                basic_paths=[basic_a, basic_b],
                genome_paths=[genome_a, genome_b],
                protein_paths=[protein_a, protein_b],
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            summary = json.loads(completed.stdout)
            self.assertEqual(summary["sample_count"], 2)

            output_dir = root / "public"
            basic_output = output_dir / IMPORTER.OUTPUT_FILENAMES["basic_stats"]
            genome_output = output_dir / IMPORTER.OUTPUT_FILENAMES["genome_busco"]
            protein_output = output_dir / IMPORTER.OUTPUT_FILENAMES["protein_busco"]
            for path in (basic_output, genome_output, protein_output):
                self.assertTrue(path.is_file())
                self.assertNotIn("/private/server/", path.read_text(encoding="utf-8"))

            with basic_output.open(encoding="utf-8", newline="") as handle:
                basic_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["assembly_unit_id"] for row in basic_rows], ["unit_a", "unit_z"])
            self.assertNotIn("genome_path", basic_rows[0])
            self.assertNotIn("gff_path", basic_rows[0])
            self.assertNotIn("protein_path", basic_rows[0])
            provenance = {row["legacy_sample"]: row for row in basic_rows}
            self.assertEqual(provenance["Old_A"]["source_basename"], basic_a.name)
            self.assertEqual(provenance["Old_A"]["source_sha256"], sha256(basic_a))
            self.assertEqual(provenance["Old_B"]["source_basename"], basic_b.name)
            self.assertEqual(provenance["Old_B"]["source_sha256"], sha256(basic_b))

            with genome_output.open(encoding="utf-8", newline="") as handle:
                genome_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertNotIn("input_path", genome_rows[0])
            self.assertNotIn("short_summary_path", genome_rows[0])

    def test_unmapped_sample_fails_before_writing_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            mapping = root / "mapping.tsv"
            write_mapping(mapping, [("Old_A", "unit_a")])
            basic, genome, protein = self.make_single_sample_inputs(root, "Unknown")
            completed = self.run_import(
                root,
                mapping=mapping,
                basic_paths=[basic],
                genome_paths=[genome],
                protein_paths=[protein],
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("unmapped legacy sample", completed.stderr)
            self.assertFalse((root / "public").exists())

    def test_duplicate_sample_across_split_files_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            mapping = root / "mapping.tsv"
            write_mapping(mapping, [("Old_A", "unit_a")])
            basic, genome, protein = self.make_single_sample_inputs(root)
            duplicate_basic = root / "basic_duplicate.tsv"
            write_tsv(duplicate_basic, IMPORTER.BASIC_STATS_COLUMNS, [basic_row("Old_A")])
            completed = self.run_import(
                root,
                mapping=mapping,
                basic_paths=[basic, duplicate_basic],
                genome_paths=[genome],
                protein_paths=[protein],
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("duplicate sample", completed.stderr)
            self.assertFalse((root / "public").exists())

    def test_schema_change_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            mapping = root / "mapping.tsv"
            write_mapping(mapping, [("Old_A", "unit_a")])
            basic, genome, protein = self.make_single_sample_inputs(root)
            reduced_columns = tuple(
                column for column in IMPORTER.BASIC_STATS_COLUMNS if column != "protein_path"
            )
            complete_row = basic_row("Old_A")
            write_tsv(
                basic,
                reduced_columns,
                [{column: complete_row[column] for column in reduced_columns}],
            )
            completed = self.run_import(
                root,
                mapping=mapping,
                basic_paths=[basic],
                genome_paths=[genome],
                protein_paths=[protein],
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("schema mismatch", completed.stderr)
            self.assertFalse((root / "public").exists())

    def test_busco_version_difference_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            mapping = root / "mapping.tsv"
            write_mapping(mapping, [("Old_A", "unit_a")])
            basic, genome, protein = self.make_single_sample_inputs(root)
            write_tsv(
                protein,
                IMPORTER.BUSCO_COLUMNS,
                [busco_row("Old_A", "proteins", version="5.7.1")],
            )
            completed = self.run_import(
                root,
                mapping=mapping,
                basic_paths=[basic],
                genome_paths=[genome],
                protein_paths=[protein],
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("incompatible BUSCO version/dataset signatures", completed.stderr)
            self.assertFalse((root / "public").exists())

    def test_public_mapping_has_exactly_the_expected_26_labels(self) -> None:
        with PUBLIC_MAPPING.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 26)
        labels = {row["legacy_sample"] for row in rows}
        units = {row["assembly_unit_id"] for row in rows}
        self.assertEqual(len(labels), 26)
        self.assertEqual(len(units), 26)
        self.assertTrue(
            {
                "Actinidia_eriantha_MDHAPA_matched",
                "Actinidia_eriantha_MHT_RefSeq",
                "Actinidia_rufa_ARU_r1.0",
                "Actinidia_rufa_Fuchu_full",
            }.issubset(labels)
        )


if __name__ == "__main__":
    unittest.main()
