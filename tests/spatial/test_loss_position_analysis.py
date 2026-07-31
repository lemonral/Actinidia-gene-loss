"""Tests for the manifest-driven multi-assembly position analysis."""

from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path

from geneloss_repro.io_utils import SchemaError, read_tsv
from geneloss_repro.spatial import analyze_loss_positions


UNITS = [
    ("act_deliciosa_a_legacy", "Actinidia deliciosa", "A"),
    ("act_deliciosa_b_legacy", "Actinidia deliciosa", "B"),
    ("act_deliciosa_c_legacy", "Actinidia deliciosa", "C"),
    ("act_deliciosa_d_legacy", "Actinidia deliciosa", "D"),
    ("act_deliciosa_e_legacy", "Actinidia deliciosa", "E"),
    ("act_deliciosa_f_legacy", "Actinidia deliciosa", "F"),
    ("act_zhejiangensis_a_legacy", "Actinidia zhejiangensis", "A"),
    ("act_zhejiangensis_b_legacy", "Actinidia zhejiangensis", "B"),
]

GFF = """##gff-version 3
Chr1\ttest\tgene\t1\t10\t.\t+\t.\tID=G1
Chr1\ttest\tgene\t21\t30\t.\t+\t.\tID=G2
Chr1\ttest\tgene\t51\t60\t.\t+\t.\tID=G3
Chr1\ttest\tgene\t91\t100\t.\t+\t.\tID=G4
"""


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(text)


class LossPositionAnalysisTest(unittest.TestCase):
    def _fixture(self, root: Path) -> dict[str, Path]:
        manifest = root / "assemblies.tsv"
        calls = root / "calls.tsv.gz"
        coordinates = root / "coordinates.tsv"
        centromeres = root / "centromeres.tsv"

        manifest_lines = [
            "assembly_unit_id\tbiological_species\thaplotype_or_subgenome\tassembly_scope\tgenome\tgff\tinclude_spatial"
        ]
        call_lines = ["target_haplotype\treference_gene_id\tclassification"]
        coordinate_lines = [
            "target_haplotype\treference_gene_id\tclassification\ttarget_chromosome\ttarget_start\ttarget_end"
        ]
        centromere_lines = [
            "assembly_unit_id\tchromosome\tcentromere_start\tcentromere_end\tevidence_source"
        ]
        positions = [1, 20, 21, 40, 41, 60, 100]
        for index, (unit, species, suffix) in enumerate(UNITS):
            genome = root / f"{unit}.fa.gz"
            gff = root / f"{unit}.gff3.gz"
            _write_gzip(genome, ">Chr1\n" + "A" * 100 + "\n")
            _write_gzip(gff, GFF)
            manifest_lines.append(
                f"{unit}\t{species}\t{suffix}\tchromosome_partition\t{genome.name}\t{gff.name}\ttrue"
            )
            gene = f"R{index + 1}"
            # Keep the final unit in scope even though this toy unit has zero
            # positive calls.  This exercises fail-closed manifest scope without
            # silently dropping an assembly unit from denominator outputs.
            classification = "not_called_loss" if index == len(UNITS) - 1 else "pseudogenized"
            call_lines.append(f"{unit}\t{gene}\t{classification}")
            if classification == "pseudogenized":
                position = positions[index]
                coordinate_lines.append(
                    f"{unit}\t{gene}\tpseudogenized\tChr1\t{position}\t{position}"
                )
            centromere_lines.append(f"{unit}\tChr1\t45\t55\tindependent_test_interval")

        manifest.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
        _write_gzip(calls, "\n".join(call_lines) + "\n")
        coordinates.write_text("\n".join(coordinate_lines) + "\n", encoding="utf-8")
        centromeres.write_text("\n".join(centromere_lines) + "\n", encoding="utf-8")
        return {
            "manifest": manifest,
            "calls": calls,
            "coordinates": coordinates,
            "centromeres": centromeres,
        }

    def test_all_named_subgenomes_primary_distances_and_legacy_mode(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fixture = self._fixture(root)
            outputs = analyze_loss_positions(
                fixture["calls"],
                fixture["coordinates"],
                fixture["manifest"],
                root / "result",
                analysis_label="primary_nonshared_pseudogenized",
                number_of_bins=5,
                centromeres_path=fixture["centromeres"],
                require_complete_centromeres=True,
                legacy_reproduction=True,
                manifest_include_column="include_spatial",
            )

            positions = read_tsv(outputs["positions"])
            self.assertEqual(len(positions), 7)
            self.assertEqual(
                {row["assembly_unit_id"] for row in positions},
                {unit for unit, _, _ in UNITS[:-1]},
            )
            self.assertEqual(
                {row["haplotype_or_subgenome"] for row in positions if "deliciosa" in row["assembly_unit_id"]},
                {"A", "B", "C", "D", "E", "F"},
            )
            self.assertTrue(
                all(row["centromere_status"] == "independently_supplied_interval" for row in positions)
            )
            normalized = [float(row["normalized_end_distance_0_end_1_center"]) for row in positions]
            self.assertTrue(all(0 <= value <= 1 for value in normalized))
            self.assertEqual(normalized[0], 0.0)
            self.assertEqual(normalized[-1], 0.0)
            bins_by_gene = {row["reference_gene_id"]: row["equal_width_bin"] for row in positions}
            self.assertEqual(
                [bins_by_gene[f"R{index}"] for index in range(1, 8)],
                ["1", "1", "2", "2", "3", "3", "5"],
            )

            equal_bins = read_tsv(outputs["equal_width_bins"])
            self.assertEqual(len(equal_bins), len(UNITS) * 5)
            self.assertEqual(
                [(row["bin_start_1based_inclusive"], row["bin_end_1based_inclusive"])
                 for row in equal_bins[:5]],
                [("1", "20"), ("21", "40"), ("41", "60"), ("61", "80"), ("81", "100")],
            )
            self.assertEqual(sum(int(row["positive_loss_fragments"]) for row in equal_bins), 7)
            self.assertEqual(sum(int(row["gff_gene_opportunities"]) for row in equal_bins), 32)
            self.assertEqual(
                {row["analysis_mode"] for row in equal_bins},
                {"primary_mutually_exclusive_equal_width"},
            )

            end_bins = read_tsv(outputs["end_distance_bins"])
            self.assertEqual(len(end_bins), len(UNITS) * 5)
            self.assertEqual(sum(int(row["positive_loss_fragments"]) for row in end_bins), 7)

            units = read_tsv(outputs["assembly_units"])
            self.assertEqual(len(units), len(UNITS))
            zero = next(row for row in units if row["assembly_unit_id"] == UNITS[-1][0])
            self.assertEqual(zero["positive_loss_fragments"], "0")
            self.assertEqual(zero["haplotype_or_subgenome"], "B")

            legacy = read_tsv(outputs["legacy_reproduction"])
            self.assertEqual(len(legacy), len(UNITS) * 5)
            self.assertEqual(
                {row["analysis_mode"] for row in legacy},
                {"manuscript_era_nested_midpoint_reproduction_only"},
            )
            self.assertEqual({row["intervals_are_mutually_exclusive"] for row in legacy}, {"false"})
            self.assertEqual({row["inferential_test_permitted"] for row in legacy}, {"false"})
            self.assertEqual(
                [(row["interval_start_legacy_coordinate"], row["interval_end_legacy_coordinate"])
                 for row in legacy[:5]],
                [("40", "60"), ("30", "70"), ("20", "80"), ("10", "90"), ("0", "100")],
            )

            metadata = json.loads(outputs["metadata"].read_text(encoding="utf-8"))
            self.assertEqual(metadata["reconciliation"]["assembly_unit_count"], 8)
            self.assertEqual(metadata["reconciliation"]["positive_call_count"], 7)
            self.assertTrue(metadata["legacy_reproduction"]["enabled"])
            self.assertIn("0 is a chromosome end", metadata["normalized_end_distance_definition"])

    def test_missing_coordinate_fails_before_writing_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fixture = self._fixture(root)
            lines = fixture["coordinates"].read_text(encoding="utf-8").splitlines()
            fixture["coordinates"].write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
            output = root / "must_not_exist"
            with self.assertRaisesRegex(SchemaError, "positive calls lack coordinates"):
                analyze_loss_positions(
                    fixture["calls"],
                    fixture["coordinates"],
                    fixture["manifest"],
                    output,
                    analysis_label="mismatch_test",
                    centromeres_path=fixture["centromeres"],
                    manifest_include_column="include_spatial",
                )
            self.assertFalse(output.exists())

    def test_centromere_scope_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fixture = self._fixture(root)
            with fixture["centromeres"].open("a", encoding="utf-8") as handle:
                handle.write("act_unknown\tChr1\t45\t55\tinvalid_scope\n")
            with self.assertRaisesRegex(SchemaError, "outside analysis scope"):
                analyze_loss_positions(
                    fixture["calls"],
                    fixture["coordinates"],
                    fixture["manifest"],
                    root / "result",
                    analysis_label="centromere_scope_test",
                    centromeres_path=fixture["centromeres"],
                    manifest_include_column="include_spatial",
                )


if __name__ == "__main__":
    unittest.main()
