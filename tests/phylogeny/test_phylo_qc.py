"""Focused tests for generic taxon and terminal-manifest phylogeny QC."""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
PHYLOGENY_SCRIPTS = PROJECT / "scripts" / "phylogeny"
sys.path.insert(0, str(PHYLOGENY_SCRIPTS))

from phylo_io import DataError, load_sequence_id_map, load_terminal_manifest  # noqa: E402


TERMINAL_HEADER = (
    "terminal_id\tbiological_species\tgrouping\tsource_fasta_stem\t"
    "canonical_tree_label\tis_root_outgroup\tinclude_species_tree\tidentity_status\n"
)


class TaxonRegistryTest(unittest.TestCase):
    def test_current_registry_contains_all_actinidia_species_as_candidates(self) -> None:
        manifest = load_terminal_manifest(PROJECT / "config" / "phylogeny" / "taxa.tsv")
        self.assertEqual(manifest.schema, "taxon_registry")
        actinidia = {
            row["biological_species"]
            for row in manifest.all_rows
            if row["biological_species"].startswith("Actinidia ")
        }
        self.assertEqual(
            actinidia,
            {
                "Actinidia arguta",
                "Actinidia chinensis",
                "Actinidia deliciosa",
                "Actinidia eriantha",
                "Actinidia hemsleyana",
                "Actinidia latifolia",
                "Actinidia longicarpa",
                "Actinidia macrosperma",
                "Actinidia polygama",
                "Actinidia reticulata",
                "Actinidia rufa",
                "Actinidia x zhejiangensis",
            },
        )
        candidate_ids = {row["sample_id"] for row in manifest.candidates}
        self.assertTrue({"act_deliciosa", "act_eriantha", "act_rufa"} <= candidate_ids)
        self.assertFalse(candidate_ids & {row["sample_id"] for row in manifest.selected})
        self.assertTrue(
            {
                "rubia_or_cinchona",
                "non_rhododendron_ericaceae",
                "actinidiaceae_sister_outgroup",
                "ampelocissus_or_nothocissus",
            }
            <= {row["sample_id"] for row in manifest.all_rows}
        )


class GenericTerminalManifestTest(unittest.TestCase):
    @staticmethod
    def write_manifest(path: Path) -> None:
        path.write_text(
            TERMINAL_HEADER
            + "alpha_hap1\tActinidia alpha\tact_alpha\talpha_hap1\tActinidia_alpha_HAP1\tfalse\ttrue\tconfirmed\n"
            + "alpha_hap2\tActinidia alpha\tact_alpha\talpha_hap2\tActinidia_alpha_HAP2\tfalse\ttrue\tconfirmed\n"
            + "outgroup\tOutgroup beta\toutgroup_beta\toutgroup\tOutgroup_beta\ttrue\ttrue\tconfirmed\n"
            + "candidate\tActinidia gamma\tact_gamma\tcandidate\tActinidia_gamma\tfalse\tcandidate\tpending_qc\n",
            encoding="utf-8",
        )

    def run_qc(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(PHYLOGENY_SCRIPTS / "phylo_qc.py"), *arguments],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_counts_and_multi_terminal_group_are_derived_from_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            manifest_path = root / "terminals.tsv"
            report = root / "metadata.tsv"
            self.write_manifest(manifest_path)

            completed = self.run_qc(
                "metadata",
                "--manifest",
                str(manifest_path),
                "--require-root",
                "--report",
                str(report),
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            text = report.read_text(encoding="utf-8")
            self.assertIn("found 3 selected terminals; count is derived from the manifest", text)
            self.assertIn("found 2 selected biological-species groups", text)
            self.assertIn("act_alpha (2)", text)
            self.assertIn("configured root outgroup: outgroup", text)
            self.assertIn("candidate rows remain outside the selected tree: candidate", text)

    def test_tree_closure_uses_selected_rows_not_candidate_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            manifest_path = root / "terminals.tsv"
            tree = root / "tree.nwk"
            report = root / "tree_qc.tsv"
            self.write_manifest(manifest_path)
            tree.write_text(
                "((Actinidia_alpha_HAP1,Actinidia_alpha_HAP2),Outgroup_beta);\n",
                encoding="utf-8",
            )

            completed = self.run_qc(
                "tree",
                "--tree",
                str(tree),
                "--manifest",
                str(manifest_path),
                "--report",
                str(report),
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            text = report.read_text(encoding="utf-8")
            self.assertIn("parsed 3 terminal labels; expected 3 from the manifest", text)
            self.assertIn("tree resolves one terminal for every configured sample", text)
            self.assertNotIn("candidate", text)

    def test_terminal_manifest_requires_explicit_species_grouping(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "bad.tsv"
            path.write_text(
                "terminal_id\tsource_fasta_stem\tcanonical_tree_label\n"
                "a\ta\tA\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(DataError, "biological-species grouping"):
                load_terminal_manifest(path)


class OrthoFinderThreeIdentifierTest(unittest.TestCase):
    def test_numeric_sequence_prefixes_resolve_through_species_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            sequence_ids = root / "SequenceIDs.txt"
            species_ids = root / "SpeciesIDs.txt"
            sequence_ids.write_text("0_0: gene_a\n1_0: gene_b\n", encoding="utf-8")
            species_ids.write_text("0: alpha.faa\n1: beta.faa\n", encoding="utf-8")
            self.assertEqual(
                load_sequence_id_map(sequence_ids, {"gene_a", "gene_b"}),
                {"gene_a": "alpha.faa", "gene_b": "beta.faa"},
            )

    def test_unknown_numeric_species_prefix_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            sequence_ids = root / "SequenceIDs.txt"
            species_ids = root / "SpeciesIDs.txt"
            sequence_ids.write_text("2_0: gene_a\n", encoding="utf-8")
            species_ids.write_text("0: alpha.faa\n", encoding="utf-8")
            with self.assertRaisesRegex(DataError, "species prefix 2"):
                load_sequence_id_map(sequence_ids, {"gene_a"})


if __name__ == "__main__":
    unittest.main()
