"""Tests for the dated-tree and CAFE5 Base publication renderer."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

try:
    import matplotlib  # noqa: F401
except ImportError:
    matplotlib = None


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "figures" / "render_phylogeny_cafe.py"
SPEC = importlib.util.spec_from_file_location("render_phylogeny_cafe", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def inputs(root: Path) -> dict[str, Path]:
    paths = {name: root / name for name in (
        "terminals.tsv", "tree.tre", "node_ages.tsv", "calibrations.tsv",
        "model.tsv", "clades.tsv", "node_tree.tre", "cafe_validation.json"
    )}
    write_tsv(paths["terminals.tsv"], MODULE.TERMINAL_COLUMNS, [
        {"terminal_id":"act_a_A","biological_species":"Actinidia arguta","grouping":"a","source_fasta_stem":"a","canonical_tree_label":"Actinidia_arguta_A","is_root_outgroup":"false","include_species_tree":"true","identity_status":"confirmed"},
        {"terminal_id":"act_b","biological_species":"Actinidia chinensis","grouping":"b","source_fasta_stem":"b","canonical_tree_label":"Actinidia_chinensis","is_root_outgroup":"false","include_species_tree":"true","identity_status":"confirmed"},
        {"terminal_id":"vitis","biological_species":"Vitis vinifera","grouping":"v","source_fasta_stem":"v","canonical_tree_label":"Vitis_vinifera","is_root_outgroup":"true","include_species_tree":"true","identity_status":"confirmed"},
    ])
    paths["tree.tre"].write_text("((Actinidia_arguta_A:2,Actinidia_chinensis:2):3,Vitis_vinifera:5);\n")
    paths["node_tree.tre"].write_text(
        "((Actinidia_arguta_A<1>,Actinidia_chinensis<2>)<4>,"
        "Vitis_vinifera<3>)<5>;\n"
    )
    write_tsv(paths["node_ages.tsv"], MODULE.NODE_AGE_COLUMNS, [
        {"node_id":"n4","descendant_tip_count":"3","descendant_tips":"Actinidia_arguta_A,Actinidia_chinensis,Vitis_vinifera","mean_ma":"5","q025_ma":"4","q975_ma":"6","chain1_mean_ma":"5","chain2_mean_ma":"5","combined_ess":"200","split_rhat":"1.01"},
        {"node_id":"n5","descendant_tip_count":"2","descendant_tips":"Actinidia_arguta_A,Actinidia_chinensis","mean_ma":"2","q025_ma":"1","q975_ma":"3","chain1_mean_ma":"2","chain2_mean_ma":"2","combined_ess":"180","split_rhat":"1.02"},
    ])
    write_tsv(paths["calibrations.tsv"], MODULE.CALIBRATION_COLUMNS, [
        {"constraint_id":"root","node_label":"root","node_id":"n4","minimum_ma":"4","maximum_ma":"6","posterior_mean_ma":"5","posterior_q025_ma":"4","posterior_q975_ma":"6","mean_inside_secondary_interval":"true"}
    ])
    write_tsv(paths["model.tsv"], MODULE.MODEL_COLUMNS, [
        {"model_id":"base_poisson","role":"primary single-rate birth-death model","family_count":"10","significant_family_count_p_lt_0.05":"2","score":"100","lambda_values":"0.01","result_file_count":"9"}
    ])
    write_tsv(paths["clades.tsv"], MODULE.CLADE_COLUMNS, [
        {"model_id":"base_poisson","taxon_or_node_id":"Actinidia_arguta_A<1>","increase":"3","decrease":"4"},
        {"model_id":"base_poisson","taxon_or_node_id":"Actinidia_chinensis<2>","increase":"5","decrease":"6"},
        {"model_id":"base_poisson","taxon_or_node_id":"Vitis_vinifera<3>","increase":"7","decrease":"8"},
        {"model_id":"base_poisson","taxon_or_node_id":"<4>","increase":"1","decrease":"1"},
    ])
    paths["cafe_validation.json"].write_text(json.dumps({
        "status":"PASS_CAFE5_BASE_VALIDATED_GAMMA_UNAVAILABLE",
        "calibration_claim":"TimeTree secondary-calibrated; not fossil-calibrated",
        "unavailable_sensitivity":{"status":"UNAVAILABLE_INITIALIZATION_FAILURE"},
    }))
    return paths


class PhylogenyCafeRendererTest(unittest.TestCase):
    def test_species_tree_labels_keep_only_zhejiangensis_lineage_suffixes(self) -> None:
        for species, expected, forbidden in (
            ("Actinidia macrosperma", r"$\mathit{A.\ macrosperma}$", "unphased"),
            ("Actinidia rufa", r"$\mathit{A.\ rufa}$", "unphased"),
            (
                "Actinidia x zhejiangensis parental lineage A",
                r"$\mathit{A.\ zhejiangensis}$ $\mathrm{A}$",
                "parental",
            ),
        ):
            observed = MODULE._species_tree_display_label({"biological_species": species})
            self.assertEqual(observed, expected)
            self.assertNotIn(forbidden, observed)

    def test_publication_rotation_places_outgroups_together_at_bottom(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            tree = Path(temporary) / "tree.tre"
            tree.write_text(
                "((Coffea_arabica_E:4,(Rhododendron_simsii:3,"
                "(Clematoclethra_scandens:2,Actinidia_eriantha_HAP1:2):1):1):1,"
                "Vitis_vinifera:5);\n"
            )
            root = MODULE.parse_newick(tree)
            clades_before = {
                MODULE.descendant_tips(node)
                for node in MODULE.walk(root)
                if not node.is_tip
            }
            MODULE.orient_publication_tree(root)
            tips = [node.label for node in MODULE.walk(root) if node.is_tip]
            clades_after = {
                MODULE.descendant_tips(node)
                for node in MODULE.walk(root)
                if not node.is_tip
            }
        self.assertEqual(
            tips,
            [
                "Actinidia_eriantha_HAP1",
                "Clematoclethra_scandens",
                "Rhododendron_simsii",
                "Coffea_arabica_E",
                "Vitis_vinifera",
            ],
        )
        self.assertEqual(clades_before, clades_after)

    def test_exact_closure_and_typography(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = inputs(Path(temporary))
            _root, rows, validation, calibrated = MODULE.prepare(
                terminals_path=paths["terminals.tsv"], dated_tree_path=paths["tree.tre"],
                node_ages_path=paths["node_ages.tsv"], calibrations_path=paths["calibrations.tsv"],
                model_summary_path=paths["model.tsv"], clade_summary_path=paths["clades.tsv"],
                cafe_node_tree_path=paths["node_tree.tre"],
                cafe_validation_path=paths["cafe_validation.json"], expected_tip_count=3,
                expected_secondary_count=1,
            )
        self.assertEqual(validation["status"], "PASS_PHYLOGENY_CAFE_PUBLICATION_BUNDLE")
        self.assertEqual(validation["maximum_root_to_tip_deviation_ma"], 0.0)
        self.assertEqual(len(rows), 4)
        terminal_rows = [row for row in rows if row["node_type"] == "terminal"]
        internal_rows = [row for row in rows if row["node_type"] == "internal"]
        self.assertEqual(len(terminal_rows), 3)
        self.assertEqual(len(internal_rows), 1)
        self.assertEqual(validation["cafe_branch_count"], 4)
        self.assertEqual(validation["cafe_internal_branch_count"], 1)
        self.assertEqual(terminal_rows[0]["upright_suffix"], "")
        self.assertIn(r"$\mathit{A.\ arguta}$", terminal_rows[0]["display_label"])
        self.assertNotIn(r"$\mathrm{A}$", terminal_rows[0]["display_label"])
        self.assertEqual(len(calibrated), 1)

    def test_tip_and_cafe_closure_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = inputs(Path(temporary))
            paths["tree.tre"].write_text("((Actinidia_arguta_A:2,Wrong_tip:2):3,Vitis_vinifera:5);\n")
            with self.assertRaisesRegex(MODULE.PhylogenyCafeError, "tip set differs"):
                MODULE.prepare(
                    terminals_path=paths["terminals.tsv"], dated_tree_path=paths["tree.tre"],
                    node_ages_path=paths["node_ages.tsv"], calibrations_path=paths["calibrations.tsv"],
                    model_summary_path=paths["model.tsv"], clade_summary_path=paths["clades.tsv"],
                    cafe_node_tree_path=paths["node_tree.tre"],
                    cafe_validation_path=paths["cafe_validation.json"], expected_tip_count=3,
                    expected_secondary_count=1,
                )

    def test_cafe_node_tree_clade_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = inputs(Path(temporary))
            paths["node_tree.tre"].write_text(
                "((Actinidia_arguta_A<1>,Vitis_vinifera<3>)<4>,"
                "Actinidia_chinensis<2>)<5>;\n"
            )
            with self.assertRaisesRegex(
                MODULE.PhylogenyCafeError,
                "descendant clades do not exactly match",
            ):
                MODULE.prepare(
                    terminals_path=paths["terminals.tsv"],
                    dated_tree_path=paths["tree.tre"],
                    node_ages_path=paths["node_ages.tsv"],
                    calibrations_path=paths["calibrations.tsv"],
                    model_summary_path=paths["model.tsv"],
                    clade_summary_path=paths["clades.tsv"],
                    cafe_node_tree_path=paths["node_tree.tre"],
                    cafe_validation_path=paths["cafe_validation.json"],
                    expected_tip_count=3,
                    expected_secondary_count=1,
                )

    @unittest.skipUnless(matplotlib is not None, "optional matplotlib is not installed")
    def test_render_complete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = inputs(root)
            bundle = MODULE.render_bundle(
                terminals_path=paths["terminals.tsv"], dated_tree_path=paths["tree.tre"],
                node_ages_path=paths["node_ages.tsv"], calibrations_path=paths["calibrations.tsv"],
                model_summary_path=paths["model.tsv"], clade_summary_path=paths["clades.tsv"],
                cafe_node_tree_path=paths["node_tree.tre"],
                cafe_validation_path=paths["cafe_validation.json"], output_dir=root / "out",
                basename="phylogeny_cafe", expected_tip_count=3, expected_secondary_count=1, dpi=90,
            )
            self.assertTrue(bundle.png.read_bytes().startswith(b"\x89PNG"))
            self.assertTrue(bundle.pdf.read_bytes().startswith(b"%PDF"))
            self.assertIn("not fossil calibrations", bundle.caption.read_text())
            self.assertIn("all 4 non-root branches", bundle.caption.read_text())


if __name__ == "__main__":
    unittest.main()
