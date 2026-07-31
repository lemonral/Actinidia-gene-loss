from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "downstream" / "prepare_manuscript_method_downstream.py"


class ManuscriptMethodDownstreamTests(unittest.TestCase):
    def test_unit_shared_partial_and_tree_events_are_separate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ledger = root / "units.tsv"
            ledger.write_text(
                "sample_id\tspecies\tploidy\tanalysis_role\tinput_scope\n"
                "ref\tReference species\tn/a\treference_callable\treference_transcript_cds\n"
                "a1\tActinidia alpha\t2x\ttarget_repertoire\twhole_genome\n"
                "a2\tActinidia alpha\t2x\ttarget_repertoire\twhole_genome\n"
                "b1\tActinidia beta\t4x\ttarget_repertoire\twhole_genome\n"
                "c1\tActinidia gamma\t6x\ttarget_repertoire\twhole_genome\n",
                encoding="utf-8",
            )
            species_map = root / "species.tsv"
            species_map.write_text(
                "assembly_unit_id\tbiological_species\tinclude\n"
                "a1\tActinidia alpha\ttrue\n"
                "a2\tActinidia alpha\ttrue\n"
                "b1\tActinidia beta\ttrue\n"
                "c1\tActinidia gamma\ttrue\n",
                encoding="utf-8",
            )
            tip_map = root / "tips.tsv"
            tip_map.write_text(
                "tree_tip\tbiological_species\tinclude\n"
                "A\tActinidia alpha\ttrue\n"
                "B\tActinidia beta\ttrue\n"
                "C\tActinidia gamma\ttrue\n"
                "O\tReference species\tfalse\n",
                encoding="utf-8",
            )
            tree = root / "tree.tre"
            tree.write_text("(((A:1,B:1):1,C:2):1,O:3);\n", encoding="utf-8")
            states = {
                "g1": {"a1": "decayed", "a2": "deleted", "b1": "decayed", "c1": "deleted"},
                "g2": {"a1": "decayed", "a2": "deleted", "b1": "retained", "c1": "retained"},
                "g3": {"a1": "decayed", "a2": "retained", "b1": "retained", "c1": "retained"},
                "g4": {"a1": "deleted", "a2": "decayed", "b1": "decayed", "c1": "retained"},
                "g5": {"a1": "deleted", "a2": "decayed", "b1": "retained", "c1": "decayed"},
            }
            matrix = root / "matrix.tsv.gz"
            with gzip.open(matrix, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "reference_gene_id", "assembly_unit_id", "manuscript_classification",
                        "manuscript_positive_loss", "manuscript_rule", "refined_decayed_cause",
                    ],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for gene, unit_states in states.items():
                    for unit, state in unit_states.items():
                        writer.writerow(
                            {
                                "reference_gene_id": gene,
                                "assembly_unit_id": unit,
                                "manuscript_classification": state,
                                "manuscript_positive_loss": str(state in {"decayed", "deleted"}).lower(),
                                "manuscript_rule": "article_threshold",
                                "refined_decayed_cause": "test",
                            }
                        )
            expression = root / "expression.tsv"
            expression.write_text(
                "reference_gene_id\tleaf_raw_count\n"
                "g1\t1\ng2\t2\ng3\t3\ng4\t4\ng5\t5\n",
                encoding="utf-8",
            )
            clusters = root / "reference.clstr"
            clusters.write_text(
                "".join(
                    f">Cluster {index}\n0 3aa, >{gene}... *\n"
                    for index, gene in enumerate(states)
                ),
                encoding="utf-8",
            )
            output = root / "out"
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--manuscript-matrix", str(matrix),
                    "--unit-ledger", str(ledger),
                    "--species-map", str(species_map),
                    "--tip-map", str(tip_map),
                    "--time-tree", str(tree),
                    "--reference-expression", str(expression),
                    "--clusters", str(clusters),
                    "--output-dir", str(output),
                    "--expected-units", "4",
                    "--expected-reference-genes", "5",
                    "--expected-lineages", "3",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            shared = (output / "shared_23_unit_genes.tsv").read_text(encoding="utf-8")
            self.assertIn("g1", shared)
            with gzip.open(output / "gene_tree_patterns.tsv.gz", "rt", encoding="utf-8") as handle:
                patterns = {row["reference_gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(patterns["g2"]["tree_pattern"], "single_terminal_branch_loss")
            self.assertEqual(patterns["g3"]["tree_pattern"], "partial_lineage_loss_only")
            self.assertEqual(patterns["g4"]["tree_pattern"], "single_internal_branch_loss")
            self.assertEqual(patterns["g5"]["tree_pattern"], "recurrent_independent_losses")
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["counts"]["shared_23_unit_genes"], 1)
            self.assertEqual(manifest["status"], "PASS_MANUSCRIPT_METHOD_DOWNSTREAM_AND_TREE_PATTERNS")


if __name__ == "__main__":
    unittest.main()
