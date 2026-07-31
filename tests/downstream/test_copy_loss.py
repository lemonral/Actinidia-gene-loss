from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "downstream" / "copy_loss.py"


class CopyLossTests(unittest.TestCase):
    def test_minimum_gene_filter_uses_reference_class_not_resolved_sample_size(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            clusters = root / "reference.clstr"
            clusters.write_text(
                ">Cluster 0\n0 3aa, >g1... *\n1 3aa, >g2... at 95%\n"
                ">Cluster 1\n0 3aa, >g3... *\n",
                encoding="utf-8",
            )
            losses = root / "loss.tsv"
            losses.write_text(
                "reference_gene_id\tassembly_unit_id\tploidy\tclassification\n"
                "g1\tu1\tdiploid\tdeleted\n"
                "g2\tu1\tdiploid\tretained\n"
                "g3\tu1\tdiploid\tretained\n"
                "g1\tu2\ttetraploid\tretained\n"
                "g3\tu2\ttetraploid\tdeleted\n",
                encoding="utf-8",
            )
            output = root / "out"
            result = subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--clusters", str(clusters),
                    "--loss-table", str(losses),
                    "--output-dir", str(output),
                    "--gene-column", "reference_gene_id",
                    "--sample-column", "assembly_unit_id",
                    "--ploidy-column", "ploidy",
                    "--class-column", "classification",
                    "--lost-values", "deleted",
                    "--min-genes", "2",
                    "--expected-copy-numbers", "2",
                    "--allow-incomplete-loss-coverage",
                    "--no-plot",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            rates = pd.read_csv(output / "copy_loss_rates_eligible.tsv", sep="\t")
            class_two = rates.loc[rates["copy_number"] == 2]
            self.assertEqual(set(class_two["sample_id"]), {"u1", "u2"})
            self.assertEqual(
                int(class_two.loc[class_two["sample_id"] == "u2", "total_genes"].iloc[0]),
                1,
            )
            self.assertTrue(class_two["passes_min_genes"].all())


if __name__ == "__main__":
    unittest.main()
