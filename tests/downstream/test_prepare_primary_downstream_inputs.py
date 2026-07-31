from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "downstream" / "prepare_primary_downstream_inputs.py"
SPEC = importlib.util.spec_from_file_location("prepare_primary_downstream_inputs", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class PreparePrimaryDownstreamInputsTests(unittest.TestCase):
    def test_builds_resolved_expression_and_missing_cluster_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cds = root / "reference.fa"
            cds.write_text(">g1\nAAA\n>g2\nCCC\n>g3\nGGG\n", encoding="utf-8")
            ledger = root / "units.tsv"
            with ledger.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(MODULE.LEDGER_FIELDS)
                writer.writerow(["ref", "Ref species", "n/a", "reference_callable", "reference_transcript_cds", "ref.fa", "ref.fa", "3"])
                writer.writerow(["u1", "Actinidia alpha", "2x", "target_repertoire", "whole_genome", "u1.fa", "u1.fa", "29"])
                writer.writerow(["u2", "Actinidia beta", "4x", "target_repertoire", "whole_genome", "u2.fa", "u2.fa", "29"])
            loss = root / "loss.tsv"
            rows = []
            classes = {"u1": ["deleted", "retained", "uncertain"], "u2": ["retained", "deleted", "retained"]}
            for unit, states in classes.items():
                for gene, state in zip(["g1", "g2", "g3"], states):
                    rows.append([gene, unit, state, "true", "test", state])
            with loss.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(MODULE.LOSS_FIELDS)
                writer.writerows(rows)
            shared = root / "shared.tsv"
            shared.write_text("reference_gene_id\ng3\n", encoding="utf-8")
            counts = root / "counts.tsv"
            counts.write_text(
                "# featureCounts\nGeneid\tLength\tleaf.bam\n"
                "g1\t3\t1\ng2\t3\t0\ng3\t3\t7\nextra\t3\t9\n",
                encoding="utf-8",
            )
            clusters = root / "reference.clstr"
            clusters.write_text(
                ">Cluster 0\n0 3aa, >g1... *\n",
                encoding="utf-8",
            )
            output = root / "out"
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT),
                    "--loss-matrix", str(loss),
                    "--unit-ledger", str(ledger),
                    "--reference-cds", str(cds),
                    "--shared-positive-genes", str(shared),
                    "--feature-counts", str(counts),
                    "--expression-sample-column", "leaf.bam",
                    "--clusters", str(clusters),
                    "--output-dir", str(output),
                    "--expected-units", "2",
                    "--expected-reference-genes", "3",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with (output / "resolved_nonshared_unit_loss_table.tsv").open(encoding="utf-8") as handle:
                resolved = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(resolved), 4)
            self.assertEqual({row["classification"] for row in resolved}, {"deleted", "retained"})
            self.assertEqual({row["ploidy"] for row in resolved}, {"diploid", "tetraploid"})
            with (output / "nonshared_reference_leaf_raw_counts.tsv").open(encoding="utf-8") as handle:
                expression = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["leaf_raw_count"] for row in expression], ["1", "0"])
            self.assertEqual(
                (output / "nonshared_cdhit_missing_reference_ids.txt").read_text(),
                "g2\n",
            )

    def test_rejects_unexpected_classification(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cds = root / "reference.fa"
            cds.write_text(">g1\nAAA\n", encoding="utf-8")
            ledger = root / "units.tsv"
            ledger.write_text(
                "\t".join(MODULE.LEDGER_FIELDS) + "\n"
                "ref\tRef\tn/a\treference_callable\treference_transcript_cds\tref.fa\tref.fa\t1\n"
                "u1\tActinidia alpha\t2x\ttarget_repertoire\twhole_genome\tu1.fa\tu1.fa\t29\n",
                encoding="utf-8",
            )
            loss = root / "loss.tsv"
            loss.write_text("\t".join(MODULE.LOSS_FIELDS) + "\ng1\tu1\tlost\ttrue\ttest\tlost\n", encoding="utf-8")
            shared = root / "shared.tsv"
            shared.write_text("reference_gene_id\ng1\n", encoding="utf-8")
            counts = root / "counts.tsv"
            counts.write_text("Geneid\tleaf\ng1\t1\n", encoding="utf-8")
            clusters = root / "ref.clstr"
            clusters.write_text(">Cluster 0\n0 3aa, >g1... *\n", encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--loss-matrix", str(loss),
                    "--unit-ledger", str(ledger), "--reference-cds", str(cds),
                    "--shared-positive-genes", str(shared),
                    "--feature-counts", str(counts), "--expression-sample-column", "leaf",
                    "--clusters", str(clusters), "--output-dir", str(root / "out"),
                    "--expected-units", "1", "--expected-reference-genes", "1",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("classification", result.stderr)


if __name__ == "__main__":
    unittest.main()
