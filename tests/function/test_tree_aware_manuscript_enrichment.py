from __future__ import annotations

import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/function/run_tree_aware_manuscript_enrichment.py"
SPEC = importlib.util.spec_from_file_location("tree_aware_enrichment", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TreeAwareManuscriptEnrichmentTests(unittest.TestCase):
    def test_foreground_membership_and_metadata_counts_close(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genes = root / "genes.tsv.gz"
            with gzip.open(genes, "wt", encoding="utf-8") as handle:
                handle.write("foreground_id\treference_gene_id\nf1\tg1\nf1\tg2\n")
            metadata = root / "metadata.tsv"
            metadata.write_text(
                "foreground_id\tanalysis_scope\tbackground_scope\tbranch_id\t"
                "descendant_lineage_count\tdescendant_lineages\tforeground_gene_count\n"
                "f1\ttree\tall_reference_genes\tb1\t1\tActinidia alpha\t2\n",
                encoding="utf-8",
            )
            foregrounds, rows = MODULE.read_foregrounds(genes, metadata)
            self.assertEqual(foregrounds, {"f1": {"g1", "g2"}})
            self.assertEqual(rows["f1"]["branch_id"], "b1")

    def test_bh_helper_is_the_frozen_enrichment_implementation(self) -> None:
        self.assertEqual(MODULE.BASE.bh_adjust([0.01, 0.2]), [0.02, 0.2])

    def test_custom_backgrounds_remain_scope_specific(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "backgrounds.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(
                    "background_scope\treference_gene_id\n"
                    "resolved__u1\tg1\n"
                    "resolved__u1\tg2\n"
                    "resolved__u2\tg2\n"
                )
            self.assertEqual(
                MODULE.read_backgrounds(path),
                {
                    "resolved__u1": {"g1", "g2"},
                    "resolved__u2": {"g2"},
                },
            )


if __name__ == "__main__":
    unittest.main()
