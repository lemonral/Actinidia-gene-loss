from __future__ import annotations

import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "function" / "prepare_species_specific_loss_foregrounds.py"
SPEC = importlib.util.spec_from_file_location("species_specific_foregrounds", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SpeciesSpecificForegroundTests(unittest.TestCase):
    def test_prepare_keeps_only_exact_single_terminal_patterns(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            patterns = root / "patterns.tsv.gz"
            with gzip.open(patterns, "wt", encoding="utf-8") as handle:
                handle.write(
                    "reference_gene_id\ttree_pattern\ttree_placement_exact\t"
                    "complete_loss_lineage_count\tpartial_loss_lineage_count\t"
                    "unknown_lineage_count\tminimum_loss_event_count\n"
                    "g1\tsingle_terminal_branch_loss\ttrue\t1\t0\t0\t1\n"
                    "g2\tsingle_terminal_branch_loss\ttrue\t1\t0\t0\t1\n"
                    "g3\trecurrent_independent_losses\ttrue\t2\t0\t0\t2\n"
                )
            events = root / "events.tsv.gz"
            with gzip.open(events, "wt", encoding="utf-8") as handle:
                handle.write(
                    "reference_gene_id\tbranch_id\tbranch_type\t"
                    "descendant_lineage_count\tdescendant_lineages\t"
                    "gene_minimum_loss_event_count\n"
                    "g1\tterminal__a\tterminal\t1\tActinidia alpha\t1\n"
                    "g2\tterminal__b\tterminal\t1\tActinidia beta\t1\n"
                    "g3\tterminal__a\tterminal\t1\tActinidia alpha\t2\n"
                )
            membership, metadata, summary = MODULE.prepare(
                patterns,
                events,
                expected_lineages=2,
                expected_genes=2,
            )
            self.assertEqual(
                {row["reference_gene_id"] for row in membership},
                {"g1", "g2"},
            )
            self.assertEqual(len(metadata), 2)
            self.assertEqual(
                sum(row["foreground_gene_count"] for row in metadata),
                2,
            )
            self.assertEqual(
                summary["status"],
                "PASS_ARTICLE_METHOD_SINGLE_TERMINAL_FOREGROUNDS",
            )


if __name__ == "__main__":
    unittest.main()
