"""Unit tests for callable-aware translated-search classification."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "gene_loss" / "run_translated_search_queue.py"
SPEC = importlib.util.spec_from_file_location("run_translated_search_queue", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def hit(subject: str, start: int, end: int, identity: str = "75", bitscore: str = "100") -> dict[str, str]:
    return {
        "qseqid": "g1", "sseqid": subject, "pident": identity, "length": "90",
        "qstart": "1", "qend": "90", "sstart": str(start), "send": str(end),
        "evalue": "1e-20", "bitscore": bitscore, "qframe": "1", "sframe": "1",
    }


def candidate(callable_value: str = "true") -> dict[str, str]:
    return {
        "reference_gene": "g1", "callable": callable_value,
        "callability_reason": "callable" if callable_value == "true" else "missing_bilateral_anchor",
        "target_chromosome": "chr1", "target_interval_start_1based": "100",
        "target_interval_end_1based": "200",
    }


class SearchTests(unittest.TestCase):
    def test_local_hit_is_uncertain_not_pseudogenized(self) -> None:
        rows, counts = MODULE.classify(
            unit="u", candidate_rows=[candidate()], hits={"g1": [hit("chr1", 120, 180)]},
            minimum_identity=50, minimum_bitscore=50, maximum_evalue=1e-5,
        )
        self.assertEqual(rows[0]["primary_state"], "uncertain_local_genomic_sequence_detected")
        self.assertEqual(rows[0]["positive_loss"], "false")
        self.assertEqual(counts["historical_decayed"], 1)

    def test_callable_no_local_hit_is_positive_deleted(self) -> None:
        rows, _ = MODULE.classify(
            unit="u", candidate_rows=[candidate()], hits={"g1": [hit("chr2", 120, 180)]},
            minimum_identity=50, minimum_bitscore=50, maximum_evalue=1e-5,
        )
        self.assertEqual(rows[0]["primary_state"], "positive_deleted")
        self.assertEqual(rows[0]["positive_loss"], "true")

    def test_non_callable_is_never_positive(self) -> None:
        rows, _ = MODULE.classify(
            unit="u", candidate_rows=[candidate("false")], hits={},
            minimum_identity=50, minimum_bitscore=50, maximum_evalue=1e-5,
        )
        self.assertEqual(rows[0]["primary_state"], "uncertain_non_callable")
        self.assertEqual(rows[0]["positive_loss"], "false")

    def test_thresholds_match_verified_historical_final_rule(self) -> None:
        self.assertTrue(MODULE.qualifies(hit("chr1", 1, 2, "50", "50"), identity=50, bitscore=50, evalue=1e-5))
        low = hit("chr1", 1, 2, "49.999", "50")
        self.assertFalse(MODULE.qualifies(low, identity=50, bitscore=50, evalue=1e-5))

    def test_empty_output_can_be_checksum_bound(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.log"
            path.touch()
            self.assertEqual(MODULE.binding(path, allow_empty=True)["bytes"], 0)
            with self.assertRaises(MODULE.SearchError):
                MODULE.binding(path)


if __name__ == "__main__":
    unittest.main()
