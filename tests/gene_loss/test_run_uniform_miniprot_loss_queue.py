"""Tests for unified Miniprot loss/disruption classification."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "gene_loss" / "run_uniform_miniprot_loss_queue.py"
SPEC = importlib.util.spec_from_file_location("run_uniform_miniprot_loss_queue", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def candidate(callable_value: str = "true") -> dict[str, str]:
    return {
        "unit": "u1",
        "reference_gene": "g1",
        "callable": callable_value,
        "callability_reason": "callable" if callable_value == "true" else "missing_bilateral_anchor",
        "target_chromosome": "chr1" if callable_value == "true" else "",
        "target_interval_start_1based": "100" if callable_value == "true" else "",
        "target_interval_end_1based": "1000" if callable_value == "true" else "",
    }


def alignment(*, frameshifts: int = 0, stops: int = 0) -> dict[str, object]:
    return {
        "query": "g1",
        "query_length": 100,
        "query_start": 0,
        "query_end": 90,
        "strand": "+",
        "target": "chr1",
        "target_start_1based": 200,
        "target_end_1based": 800,
        "query_coverage": 0.9,
        "identity": 0.8,
        "alignment_score": 200,
        "frameshifts": frameshifts,
        "stops": stops,
        "mapq": 60,
    }


class UniformMiniprotTests(unittest.TestCase):
    def classify(self, value: dict[str, object] | None, callable_value: str = "true") -> dict[str, str]:
        return MODULE.classify_candidate(
            candidate(callable_value),
            value,
            minimum_query_coverage=0.5,
            minimum_identity=0.5,
            minimum_alignment_score=50,
        )

    def test_paf_parser_reads_disruptive_tags(self) -> None:
        raw = "g1\t100\t0\t90\t+\tchr1\t2000\t199\t800\t210\t270\t60\tAS:i:200\tfs:i:1\tst:i:2\n"
        parsed = MODULE.parse_paf_line(raw, 1)
        self.assertEqual(parsed["target_start_1based"], 200)
        self.assertEqual(parsed["frameshifts"], 1)
        self.assertEqual(parsed["stops"], 2)

    def test_no_local_alignment_is_deleted(self) -> None:
        self.assertEqual(self.classify(None)["classification"], "deleted")

    def test_frameshift_is_positive_pseudogenized(self) -> None:
        row = self.classify(alignment(frameshifts=1))
        self.assertEqual(row["classification"], "pseudogenized")
        self.assertEqual(row["positive_loss"], "true")

    def test_internal_stop_is_positive_pseudogenized(self) -> None:
        self.assertEqual(self.classify(alignment(stops=1))["classification"], "pseudogenized")

    def test_sequence_without_disruption_remains_uncertain(self) -> None:
        self.assertEqual(self.classify(alignment())["classification"], "uncertain")

    def test_non_callable_remains_uncertain(self) -> None:
        row = self.classify(None, callable_value="false")
        self.assertEqual(row["classification"], "uncertain")
        self.assertEqual(row["positive_loss"], "false")

    def test_best_alignment_must_be_contained_in_expected_interval(self) -> None:
        inside = alignment()
        outside = {**alignment(), "alignment_score": 999, "target_start_1based": 50}
        self.assertEqual(MODULE.best_local_alignment(candidate(), [outside, inside]), inside)


if __name__ == "__main__":
    unittest.main()
