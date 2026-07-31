"""Focused tests for deterministic BUSCO v5 short-summary parsing."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "qc" / "collect_busco.py"
SPEC = importlib.util.spec_from_file_location("collect_busco_test", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class CollectBuscoCountTest(unittest.TestCase):
    def parse_with_count_lines(self, *count_lines: str):
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = (
                Path(temporary_directory)
                / "short_summary.specific.embryophyta_odb10.act_test.txt"
            )
            path.write_text(
                "\n".join(
                    (
                        "# BUSCO version is: 5.8.2",
                        "# The lineage dataset is: embryophyta_odb10 "
                        "(Creation date: 2024-01-08, number of BUSCOs: 1614)",
                        "# Summarized benchmarking in BUSCO notation for file /tmp/input.fa",
                        "# BUSCO was run in mode: euk_genome_min",
                        "C:97.3%[S:70.0%,D:27.3%],F:1.6%,M:1.1%,n:1614",
                        *count_lines,
                        "",
                    )
                ),
                encoding="utf-8",
            )
            return MODULE.parse_short_summary(path)

    def test_complete_count_accepts_busco_internal_stop_suffix(self) -> None:
        summary = self.parse_with_count_lines(
            "1570 Complete BUSCOs (C) (of which 54 contain internal stop codons)",
            "1129 Complete and single-copy BUSCOs (S)",
            "441 Complete and duplicated BUSCOs (D)",
            "26 Fragmented BUSCOs (F)",
            "18 Missing BUSCOs (M)",
        )
        self.assertEqual(
            (summary.C_count, summary.S_count, summary.D_count, summary.F_count, summary.M_count),
            ("1570", "1129", "441", "26", "18"),
        )

    def test_arbitrary_suffix_is_not_accepted(self) -> None:
        summary = self.parse_with_count_lines("1570 Complete BUSCOs (C) unreviewed text")
        self.assertEqual(summary.C_count, "")

    def test_internal_stop_suffix_is_not_accepted_for_noncomplete_row(self) -> None:
        summary = self.parse_with_count_lines(
            "1129 Complete and single-copy BUSCOs (S) "
            "(of which 54 contain internal stop codons)"
        )
        self.assertEqual(summary.S_count, "")


if __name__ == "__main__":
    unittest.main()
