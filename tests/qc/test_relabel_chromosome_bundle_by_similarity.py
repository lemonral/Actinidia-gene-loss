"""Contracts for the simplified HY4A chromosome naming workflow."""

from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "qc" / "relabel_chromosome_bundle_by_similarity.py"
SPEC = importlib.util.spec_from_file_location("relabel_chromosome_bundle_by_similarity", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class RelabelTests(unittest.TestCase):
    def test_map_requires_complete_chr_bijection(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "map.tsv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=MODULE.LABEL_COLUMNS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for index in range(1, 30):
                    writer.writerow(
                        {
                            "query_chromosome": f"p{index}",
                            "final_chromosome": f"Chr{index:02d}",
                            "coordinate_reference": "act_chinensis_hongyang_v4_hy4a",
                            "assignment_method": "global_one_to_one_maximum_nucleotide_similarity",
                            "assigned_score": "1",
                            "reciprocal_coverage": "0.01",
                            "orientation_to_hy4a": "+",
                            "hy4p_and_jcvi_agree": "true",
                            "strict_homology_gates_pass": "false",
                            "confidence_flag": "SUPPORTED",
                        }
                    )
            labels, confidence = MODULE.read_label_map(path)
            self.assertEqual(labels["p1"], "Chr01")
            self.assertEqual(confidence["p29"], "SUPPORTED")


if __name__ == "__main__":
    unittest.main()
