from __future__ import annotations

import tempfile
import unittest
from decimal import Decimal
from pathlib import Path

from geneloss_repro.nucleotide_matrix import (
    NucleotideMatrixError,
    build_nucleotide_rows,
    read_role_normalized_paf,
)


def paf_row(
    query: str,
    query_length: int,
    query_start: int,
    query_end: int,
    strand: str,
    target: str,
    target_length: int,
    target_start: int,
    target_end: int,
    *,
    matching: int = 18_000,
    block: int = 20_000,
    mapq: int = 60,
    de: str = "0.05",
    tp: str = "P",
) -> str:
    return "\t".join(
        (
            query,
            str(query_length),
            str(query_start),
            str(query_end),
            strand,
            target,
            str(target_length),
            str(target_start),
            str(target_end),
            str(matching),
            str(block),
            str(mapq),
            f"tp:A:{tp}",
            f"de:f:{de}",
            f"cg:Z:{block}M",
            f"cs:Z::{block}",
        )
    )


class NucleotideMatrixTests(unittest.TestCase):
    def test_bidirectional_interval_unions_are_not_raw_row_sums(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            forward = root / "forward.paf"
            reverse = root / "reverse.paf"
            forward.write_text(
                paf_row("PubChr01", 2_000_000, 0, 20_000, "+", "Chr01A", 2_100_000, 0, 20_000)
                + "\n"
                + paf_row(
                    "PubChr01",
                    2_000_000,
                    10_000,
                    30_000,
                    "+",
                    "Chr01A",
                    2_100_000,
                    10_000,
                    30_000,
                )
                + "\n",
                encoding="utf-8",
            )
            reverse.write_text(
                paf_row(
                    "Chr01A",
                    2_100_000,
                    20_000,
                    40_000,
                    "+",
                    "PubChr01",
                    2_000_000,
                    20_000,
                    40_000,
                )
                + "\n",
                encoding="utf-8",
            )
            target_lengths = {"PubChr01": 2_000_000}
            reference_lengths = {"Chr01A": 2_100_000}
            first, first_audit = read_role_normalized_paf(
                forward,
                query_role="target",
                target_lengths=target_lengths,
                reference_lengths=reference_lengths,
            )
            second, second_audit = read_role_normalized_paf(
                reverse,
                query_role="reference",
                target_lengths=target_lengths,
                reference_lengths=reference_lengths,
            )
            rows, orientation = build_nucleotide_rows(
                target_lengths=target_lengths,
                reference_lengths=reference_lengths,
                canonical_by_reference={"Chr01A": "Chr01"},
                forward_records=first,
                reverse_records=second,
            )
            self.assertEqual(first_audit.retained_rows, 2)
            self.assertEqual(second_audit.retained_rows, 1)
            self.assertEqual(rows[0]["query_covered_bp"], "40000")
            self.assertEqual(rows[0]["reference_covered_bp"], "40000")
            self.assertEqual(rows[0]["matching_bases"], "54000")
            self.assertEqual(rows[0]["orientation"], "mixed")
            self.assertEqual(orientation[0]["dominant_orientation"], "+")
            self.assertEqual(orientation[0]["automatic_orientation_gate"], "false")

    def test_orientation_gate_uses_dominant_matching_support(self) -> None:
        from geneloss_repro.nucleotide_matrix import PafRecord

        records = [
            PafRecord(
                "PubChr01",
                2_000_000,
                0,
                1_200_000,
                "Chr01A",
                2_100_000,
                0,
                1_200_000,
                "-",
                1_100_000,
                1_200_000,
                Decimal("0.05"),
            ),
            PafRecord(
                "PubChr01",
                2_000_000,
                1_200_000,
                1_300_000,
                "Chr01A",
                2_100_000,
                1_200_000,
                1_300_000,
                "+",
                100_000,
                100_000,
                Decimal("0.05"),
            ),
        ]
        rows, audit = build_nucleotide_rows(
            target_lengths={"PubChr01": 2_000_000},
            reference_lengths={"Chr01A": 2_100_000},
            canonical_by_reference={"Chr01A": "Chr01"},
            forward_records=records,
            reverse_records=[],
        )
        self.assertEqual(rows[0]["orientation"], "-")
        self.assertEqual(audit[0]["dominant_fraction"], "0.916666666666667")
        self.assertEqual(audit[0]["automatic_orientation_gate"], "true")

    def test_unaligned_cartesian_cell_is_explicit(self) -> None:
        rows, _ = build_nucleotide_rows(
            target_lengths={"PubChr01": 100_000, "PubChr02": 100_000},
            reference_lengths={"Chr01A": 100_000, "Chr02A": 100_000},
            canonical_by_reference={"Chr01A": "Chr01", "Chr02A": "Chr02"},
            forward_records=[],
            reverse_records=[],
        )
        self.assertEqual(len(rows), 4)
        self.assertEqual({row["score"] for row in rows}, {"0"})
        self.assertEqual({row["weighted_divergence"] for row in rows}, {"1"})
        self.assertEqual({row["orientation"] for row in rows}, {"none"})

    def test_missing_required_paf_tag_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "broken.paf"
            path.write_text(
                "\t".join(
                    paf_row(
                        "PubChr01", 20_000, 0, 20_000, "+", "Chr01A", 20_000, 0, 20_000
                    ).split("\t")[:-1]
                )
                + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(NucleotideMatrixError, "missing required PAF tags"):
                read_role_normalized_paf(
                    path,
                    query_role="target",
                    target_lengths={"PubChr01": 20_000},
                    reference_lengths={"Chr01A": 20_000},
                )


if __name__ == "__main__":
    unittest.main()
