"""Regression tests for strand-aware GFF3 CDS phase handling."""

from __future__ import annotations

import unittest

from geneloss_repro.annotation import build_spliced_cds, translate_standard
from geneloss_repro.gff import GffFeature, Transcript
from geneloss_repro.io_utils import SchemaError


def feature(sequence_id: str, start: int, end: int, strand: str, phase: str) -> GffFeature:
    return GffFeature(
        sequence_id=sequence_id,
        feature_type="CDS",
        start=start,
        end=end,
        strand=strand,
        phase=phase,
        attributes={"Parent": "tx1"},
        line_number=1,
    )


class AnnotationPhaseTest(unittest.TestCase):
    def test_internal_phase_bases_are_retained_on_plus_strand(self) -> None:
        genome = {"Chr1": "ATGGNNNNNCCAAA"}
        transcript = Transcript(
            "gene1",
            "tx1",
            "Chr1",
            1,
            14,
            "+",
            (feature("Chr1", 1, 4, "+", "0"), feature("Chr1", 10, 14, "+", "2")),
        )
        cds = build_spliced_cds(transcript, genome)
        self.assertEqual(cds, "ATGGCCAAA")
        self.assertEqual(translate_standard(cds), "MAK")

    def test_internal_phase_bases_are_retained_on_minus_strand(self) -> None:
        genome = {"Chr1": "TTTGGNNNNCCAT"}
        transcript = Transcript(
            "gene1",
            "tx1",
            "Chr1",
            1,
            13,
            "-",
            (feature("Chr1", 1, 5, "-", "2"), feature("Chr1", 10, 13, "-", "0")),
        )
        cds = build_spliced_cds(transcript, genome)
        self.assertEqual(cds, "ATGGCCAAA")
        self.assertEqual(translate_standard(cds), "MAK")

    def test_initial_phase_is_applied_once_for_partial_cds(self) -> None:
        genome = {"Chr1": "AATGAAA"}
        transcript = Transcript(
            "gene1",
            "tx1",
            "Chr1",
            1,
            7,
            "+",
            (feature("Chr1", 1, 7, "+", "1"),),
        )
        cds = build_spliced_cds(transcript, genome)
        self.assertEqual(cds, "ATGAAA")
        self.assertEqual(translate_standard(cds), "MK")

    def test_inconsistent_phase_chain_fails_closed(self) -> None:
        genome = {"Chr1": "ATGGNNNNNCCAAA"}
        transcript = Transcript(
            "gene1",
            "tx1",
            "Chr1",
            1,
            14,
            "+",
            (feature("Chr1", 1, 4, "+", "0"), feature("Chr1", 10, 14, "+", "1")),
        )
        with self.assertRaisesRegex(SchemaError, "inconsistent CDS phase chain"):
            build_spliced_cds(transcript, genome)


if __name__ == "__main__":
    unittest.main()
