from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from geneloss_repro.chromosome_harmonization import (
    HarmonizationError,
    build_actions,
    read_fasta,
    transform_genome,
    transform_gff,
    validate_cds_protein_closure,
    validate_sequence_closure,
    write_fasta,
)


class ChromosomeHarmonizationTests(unittest.TestCase):
    def test_relabel_flip_and_sequence_closure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome_path = root / "source.fa.gz"
            with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
                handle.write(">Pub01\nATGAAATAG\n>Pub02\nCTATTTCAT\n")
            gff = root / "source.gff3"
            gff.write_text(
                "##gff-version 3\n"
                "##sequence-region Pub01 1 9\n"
                "##sequence-region Pub02 1 9\n"
                "Pub01\tx\tgene\t1\t9\t.\t+\t.\tID=g1\n"
                "Pub01\tx\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
                "Pub01\tx\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1\n"
                "Pub02\tx\tgene\t1\t9\t.\t-\t.\tID=g2\n"
                "Pub02\tx\tmRNA\t1\t9\t.\t-\t.\tID=t2;Parent=g2\n"
                "Pub02\tx\tCDS\t1\t9\t.\t-\t0\tID=c2;Parent=t2\n",
                encoding="utf-8",
            )
            genome = read_fasta(genome_path)
            actions = build_actions(
                genome,
                {"Pub01": "Chr02", "Pub02": "Chr01"},
                {"Pub01": "+", "Pub02": "-"},
            )
            transformed = transform_genome(genome, actions)
            transformed_fasta = root / "transformed.fa"
            transformed_gff = root / "transformed.gff3"
            write_fasta(transformed_fasta, transformed)
            audit = transform_gff(gff, transformed_gff, actions)
            validate_sequence_closure(
                source_genome=genome,
                transformed_genome=read_fasta(transformed_fasta),
                actions=actions,
            )
            closure = validate_cds_protein_closure(
                source_genome=genome,
                source_gff=gff,
                transformed_genome=read_fasta(transformed_fasta),
                transformed_gff=transformed_gff,
                expected_cds={"t1": "ATGAAATAG", "t2": "ATGAAATAG"},
                expected_proteins={"t1": "MK", "t2": "MK"},
            )
            self.assertEqual(transformed["Chr01"], "ATGAAATAG")
            self.assertEqual(audit["reversed_feature_rows"], 3)
            self.assertEqual(closure["exact_cds_matches"], 2)
            text = transformed_gff.read_text(encoding="utf-8")
            self.assertIn("Chr01\tx\tmRNA\t1\t9\t.\t+\t.\tID=t2;Parent=g2", text)

    def test_mixed_orientation_fails_closed(self) -> None:
        with self.assertRaisesRegex(HarmonizationError, "manual review"):
            build_actions(
                {"Pub01": "A"}, {"Pub01": "Chr01"}, {"Pub01": "mixed"}
            )

    def test_deterministic_gzip_fasta_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first.fa.gz"
            second = root / "second.fa.gz"
            records = {"Chr02": "TTT", "Chr01": "AAA"}
            write_fasta(first, records)
            write_fasta(second, records)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(read_fasta(first), records)


if __name__ == "__main__":
    unittest.main()
