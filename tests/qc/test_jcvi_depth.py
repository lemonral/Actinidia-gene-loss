"""Focused tests for generic JCVI gene-depth reconstruction."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
import unittest
from pathlib import Path


QC_SCRIPTS = Path(__file__).resolve().parents[2] / "scripts" / "qc"
sys.path.insert(0, str(QC_SCRIPTS))

from jcvi_depth import ReconstructionError, load_bed, summarize_depth  # noqa: E402
from summarize_jcvi_depth import DepthSummaryError, run  # noqa: E402


class JcviDepthHelpersTest(unittest.TestCase):
    def test_legacy_half_open_endpoint_is_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            bed_path = Path(temporary_directory) / "reference.bed"
            bed_path.write_text(
                "Chr1\t30\t39\tR4\n"
                "Chr1\t0\t9\tR1\n"
                "Chr1\t20\t29\tR3\n"
                "Chr1\t10\t19\tR2\n",
                encoding="utf-8",
            )
            summary = summarize_depth(
                [[("R1", "Q1"), ("R3", "Q3")]], load_bed(bed_path), 0
            )
            self.assertEqual(summary["total"], 4)
            self.assertEqual(summary["nonzero"], 2)
            self.assertEqual(summary["zero"], 2)
            self.assertEqual(summary["coverage"], 50.0)

    def test_cross_sequence_anchor_block_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            bed_path = Path(temporary_directory) / "reference.bed"
            bed_path.write_text(
                "Chr1\t0\t9\tR1\nChr2\t0\t9\tR2\n", encoding="utf-8"
            )
            with self.assertRaises(ReconstructionError):
                summarize_depth(
                    [[("R1", "Q1"), ("R2", "Q2")]], load_bed(bed_path), 0
                )


class JcviDepthIdentifierPolicyTest(unittest.TestCase):
    @staticmethod
    def _write_inputs(root: Path, reference_bed_only: bool) -> argparse.Namespace:
        reference_protein = root / "reference.fa"
        query_protein = root / "query.fa"
        reference_bed = root / "reference.bed"
        query_bed = root / "query.bed"
        anchors = root / "comparison.anchors"
        depthfile = root / "comparison.depth"

        reference_protein.write_text(
            ">R1\nM\n>R2\nM\n>R3\nM\n>R4\nM\n", encoding="utf-8"
        )
        query_protein.write_text(
            ">Q1\nM\n>Q2\nM\n>Q3\nM\n>Q4\nM\n", encoding="utf-8"
        )
        reference_rows = (
            "Chr1\t0\t9\tR1\n"
            "Chr1\t10\t19\tR2\n"
            "Chr1\t20\t29\tR3\n"
            "Chr1\t30\t39\tR4\n"
        )
        if reference_bed_only:
            reference_rows += "Chr1\t40\t49\tRX\n"
        reference_bed.write_text(reference_rows, encoding="utf-8")
        query_bed.write_text(
            "ChrQ\t0\t9\tQ1\n"
            "ChrQ\t10\t19\tQ2\n"
            "ChrQ\t20\t29\tQ3\n"
            "ChrQ\t30\t39\tQ4\n",
            encoding="utf-8",
        )
        anchors.write_text("###\nR1\tQ1\nR3\tQ3\n", encoding="utf-8")
        depth_rows = (
            "R1\t1\nR2\t1\nR3\t0\nR4\t0\n"
            + ("RX\t0\n" if reference_bed_only else "")
            + "Q1\t1\nQ2\t1\nQ3\t0\nQ4\t0\n"
        )
        depthfile.write_text(depth_rows, encoding="utf-8")
        return argparse.Namespace(
            sample="Actinidia_test_HAP1",
            display_name="A. test HAP1",
            accession="TEST0001",
            reference_protein=reference_protein,
            reference_bed=reference_bed,
            query_protein=query_protein,
            query_bed=query_bed,
            anchors=anchors,
            depthfile=depthfile,
            output_tsv=root / "summary.tsv",
            output_json=root / "summary.json",
            minimum_block_size=2,
            allowed_reference_bed_only_ids=None,
        )

    def test_default_requires_exact_reference_bed_fasta_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            args = self._write_inputs(Path(temporary_directory), reference_bed_only=True)
            with self.assertRaisesRegex(DepthSummaryError, "BED/FASTA ID contract"):
                run(args)

    def test_explicit_allowlist_must_exactly_match_bed_only_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            args = self._write_inputs(root, reference_bed_only=True)
            allowlist = root / "reference_bed_only.txt"
            allowlist.write_text("# audited exception\nRX\n", encoding="utf-8")
            args.allowed_reference_bed_only_ids = allowlist
            run(args)

            payload = json.loads(args.output_json.read_text(encoding="utf-8"))
            self.assertEqual(payload["reference_bed_only_ids"], ["RX"])
            self.assertEqual(payload["allowed_reference_bed_only_ids"], ["RX"])
            self.assertTrue(payload["summary"]["reference_bed_only_allowlist_used"])
            self.assertEqual(payload["summary"]["reference_bed_rows"], 5)
            self.assertEqual(payload["summary"]["reference_nonzero_depth_gene_indices"], 2)


if __name__ == "__main__":
    unittest.main()
