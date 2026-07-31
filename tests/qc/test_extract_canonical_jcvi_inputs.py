from __future__ import annotations

import argparse
import json
import tempfile
import unittest
from pathlib import Path

from scripts.qc.extract_canonical_jcvi_inputs import CanonicalInputError, run


class CanonicalJcviInputTests(unittest.TestCase):
    def test_exact_t1_gene_closure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            protein = root / "publisher.faa"
            protein.write_text(
                ">g1.t1\nMAAA\n>g1.t2\nMBBB\n>g2.t1\nMCCC\n", encoding="utf-8"
            )
            gff = root / "publisher.gff3"
            gff.write_text(
                "##gff-version 3\n"
                "Chr01\tx\tgene\t1\t30\t.\t+\t.\tID=g1\n"
                "Chr01\tx\tmRNA\t1\t20\t.\t+\t.\tID=g1.t1;Parent=g1\n"
                "Chr01\tx\tmRNA\t2\t30\t.\t+\t.\tID=g1.t2;Parent=g1\n"
                "Chr02\tx\tgene\t5\t40\t.\t-\t.\tID=g2\n"
                "Chr02\tx\tmRNA\t5\t40\t.\t-\t.\tID=g2.t1;Parent=g2\n",
                encoding="utf-8",
            )
            output = root / "out"
            result = run(
                argparse.Namespace(
                    sample_id="test",
                    protein=protein,
                    gff=gff,
                    canonical_id_regex=r"\.t1$",
                    output_dir=output,
                )
            )
            self.assertEqual(result, output.resolve())
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["counts"]["selected_proteins"], 2)
            self.assertEqual(
                (output / "test.canonical.bed").read_text().splitlines(),
                ["Chr01\t0\t20\tg1.t1\t0\t+", "Chr02\t4\t40\tg2.t1\t0\t-"],
            )

    def test_incomplete_gene_closure_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            protein = root / "publisher.faa"
            protein.write_text(">g1.t1\nMAAA\n", encoding="utf-8")
            gff = root / "publisher.gff3"
            gff.write_text(
                "Chr01\tx\tgene\t1\t30\t.\t+\t.\tID=g1\n"
                "Chr01\tx\tmRNA\t1\t20\t.\t+\t.\tID=g1.t1;Parent=g1\n"
                "Chr02\tx\tgene\t1\t30\t.\t+\t.\tID=g2\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(CanonicalInputError, "complete GFF gene set"):
                run(
                    argparse.Namespace(
                        sample_id="test",
                        protein=protein,
                        gff=gff,
                        canonical_id_regex=r"\.t1$",
                        output_dir=root / "out",
                    )
                )


if __name__ == "__main__":
    unittest.main()
