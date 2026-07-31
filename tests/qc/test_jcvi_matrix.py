from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from geneloss_repro.jcvi_matrix import (
    JcviMatrixError,
    build_jcvi_rows,
    read_bed,
    read_normalized_anchor_pairs,
    relabel_reference_bed_from_canonical_truth,
    require_bed_gff_identity,
    require_bed_protein_identity,
)


class JcviMatrixTests(unittest.TestCase):
    def _inputs(self, root: Path):
        target_bed = root / "target.bed"
        reference_bed = root / "reference.bed"
        target_protein = root / "target.faa"
        target_gff = root / "target.gff3"
        target_bed.write_text(
            "Pub01\t0\t10\tt1\t0\t+\nPub01\t20\t30\tt2\t0\t-\n"
            "Pub02\t0\t10\tt3\t0\t+\n",
            encoding="utf-8",
        )
        reference_bed.write_text(
            "Ref01\t0\t10\tr1\t0\t+\nRef01\t20\t30\tr2\t0\t-\n"
            "Ref02\t0\t10\tr3\t0\t+\n",
            encoding="utf-8",
        )
        target_protein.write_text(">t1\nM\n>t2\nM\n>t3\nM\n", encoding="utf-8")
        target_gff.write_text(
            "##gff-version 3\n"
            "Pub01\tx\tmRNA\t1\t10\t.\t+\t.\tID=t1;Parent=g1\n"
            "Pub01\tx\tmRNA\t21\t30\t.\t-\t.\tID=t2;Parent=g2\n"
            "Pub02\tx\tmRNA\t1\t10\t.\t+\t.\tID=t3;Parent=g3\n",
            encoding="utf-8",
        )
        return target_bed, reference_bed, target_protein, target_gff

    def test_bidirectional_pairs_are_role_normalized_and_unioned(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target_path, reference_path, protein, gff = self._inputs(root)
            target = read_bed(target_path)
            reference = read_bed(reference_path)
            require_bed_protein_identity(target, protein, label="target")
            require_bed_gff_identity(target, gff)
            forward = root / "forward.anchors"
            reverse = root / "reverse.anchors"
            forward.write_text("###\nr1\tt1\nr2\tt2\n", encoding="utf-8")
            reverse.write_text("###\nt1\tr1\nt3\tr3\n", encoding="utf-8")
            forward_pairs, _ = read_normalized_anchor_pairs(
                forward, first_bed=reference, second_bed=target, first_role="reference"
            )
            reverse_pairs, _ = read_normalized_anchor_pairs(
                reverse, first_bed=target, second_bed=reference, first_role="target"
            )
            rows, audit = build_jcvi_rows(
                target_bed=target,
                reference_bed=reference,
                canonical_by_reference={"Ref01": "Chr01", "Ref02": "Chr02"},
                forward_pairs=forward_pairs,
                reverse_pairs=reverse_pairs,
            )
            self.assertEqual(len(rows), 4)
            self.assertEqual(audit["bidirectional_union_pairs"], 3)
            self.assertEqual(audit["bidirectional_intersection_pairs"], 1)
            first = rows[0]
            self.assertEqual(first["query_anchored_genes"], "2")
            self.assertEqual(first["reference_anchored_genes"], "2")
            self.assertEqual(first["unique_anchor_pairs"], "2")
            self.assertEqual(first["score"], "1")

    def test_duplicate_raw_pair_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target_path, reference_path, _, _ = self._inputs(root)
            target = read_bed(target_path)
            reference = read_bed(reference_path)
            anchors = root / "duplicate.anchors"
            anchors.write_text("r1\tt1\nr1\tt1\n", encoding="utf-8")
            with self.assertRaisesRegex(JcviMatrixError, "duplicate raw anchor pair"):
                read_normalized_anchor_pairs(
                    anchors,
                    first_bed=reference,
                    second_bed=target,
                    first_role="reference",
                )

    def test_exact_canonical_labels_map_to_haplome_reference_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, reference_path, _, _ = self._inputs(root)
            reference = read_bed(reference_path)
            relabelled, aliases = relabel_reference_bed_from_canonical_truth(
                reference, {"Chr01A": "Ref01", "Chr02A": "Ref02"}
            )
            self.assertEqual(relabelled.chromosomes, ("Chr01A", "Chr02A"))
            self.assertEqual(relabelled.gene_to_chromosome["r1"], "Chr01A")
            self.assertEqual(aliases, {"Ref01": "Chr01A", "Ref02": "Chr02A"})

    def test_nonexact_reference_alias_set_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, reference_path, _, _ = self._inputs(root)
            reference = read_bed(reference_path)
            with self.assertRaisesRegex(JcviMatrixError, "exact canonical-label set"):
                relabel_reference_bed_from_canonical_truth(
                    reference, {"Chr01A": "Ref01", "Chr02A": "Other"}
                )


if __name__ == "__main__":
    unittest.main()
