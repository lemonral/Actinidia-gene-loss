"""Focused tests for the fail-closed chromosome-scope materializer."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[2] / "scripts" / "qc" / "materialize_chromosome_scope.py"
SPEC = importlib.util.spec_from_file_location("materialize_chromosome_scope", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class Fixture:
    def __init__(self, root: Path, *, gzip_inputs: bool = False) -> None:
        self.root = root
        self.genome = root / ("source.fa.gz" if gzip_inputs else "source.fa")
        self.gff = root / ("source.gff3.gz" if gzip_inputs else "source.gff3")
        self.mapping = root / "seqids.tsv"
        self.output = root / "materialized"
        genome_text = ">old_chr_1 source chromosome\nACGTACGTAA\n>unplaced_7\nNNNNNN\n"
        gff_text = (
            "##gff-version 3\n"
            "##sequence-region old_chr_1 1 10\n"
            "##sequence-region unplaced_7 1 6\n"
            "old_chr_1\ttest\tgene\t1\t10\t.\t+\t.\tID=gene1\n"
            "old_chr_1\ttest\tmRNA\t1\t10\t.\t+\t.\tID=tx1;Parent=gene1\n"
            "old_chr_1\ttest\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=tx1\n"
            "old_chr_1\ttest\tCDS\t7\t9\t.\t+\t0\tID=cds1;Parent=tx1\n"
            "unplaced_7\ttest\tgene\t1\t6\t.\t+\t.\tID=scaffold_gene\n"
        )
        if gzip_inputs:
            with gzip.open(self.genome, "wt", encoding="utf-8") as handle:
                handle.write(genome_text)
            with gzip.open(self.gff, "wt", encoding="utf-8") as handle:
                handle.write(gff_text)
        else:
            self.genome.write_text(genome_text, encoding="utf-8")
            self.gff.write_text(gff_text, encoding="utf-8")
        self.mapping.write_text(
            "source_seqid\tcanonical_seqid\nold_chr_1\tChr01\n", encoding="utf-8"
        )

    def run(self):
        return MODULE.materialize_chromosome_scope(
            genome_path=self.genome,
            gff_path=self.gff,
            map_path=self.mapping,
            output_dir=self.output,
            prefix="test_assembly",
        )


class ChromosomeScopeTests(unittest.TestCase):
    def test_filters_unmatched_scaffold_renames_and_audits(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary), gzip_inputs=True)
            genome_before = fixture.genome.read_bytes()
            gff_before = fixture.gff.read_bytes()
            result = fixture.run()
            self.assertEqual(result.retained_sequences, 1)
            self.assertEqual(result.excluded_sequences, 1)
            self.assertEqual(result.retained_features, 4)
            self.assertEqual(result.excluded_features, 1)
            with gzip.open(
                fixture.output / "test_assembly.genome.fa.gz", "rt", encoding="utf-8"
            ) as handle:
                genome = handle.read()
            self.assertIn(">Chr01 source chromosome", genome)
            self.assertNotIn("unplaced_7", genome)
            with gzip.open(
                fixture.output / "test_assembly.annotation.gff3.gz", "rt", encoding="utf-8"
            ) as handle:
                gff = handle.read()
            self.assertIn("##sequence-region Chr01 1 10", gff)
            self.assertIn("Chr01\ttest\tCDS", gff)
            self.assertNotIn("unplaced_7", gff)
            with (fixture.output / "audit" / "sequence_scope.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["scope"] for row in rows], ["retained", "excluded"])
            self.assertEqual(rows[1]["length_bp"], "6")
            self.assertEqual(rows[1]["gff_feature_count"], "1")
            validation = json.loads(
                (fixture.output / "audit" / "validation.json").read_text(encoding="utf-8")
            )
            self.assertEqual(validation["status"], "PASS")
            self.assertFalse(validation["policy"]["raw_inputs_modified"])
            self.assertEqual(fixture.genome.read_bytes(), genome_before)
            self.assertEqual(fixture.gff.read_bytes(), gff_before)
            self.assertTrue((fixture.output / "checksums.tsv").is_file())

    def test_multipart_cds_id_is_valid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            result = fixture.run()
            self.assertEqual(result.retained_features, 4)
            self.assertTrue(fixture.output.is_dir())

    def test_multi_parent_multipart_cds_is_valid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            text = fixture.gff.read_text(encoding="utf-8")
            text = text.replace(
                "old_chr_1\ttest\tmRNA\t1\t10\t.\t+\t.\tID=tx1;Parent=gene1\n",
                "old_chr_1\ttest\tmRNA\t1\t10\t.\t+\t.\tID=tx1;Parent=gene1\n"
                "old_chr_1\ttest\tmRNA\t1\t10\t.\t+\t.\tID=tx2;Parent=gene1\n",
            ).replace("Parent=tx1", "Parent=tx1,tx2")
            fixture.gff.write_text(text, encoding="utf-8")
            result = fixture.run()
            self.assertEqual(result.retained_features, 5)
            self.assertTrue(fixture.output.is_dir())

    def test_reproducible_gzip_payloads_are_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first_root = root / "first"
            second_root = root / "second"
            first_root.mkdir()
            second_root.mkdir()
            first = Fixture(first_root, gzip_inputs=True)
            second = Fixture(second_root, gzip_inputs=True)
            first.run()
            second.run()
            for relative_path in (
                "test_assembly.genome.fa.gz",
                "test_assembly.annotation.gff3.gz",
            ):
                self.assertEqual(
                    (first.output / relative_path).read_bytes(),
                    (second.output / relative_path).read_bytes(),
                )

    def test_explicit_different_genome_and_gff_seqids_are_paired(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.gff.write_text(
                fixture.gff.read_text(encoding="utf-8").replace(
                    "old_chr_1", "annotation_chr_1"
                ),
                encoding="utf-8",
            )
            fixture.mapping.write_text(
                "genome_seqid\tgff_seqid\tcanonical_seqid\n"
                "old_chr_1\tannotation_chr_1\tChr01\n",
                encoding="utf-8",
            )
            result = fixture.run()
            self.assertEqual(result.retained_sequences, 1)
            self.assertEqual(result.retained_features, 4)
            with gzip.open(
                fixture.output / "test_assembly.annotation.gff3.gz",
                "rt",
                encoding="utf-8",
            ) as handle:
                gff = handle.read()
            self.assertIn("##sequence-region Chr01 1 10", gff)
            self.assertIn("Chr01\ttest\tgene", gff)
            self.assertNotIn("annotation_chr_1", gff)
            validation = json.loads(
                (fixture.output / "audit" / "validation.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                validation["policy"]["seqid_map_schema"],
                "explicit_genome_gff_canonical_v1",
            )
            self.assertFalse(
                validation["policy"]["seqid_correspondence_inferred_by_order"]
            )

    def test_mapped_gff_seqid_without_features_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.mapping.write_text(
                "genome_seqid\tgff_seqid\tcanonical_seqid\n"
                "old_chr_1\tmissing_annotation_chr\tChr01\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                MODULE.ChromosomeScopeError, "without any GFF3 feature row"
            ):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_explicit_blank_gff_seqid_retains_featureless_genome_record(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.gff.write_text(
                fixture.gff.read_text(encoding="utf-8").replace(
                    "unplaced_7\ttest\tgene\t1\t6\t.\t+\t.\tID=scaffold_gene\n",
                    "",
                ),
                encoding="utf-8",
            )
            fixture.mapping.write_text(
                "genome_seqid\tgff_seqid\tcanonical_seqid\n"
                "old_chr_1\told_chr_1\tChr01\n"
                "unplaced_7\t\tFeaturelessContig07\n",
                encoding="utf-8",
            )

            result = fixture.run()

            self.assertEqual(result.retained_sequences, 2)
            self.assertEqual(result.excluded_sequences, 0)
            self.assertEqual(result.retained_features, 4)
            self.assertEqual(result.excluded_features, 0)
            with gzip.open(
                fixture.output / "test_assembly.genome.fa.gz", "rt", encoding="utf-8"
            ) as handle:
                genome = handle.read()
            self.assertIn(">Chr01 source chromosome", genome)
            self.assertIn(">FeaturelessContig07", genome)
            with gzip.open(
                fixture.output / "test_assembly.annotation.gff3.gz",
                "rt",
                encoding="utf-8",
            ) as handle:
                gff = handle.read()
            self.assertNotIn("FeaturelessContig07", gff)
            self.assertNotIn("unplaced_7", gff)

            with (fixture.output / "audit" / "sequence_scope.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            by_genome = {row["genome_seqid"]: row for row in rows}
            featureless = by_genome["unplaced_7"]
            self.assertEqual(featureless["scope"], "retained")
            self.assertEqual(featureless["gff_seqid"], "")
            self.assertEqual(featureless["canonical_seqid"], "FeaturelessContig07")
            self.assertEqual(featureless["gff_feature_count"], "0")

            validation = json.loads(
                (fixture.output / "audit" / "validation.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(validation["counts"]["map_rows"], 2)
            self.assertEqual(validation["counts"]["nonempty_mapped_gff_seqids"], 1)
            self.assertEqual(
                validation["counts"][
                    "retained_genome_sequences_declared_without_gff"
                ],
                1,
            )
            self.assertEqual(
                validation["checks"][
                    "retained_genome_seqids_declared_without_gff_have_zero_feature_rows"
                ],
                "PASS",
            )

    def test_feature_on_genome_seqid_declared_without_gff_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.mapping.write_text(
                "genome_seqid\tgff_seqid\tcanonical_seqid\n"
                "old_chr_1\told_chr_1\tChr01\n"
                "unplaced_7\t\tFeaturelessContig07\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                MODULE.ChromosomeScopeError, "declared without GFF features"
            ):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_featureless_genome_and_mapped_gff_cross_domain_collisions_are_rejected(
        self,
    ) -> None:
        map_bodies = (
            "unplaced_7\t\tFeaturelessContig07\n"
            "old_chr_1\tunplaced_7\tChr01\n",
            "old_chr_1\tunplaced_7\tChr01\n"
            "unplaced_7\t\tFeaturelessContig07\n",
        )
        for map_body in map_bodies:
            with self.subTest(map_body=map_body):
                with tempfile.TemporaryDirectory() as temporary:
                    fixture = Fixture(Path(temporary))
                    fixture.mapping.write_text(
                        "genome_seqid\tgff_seqid\tcanonical_seqid\n" + map_body,
                        encoding="utf-8",
                    )
                    with self.assertRaisesRegex(
                        MODULE.ChromosomeScopeError,
                        "declared without GFF features|without GFF features because",
                    ):
                        fixture.run()
                    self.assertFalse(fixture.output.exists())

    def test_duplicate_nonempty_explicit_gff_seqid_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.mapping.write_text(
                "genome_seqid\tgff_seqid\tcanonical_seqid\n"
                "old_chr_1\told_chr_1\tChr01\n"
                "unplaced_7\told_chr_1\tChr02\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "duplicate gff_seqid"):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_duplicate_canonical_after_featureless_row_is_rejected_cleanly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.mapping.write_text(
                "genome_seqid\tgff_seqid\tcanonical_seqid\n"
                "unplaced_7\t\tChr01\n"
                "old_chr_1\told_chr_1\tChr01\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "already mapped"):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_whitespace_only_explicit_gff_seqid_is_not_treated_as_blank(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.mapping.write_text(
                "genome_seqid\tgff_seqid\tcanonical_seqid\n"
                "old_chr_1\t   \tChr01\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "whitespace-only"):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_legacy_map_does_not_gain_featureless_semantics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.gff.write_text(
                "##gff-version 3\n"
                "unplaced_7\ttest\tgene\t1\t6\t.\t+\t.\tID=scaffold_gene\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                MODULE.ChromosomeScopeError, "without any GFF3 feature row"
            ):
                fixture.run()
            self.assertFalse(fixture.output.exists())

        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.mapping.write_text(
                "source_seqid\tcanonical_seqid\n\tChr01\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(
                MODULE.ChromosomeScopeError, "genome_seqid must be one non-empty"
            ):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_out_of_bounds_feature_publishes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.gff.write_text(
                fixture.gff.read_text(encoding="utf-8").replace(
                    "old_chr_1\ttest\tgene\t1\t10", "old_chr_1\ttest\tgene\t1\t11"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "outside"):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_map_source_missing_from_genome_publishes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            with fixture.mapping.open("a", encoding="utf-8") as handle:
                handle.write("absent_chr\tChr02\n")
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "absent from the genome"):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_broken_parent_publishes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.gff.write_text(
                fixture.gff.read_text(encoding="utf-8").replace(
                    "ID=tx1;Parent=gene1", "ID=tx1;Parent=missing_gene"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "has no matching ID"):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_duplicate_map_target_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            with fixture.mapping.open("a", encoding="utf-8") as handle:
                handle.write("unplaced_7\tChr01\n")
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "already mapped"):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_duplicate_genome_seqid_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            with fixture.genome.open("a", encoding="utf-8") as handle:
                handle.write(">old_chr_1 duplicate\nACGT\n")
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "duplicate FASTA seqid"):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_unknown_gff_seqid_and_embedded_fasta_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            with fixture.gff.open("a", encoding="utf-8") as handle:
                handle.write("unknown\ttest\tgene\t1\t2\t.\t+\t.\tID=unknown_gene\n")
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "no genome sequence"):
                fixture.run()
            self.assertFalse(fixture.output.exists())
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            with fixture.gff.open("a", encoding="utf-8") as handle:
                handle.write("##FASTA\n>old_chr_1\nACGT\n")
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "embedded GFF3 FASTA"):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_existing_output_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            fixture.output.mkdir()
            marker = fixture.output / "keep.txt"
            marker.write_text("unchanged\n", encoding="utf-8")
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "refusing overwrite"):
                fixture.run()
            self.assertEqual(marker.read_text(encoding="utf-8"), "unchanged\n")

    def test_dangling_output_symlink_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary))
            missing_target = fixture.root / "missing-target"
            fixture.output.symlink_to(missing_target, target_is_directory=True)
            with self.assertRaisesRegex(MODULE.ChromosomeScopeError, "refusing overwrite"):
                fixture.run()
            self.assertTrue(fixture.output.is_symlink())
            self.assertEqual(fixture.output.readlink(), missing_target)


if __name__ == "__main__":
    unittest.main()
