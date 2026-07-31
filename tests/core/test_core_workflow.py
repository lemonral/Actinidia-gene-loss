"""Standard-library smoke tests for the refactored core gene-loss path."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from geneloss_repro.annotation import extract_annotation
from geneloss_repro.blast import classify_tblastx, iter_blast_hits, summarize_classification
from geneloss_repro.io_utils import concatenate_tsv, read_tsv
from geneloss_repro.master import build_loss_master
from geneloss_repro.spatial import spatial_summary
from geneloss_repro.synorth import call_candidates, normalize_synorth


REFERENCE_GFF = """##gff-version 3
ChrR\tunit\tgene\t1\t9\t.\t+\t.\tID=R1
ChrR\tunit\tmRNA\t1\t9\t.\t+\t.\tID=R1.t1;Parent=R1
ChrR\tunit\tCDS\t1\t9\t.\t+\t0\tID=R1.cds;Parent=R1.t1
ChrR\tunit\tgene\t20\t28\t.\t+\t.\tID=R2
ChrR\tunit\tmRNA\t20\t28\t.\t+\t.\tID=R2.t1;Parent=R2
ChrR\tunit\tCDS\t20\t28\t.\t+\t0\tID=R2.cds;Parent=R2.t1
ChrR\tunit\tgene\t40\t48\t.\t+\t.\tID=R3
ChrR\tunit\tmRNA\t40\t48\t.\t+\t.\tID=R3.t1;Parent=R3
ChrR\tunit\tCDS\t40\t48\t.\t+\t0\tID=R3.cds;Parent=R3.t1
ChrR\tunit\tgene\t60\t68\t.\t+\t.\tID=R4
ChrR\tunit\tmRNA\t60\t68\t.\t+\t.\tID=R4.t1;Parent=R4
ChrR\tunit\tCDS\t60\t68\t.\t+\t0\tID=R4.cds;Parent=R4.t1
"""

TARGET_GFF = """##gff-version 3
ChrT\tunit\tgene\t1\t20\t.\t+\t.\tID=T1
ChrT\tunit\tgene\t25\t45\t.\t+\t.\tID=T2
ChrT\tunit\tgene\t50\t70\t.\t+\t.\tID=T3
ChrT\tunit\tgene\t75\t95\t.\t+\t.\tID=T4
"""


class CoreWorkflowTest(unittest.TestCase):
    def test_end_to_end_12_column_and_spatial(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            genome = root / "reference.fa"
            genome.write_text(">ChrR\n" + "ATG" * 30 + "\n", encoding="utf-8")
            gff = root / "reference.gff3"
            gff.write_text(REFERENCE_GFF, encoding="utf-8")
            annotation = extract_annotation(genome, gff, root / "annotation", "Reference")
            self.assertEqual(len(read_tsv(annotation["coords"])), 4)

            raw_synorth = root / "raw.synorth.txt"
            raw_synorth.write_text(
                "R1.t1\tChrR\t1\t9\tT1\tChrT\t1\t20\t1e-40\t1|1\t1|1\t+\tBest_Hit\n"
                "R3.t1\tChrR\t40\t48\tT3\tChrT\t50\t70\t1e-30\t1|1\t1|1\t+\tBest_Hit\n",
                encoding="utf-8",
            )
            anchors = root / "anchors.tsv"
            _, side = normalize_synorth(raw_synorth, annotation["coords"], "Target_A", anchors)
            self.assertEqual(side, "first")
            candidates = root / "candidates.tsv"
            called, retained = call_candidates(annotation["coords"], anchors, candidates, flank_genes=1, mode="bracketed")
            self.assertEqual([row["reference_gene"] for row in called], ["R2.t1"])
            self.assertEqual(len(retained), 2)

            blast = root / "hits.12.txt"
            blast.write_text("R2.t1\tChrT\t70\t90\t1\t0\t1\t90\t65\t25\t1e-20\t100\n", encoding="utf-8")
            classification = root / "classification.tsv"
            schema = root / "schema.tsv"
            classified = classify_tblastx(candidates, blast, classification, schema, blast_schema="blast12")
            self.assertEqual(classified[0]["classification"], "pseudogenized")
            self.assertEqual(classified[0]["best_subject_start"], 25)
            self.assertEqual(classified[0]["best_subject_end"], 65)

            retained_path = candidates.with_name(candidates.stem + ".retained_anchors.tsv")
            master = root / "gene_loss_master.tsv"
            master_rows = build_loss_master(
                annotation["coords"], classification, retained_path, master,
                noncandidate_class="retained_by_synorth", run_id="test_master",
            )
            self.assertEqual(len(master_rows), 4)
            self.assertEqual(sum(row["rate_eligible"] == "true" for row in master_rows), 4)
            master_by_gene = {row["reference_gene_id"]: row for row in master_rows}
            self.assertEqual(master_by_gene["R2.t1"]["classification"], "pseudogenized")
            self.assertEqual(master_by_gene["R4.t1"]["classification"], "retained_by_synorth")
            self.assertEqual(master_by_gene["R4.t1"]["run_id"], "test_master")

            summary = root / "summary.tsv"
            summarize_classification(classification, annotation["coords"], summary)
            summary_row = read_tsv(summary)[0]
            self.assertEqual(summary_row["pseudogenized_count"], "1")
            self.assertEqual(summary_row["reference_gene_count"], "4")

            target_gff = root / "target.gff3"
            target_gff.write_text(TARGET_GFF, encoding="utf-8")
            lengths = root / "lengths.tsv"
            lengths.write_text("chromosome\tlength\nChrT\t100\n", encoding="utf-8")
            outputs = spatial_summary(classification, target_gff, root / "spatial", chromosome_lengths_path=lengths)
            inter = read_tsv(outputs["inter"])
            self.assertEqual(sum(int(row["pseudogenized_reference_genes"]) for row in inter), 1)
            intra = read_tsv(outputs["intra"])
            self.assertEqual(len(intra), 5)
            self.assertEqual({row["bin_mode"] for row in intra}, {"equal-width"})

    def test_legacy_six_column_schema_is_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "hits.6.txt"
            path.write_text("Q1\tChr1\t55\t120\t88\t1e-20\n", encoding="utf-8")
            hits, detected = iter_blast_hits(path, schema="legacy6-bitscore-evalue")
            parsed = list(hits)
            self.assertEqual(detected, "legacy6-bitscore-evalue")
            self.assertEqual(parsed[0].bitscore, 88.0)
            self.assertEqual(parsed[0].evalue, 1e-20)

    def test_candidate_absent_from_actual_query_fasta_is_uncertain(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            candidates = root / "candidates.tsv"
            candidates.write_text("target_sample\treference_gene\nTarget_A\tQ_missing\n", encoding="utf-8")
            query_fasta = root / "queries.fa"
            query_fasta.write_text(">Different_query\nATG\n", encoding="utf-8")
            blast = root / "hits.txt"
            blast.write_text("Q_missing\tChr1\t99\t100\t0\t0\t1\t100\t1\t100\t1e-30\t300\n", encoding="utf-8")
            output = root / "classification.tsv"
            rows = classify_tblastx(candidates, blast, output, root / "schema.tsv", blast_schema="blast12", query_fasta_path=query_fasta)
            self.assertEqual(rows[0]["classification"], "uncertain")
            self.assertEqual(rows[0]["decision_reason"], "candidate_absent_from_tblastx_query_fasta")

    def test_concat_rejects_mixed_schema(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            first = root / "first.tsv"
            second = root / "second.tsv"
            first.write_text("sample_id\tvalue\nA\t1\n", encoding="utf-8")
            second.write_text("sample_id\tvalue\nB\t2\n", encoding="utf-8")
            output = root / "combined.tsv"
            concatenate_tsv([first, second], output)
            self.assertEqual(len(read_tsv(output)), 2)
            second.write_text("sample_id\tdifferent_value\nB\t2\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                concatenate_tsv([first, second], output)


if __name__ == "__main__":
    unittest.main()
