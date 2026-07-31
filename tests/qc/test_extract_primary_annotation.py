"""Tests for deterministic, atomic primary annotation standardization."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from geneloss_repro.primary_annotation import (
    PrimaryAnnotationError,
    parse_canonical_rule,
    standardize_primary_annotation,
)


def replace_segment(sequence: list[str], start: int, content: str) -> None:
    sequence[start - 1:start - 1 + len(content)] = list(content)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier = ""
    parts: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            if identifier:
                records[identifier] = "".join(parts)
            identifier = raw[1:].split()[0]
            parts = []
        else:
            parts.append(raw.strip())
    if identifier:
        records[identifier] = "".join(parts)
    return records


class MultiIsoformFixture:
    def __init__(self, root: Path, *, gzip_inputs: bool = True) -> None:
        self.root = root
        self.output = root / "primary"
        self.genome = root / ("scope.fa.gz" if gzip_inputs else "scope.fa")
        self.gff = root / ("scope.gff3.gz" if gzip_inputs else "scope.gff3")
        sequence = list("N" * 180)
        replace_segment(sequence, 1, "ATGAAACCCGGGTAA")
        replace_segment(sequence, 21, "ATGAAATAA")
        replace_segment(sequence, 50, "ATGAAATAACCCGGGTAA")
        replace_segment(sequence, 80, "ATGCCCAAATAA")
        replace_segment(sequence, 110, "TTAGGGCAT")
        replace_segment(sequence, 130, "ATGG")
        replace_segment(sequence, 140, "CCAAA")
        fasta = ">Chr01 chromosome scope\n" + "".join(sequence) + "\n"
        gff = (
            "##gff-version 3\n"
            "Chr01\ttest\tgene\t1\t40\t.\t+\t.\tID=gene1\n"
            "Chr01\ttest\tmRNA\t1\t15\t.\t+\t.\tID=tx_long;Parent=gene1\n"
            "Chr01\ttest\tCDS\t1\t15\t.\t+\t0\tParent=tx_long\n"
            "Chr01\ttest\tmRNA\t21\t29\t.\t+\t.\tID=tx_canonical;Parent=gene1;tag=canonical\n"
            "Chr01\ttest\tCDS\t21\t29\t.\t+\t0\tParent=tx_canonical\n"
            "Chr01\ttest\tgene\t50\t100\t.\t+\t.\tID=gene2\n"
            "Chr01\ttest\tmRNA\t50\t67\t.\t+\t.\tID=tx_invalid_long;Parent=gene2\n"
            "Chr01\ttest\tCDS\t50\t67\t.\t+\t0\tParent=tx_invalid_long\n"
            "Chr01\ttest\tmRNA\t80\t91\t.\t+\t.\tID=tx_valid_short;Parent=gene2\n"
            "Chr01\ttest\tCDS\t80\t91\t.\t+\t0\tParent=tx_valid_short\n"
            "Chr01\ttest\tgene\t110\t118\t.\t-\t.\tID=gene3\n"
            "Chr01\ttest\tmRNA\t110\t118\t.\t-\t.\tID=tx_minus;Parent=gene3\n"
            "Chr01\ttest\tCDS\t110\t118\t.\t-\t0\tParent=tx_minus\n"
            "Chr01\ttest\tgene\t130\t144\t.\t+\t.\tID=gene4\n"
            "Chr01\ttest\tmRNA\t130\t144\t.\t+\t.\tID=tx_phased;Parent=gene4\n"
            "Chr01\ttest\tCDS\t130\t133\t.\t+\t0\tParent=tx_phased\n"
            "Chr01\ttest\tCDS\t140\t144\t.\t+\t2\tParent=tx_phased\n"
            "Chr01\ttest\tgene\t150\t170\t.\t+\t.\tID=gene5\n"
            "Chr01\ttest\tmRNA\t150\t170\t.\t+\t.\tID=tx_noncoding;Parent=gene5\n"
        )
        if gzip_inputs:
            with gzip.open(self.genome, "wt", encoding="utf-8") as handle:
                handle.write(fasta)
            with gzip.open(self.gff, "wt", encoding="utf-8") as handle:
                handle.write(gff)
        else:
            self.genome.write_text(fasta, encoding="utf-8")
            self.gff.write_text(gff, encoding="utf-8")

    def run(self, **kwargs):
        return standardize_primary_annotation(
            self.genome,
            self.gff,
            self.output,
            "act_test_hap1",
            canonical_rules=(parse_canonical_rule("tag=canonical"),),
            gffread="none",
            **kwargs,
        )


def write_single_fixture(root: Path, *, phase: str = "0") -> tuple[Path, Path]:
    genome = root / "single.fa"
    gff = root / "single.gff3"
    genome.write_text(">Chr01\nATGAAATAA\n", encoding="utf-8")
    gff.write_text(
        "##gff-version 3\n"
        "Chr01\ttest\tgene\t1\t9\t.\t+\t.\tID=gene1\n"
        "Chr01\ttest\tmRNA\t1\t9\t.\t+\t.\tID=tx1;Parent=gene1\n"
        f"Chr01\ttest\tCDS\t1\t9\t.\t+\t{phase}\tParent=tx1\n",
        encoding="utf-8",
    )
    return genome, gff


def write_fake_gffread(
    root: Path,
    *,
    matching: bool,
    reject_feature: str | None = None,
    extra_output_id: str | None = None,
    output_id: str = "tx1",
) -> Path:
    executable = root / ("fake_gffread_match.py" if matching else "fake_gffread_mismatch.py")
    cds = "ATGAAATAA" if matching else "ATGCCCTAA"
    protein = "MK*" if matching else "MP*"
    cds_text = f">{output_id}\n{cds}\n"
    protein_text = f">{output_id}\n{protein}\n"
    if extra_output_id is not None:
        cds_text += f">{extra_output_id}\nATGAAATAA\n"
        protein_text += f">{extra_output_id}\nMK*\n"
    feature_check = ""
    if reject_feature is not None:
        feature_check = (
            "gff_text = pathlib.Path(sys.argv[1]).read_text(encoding='utf-8')\n"
            f"if '\\t{reject_feature}\\t' in gff_text:\n"
            "    raise SystemExit(91)\n"
        )
    executable.write_text(
        f"#!{sys.executable}\n"
        "import pathlib, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('gffread fake-1.0')\n"
        "    raise SystemExit(0)\n"
        f"{feature_check}"
        "cds_path = pathlib.Path(sys.argv[sys.argv.index('-x') + 1])\n"
        "protein_path = pathlib.Path(sys.argv[sys.argv.index('-y') + 1])\n"
        f"cds_path.write_text({cds_text!r}, encoding='utf-8')\n"
        f"protein_path.write_text({protein_text!r}, encoding='utf-8')\n",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    return executable


class PrimaryAnnotationTests(unittest.TestCase):
    def test_cli_help_works_outside_repository_without_pythonpath(self) -> None:
        script = Path(__file__).parents[2] / "scripts" / "qc" / "extract_primary_annotation.py"
        environment = dict(os.environ)
        environment.pop("PYTHONPATH", None)
        with tempfile.TemporaryDirectory() as temporary:
            completed = subprocess.run(
                [sys.executable, str(script), "--help"],
                cwd=temporary,
                env=environment,
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--require-gffread", completed.stdout)
        self.assertIn("--gene-as-transcript", completed.stdout)

    def test_validates_before_selection_and_prefers_configured_canonical_tag(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = MultiIsoformFixture(Path(temporary), gzip_inputs=True)
            result = fixture.run()
            self.assertEqual(result.source_gene_count, 5)
            self.assertEqual(result.selected_gene_count, 4)
            self.assertEqual(result.invalid_coding_gene_count, 0)
            proteins = read_fasta(fixture.output / "act_test_hap1.protein.faa")
            self.assertEqual(
                proteins,
                {
                    "tx_canonical": "MK",
                    "tx_valid_short": "MPK",
                    "tx_minus": "MP",
                    "tx_phased": "MAK",
                },
            )
            with (fixture.output / "act_test_hap1.transcript_audit.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                audit = {row["transcript_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(audit["tx_long"]["disposition"], "valid_not_selected")
            self.assertEqual(audit["tx_canonical"]["disposition"], "selected")
            self.assertEqual(audit["tx_invalid_long"]["validation_status"], "internal_stop_codon")
            self.assertEqual(audit["tx_valid_short"]["disposition"], "selected")
            self.assertEqual(audit["tx_noncoding"]["disposition"], "noncoding_or_no_CDS")
            primary_gff = (fixture.output / "act_test_hap1.primary.gff3").read_text(encoding="utf-8")
            self.assertIn("Chr01\tprimary_isoform_standardizer\tgene\t1\t40", primary_gff)
            manifest_text = (fixture.output / "run_manifest.json").read_text(encoding="utf-8")
            manifest = json.loads(manifest_text)
            self.assertEqual(manifest["gffread_comparison"]["status"], "NOT_RUN_EXPLICITLY_DISABLED")
            self.assertEqual(manifest["publication_gate"], "BLOCKED_GFFREAD_NOT_RUN")
            self.assertNotIn(str(Path(temporary)), manifest_text)
            with (fixture.output / "checksums.tsv").open(encoding="utf-8", newline="") as handle:
                checksums = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(checksums)
            self.assertTrue(all((fixture.output / row["file"]).is_file() for row in checksums))

    def test_all_invalid_coding_gene_fails_without_publishing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome, gff = write_single_fixture(root)
            genome.write_text(">Chr01\nATGTAACCC\n", encoding="utf-8")
            output = root / "output"
            with self.assertRaisesRegex(PrimaryAnnotationError, "no valid coding transcript"):
                standardize_primary_annotation(genome, gff, output, "sample", gffread="none")
            self.assertFalse(output.exists())
            self.assertFalse(list(root.glob(".output.staging.*")))

    def test_final_tie_break_uses_transcript_id_not_gff_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome = root / "tie.fa"
            gff = root / "tie.gff3"
            genome.write_text(">Chr01\nATGAAATAANATGAAATAA\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "Chr01\ttest\tgene\t1\t19\t.\t+\t.\tID=gene1\n"
                "Chr01\ttest\tmRNA\t1\t9\t.\t+\t.\tID=tx_z;Parent=gene1\n"
                "Chr01\ttest\tCDS\t1\t9\t.\t+\t0\tParent=tx_z\n"
                "Chr01\ttest\tmRNA\t11\t19\t.\t+\t.\tID=tx_a;Parent=gene1\n"
                "Chr01\ttest\tCDS\t11\t19\t.\t+\t0\tParent=tx_a\n",
                encoding="utf-8",
            )
            output = root / "output"
            standardize_primary_annotation(genome, gff, output, "sample", gffread="none")
            self.assertEqual(list(read_fasta(output / "sample.protein.faa")), ["tx_a"])

    def test_explicit_omit_policy_retains_audited_invalid_gene(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = MultiIsoformFixture(root, gzip_inputs=False)
            with fixture.gff.open("a", encoding="utf-8") as handle:
                handle.write(
                    "Chr01\ttest\tgene\t171\t179\t.\t+\t.\tID=gene6\n"
                    "Chr01\ttest\tmRNA\t171\t179\t.\t+\t.\tID=tx_bad_only;Parent=gene6\n"
                    "Chr01\ttest\tCDS\t171\t179\t.\t+\t0\tParent=tx_bad_only\n"
                )
            # NNN translates to X and is valid-but-flagged; make the only model
            # frame-invalid instead so omission is unambiguous.
            text = fixture.gff.read_text(encoding="utf-8").replace(
                "CDS\t171\t179", "CDS\t171\t178"
            )
            fixture.gff.write_text(text, encoding="utf-8")
            result = fixture.run(invalid_coding_gene_policy="omit")
            self.assertEqual(result.invalid_coding_gene_count, 1)
            with (fixture.output / "act_test_hap1.gene_audit.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["gene6"]["status"], "no_valid_coding_transcript")

    def test_missing_phase_requires_explicit_zero_compatibility_policy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome, gff = write_single_fixture(root, phase=".")
            strict_output = root / "strict"
            with self.assertRaisesRegex(PrimaryAnnotationError, "no valid coding transcript"):
                standardize_primary_annotation(genome, gff, strict_output, "sample", gffread="none")
            self.assertFalse(strict_output.exists())
            compatible_output = root / "compatible"
            standardize_primary_annotation(
                genome,
                gff,
                compatible_output,
                "sample",
                missing_phase_policy="zero",
                gffread="none",
            )
            with (compatible_output / "sample.primary_isoforms.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertIn("missing_phase_treated_as_zero", row["QC_flags"])

    def test_matching_gffread_is_required_and_recorded_without_its_path(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome, gff = write_single_fixture(root)
            executable = write_fake_gffread(root, matching=True)
            output = root / "output"
            result = standardize_primary_annotation(
                genome,
                gff,
                output,
                "sample",
                gffread=str(executable),
                require_gffread=True,
            )
            self.assertEqual(result.gffread_status, "PASS")
            with (output / "sample.gffread_comparison.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertTrue(all(row["status"] == "PASS" for row in rows))
            manifest_text = (output / "run_manifest.json").read_text(encoding="utf-8")
            self.assertNotIn(str(executable), manifest_text)
            manifest = json.loads(manifest_text)
            self.assertEqual(manifest["gffread_comparison"]["version"], "gffread fake-1.0")
            self.assertEqual(
                manifest["gffread_comparison"]["comparison_annotation_scope"],
                "selected_primary_GFF3_gene_mRNA_CDS_only",
            )
            self.assertEqual(
                manifest["gffread_comparison"]["comparison_GFF3_file_name"],
                "sample.primary.gff3",
            )
            self.assertFalse(
                manifest["gffread_comparison"]["source_full_GFF3_passed_to_gffread"]
            )
            self.assertTrue(
                manifest["gffread_comparison"]["exact_selected_CDS_ID_set"]
            )
            self.assertTrue(
                manifest["gffread_comparison"]["exact_selected_protein_ID_set"]
            )
            self.assertEqual(manifest["publication_gate"], "PASS")

    def test_unrelated_malformed_start_codon_is_not_passed_to_gffread(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome, gff = write_single_fixture(root)
            genome.write_text(">Chr01\nATGAAATAA" + "N" * 21 + "\n", encoding="utf-8")
            with gff.open("a", encoding="utf-8") as handle:
                # This non-CDS row lies outside the declared mRNA and has an
                # invalid 21-nt codon span. It is retained in the publisher
                # input audit but must not enter selected-CDS extraction.
                handle.write(
                    "Chr01\ttest\tstart_codon\t10\t30\t.\t+\t.\tParent=tx1\n"
                )
            executable = write_fake_gffread(
                root, matching=True, reject_feature="start_codon"
            )
            output = root / "output"
            result = standardize_primary_annotation(
                genome,
                gff,
                output,
                "sample",
                gffread=str(executable),
                require_gffread=True,
            )
            self.assertEqual(result.gffread_status, "PASS")
            primary_gff_path = output / "sample.primary.gff3"
            primary_gff = primary_gff_path.read_text(encoding="utf-8")
            self.assertNotIn("\tstart_codon\t", primary_gff)
            self.assertEqual(
                {
                    line.split("\t")[2]
                    for line in primary_gff.splitlines()
                    if not line.startswith("#")
                },
                {"gene", "mRNA", "CDS"},
            )
            manifest_text = (output / "run_manifest.json").read_text(encoding="utf-8")
            manifest = json.loads(manifest_text)
            comparison = manifest["gffread_comparison"]
            self.assertEqual(manifest["workflow_version"], "1.2.0")
            self.assertEqual(
                comparison["comparison_annotation_scope"],
                "selected_primary_GFF3_gene_mRNA_CDS_only",
            )
            self.assertEqual(comparison["comparison_GFF3_file_name"], primary_gff_path.name)
            self.assertEqual(
                comparison["comparison_GFF3_sha256"],
                hashlib.sha256(primary_gff_path.read_bytes()).hexdigest(),
            )
            self.assertFalse(comparison["source_full_GFF3_passed_to_gffread"])
            self.assertNotIn(str(root), manifest_text)

    def test_gffread_extra_output_id_rejects_atomic_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome, gff = write_single_fixture(root)
            executable = write_fake_gffread(
                root, matching=True, extra_output_id="unexpected_tx"
            )
            output = root / "output"
            with self.assertRaisesRegex(
                PrimaryAnnotationError,
                "ID set does not exactly equal the selected-primary set.*extra=1",
            ):
                standardize_primary_annotation(
                    genome,
                    gff,
                    output,
                    "sample",
                    gffread=str(executable),
                    require_gffread=True,
                )
            self.assertFalse(output.exists())
            self.assertFalse(list(root.glob(".output.staging.*")))

    def test_gffread_sequence_difference_rejects_atomic_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome, gff = write_single_fixture(root)
            executable = write_fake_gffread(root, matching=False)
            output = root / "output"
            with self.assertRaisesRegex(PrimaryAnnotationError, "gffread comparison failed"):
                standardize_primary_annotation(
                    genome,
                    gff,
                    output,
                    "sample",
                    gffread=str(executable),
                    require_gffread=True,
                )
            self.assertFalse(output.exists())
            self.assertFalse(list(root.glob(".output.staging.*")))

    def test_existing_output_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome, gff = write_single_fixture(root)
            output = root / "output"
            output.mkdir()
            sentinel = output / "keep.txt"
            sentinel.write_text("keep\n", encoding="utf-8")
            with self.assertRaisesRegex(PrimaryAnnotationError, "already exists"):
                standardize_primary_annotation(genome, gff, output, "sample", gffread="none")
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep\n")

    def test_gene_as_transcript_requires_explicit_opt_in(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome = root / "gene_only.fa"
            gff = root / "gene_only.gff3"
            output = root / "output"
            genome.write_text(">Chr01\nATGAAATAA\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "Chr01\ttest\tgene\t1\t9\t.\t+\t.\tID=ARU_gene1;biotype=protein_coding\n"
                "Chr01\ttest\tCDS\t1\t9\t.\t+\t0\tParent=ARU_gene1\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(PrimaryAnnotationError, "no transcript rows found"):
                standardize_primary_annotation(genome, gff, output, "sample", gffread="none")
            self.assertFalse(output.exists())
            self.assertFalse(list(root.glob(".output.staging.*")))

    def test_gene_as_transcript_preserves_ids_audits_noncoding_pseudogene_and_passes_gffread(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome = root / "gene_only.fa"
            gff = root / "gene_only.gff3"
            output = root / "output"
            genome.write_text(">Chr01\nATGAAATAA" + "N" * 21 + "\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "Chr01\ttest\tgene\t1\t9\t.\t+\t.\tID=ARU_gene1;biotype=protein_coding\n"
                "Chr01\ttest\tCDS\t1\t9\t.\t+\t0\tParent=ARU_gene1\n"
                "Chr01\ttest\tpseudogene\t20\t29\t.\t+\t.\tID=ARU_dead1\n",
                encoding="utf-8",
            )
            executable = write_fake_gffread(
                root,
                matching=True,
                reject_feature="gene",
                output_id="ARU_gene1",
            )
            result = standardize_primary_annotation(
                genome,
                gff,
                output,
                "sample",
                gene_as_transcript=True,
                gffread=str(executable),
                require_gffread=True,
            )
            self.assertEqual(result.source_gene_count, 2)
            self.assertEqual(result.selected_gene_count, 1)
            self.assertEqual(result.invalid_coding_gene_count, 0)
            self.assertEqual(read_fasta(output / "sample.protein.faa"), {"ARU_gene1": "MK"})

            with (output / "sample.primary_isoforms.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                selected = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(selected["gene_id"], "ARU_gene1")
            self.assertEqual(selected["transcript_id"], "ARU_gene1")
            self.assertIn("gene_as_transcript", selected["QC_flags"])

            with (output / "sample.gene_audit.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                gene_audit = {
                    row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(gene_audit["ARU_gene1"]["status"], "selected")
            self.assertEqual(gene_audit["ARU_dead1"]["status"], "no_CDS_transcript")
            self.assertEqual(gene_audit["ARU_dead1"]["invalid_coding_transcript_count"], "0")

            with (output / "sample.transcript_audit.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                transcript_audit = {
                    row["transcript_id"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(transcript_audit["ARU_dead1"]["validation_status"], "no_CDS")
            self.assertEqual(transcript_audit["ARU_dead1"]["QC_flags"], "gene_as_transcript")
            self.assertEqual(
                transcript_audit["ARU_dead1"]["disposition"],
                "noncoding_or_no_CDS",
            )

            primary_gff = (output / "sample.primary.gff3").read_text(encoding="utf-8")
            feature_rows = [line.split("\t") for line in primary_gff.splitlines() if line and not line.startswith("#")]
            self.assertEqual([row[2] for row in feature_rows], ["mRNA", "CDS"])
            self.assertIn("ID=ARU_gene1;gene_as_transcript=true", feature_rows[0][8])
            self.assertEqual(feature_rows[1][8], "Parent=ARU_gene1")

            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["policy"]["gene_as_transcript_requested"])
            self.assertEqual(manifest["policy"]["annotation_graph_mode"], "gene_as_transcript")
            self.assertEqual(manifest["counts"]["source_transcript_rows"], 0)
            self.assertEqual(manifest["counts"]["synthesized_gene_as_transcripts"], 2)
            self.assertEqual(manifest["counts"]["selected_transcripts"], 1)
            self.assertEqual(manifest["publication_gate"], "PASS")
            self.assertEqual(
                manifest["gffread_comparison"]["comparison_annotation_scope"],
                "selected_primary_GFF3_top_level_mRNA_CDS_only_gene_as_transcript",
            )
            self.assertTrue(manifest["gffread_comparison"]["exact_selected_CDS_ID_set"])
            self.assertTrue(manifest["gffread_comparison"]["exact_selected_protein_ID_set"])

    def test_gene_as_transcript_rejects_mixed_declared_and_direct_parent_graphs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome = root / "mixed.fa"
            gff = root / "mixed.gff3"
            output = root / "output"
            genome.write_text(">Chr01\nATGAAATAANATGAAATAA\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "Chr01\ttest\tgene\t1\t9\t.\t+\t.\tID=gene_direct\n"
                "Chr01\ttest\tCDS\t1\t9\t.\t+\t0\tParent=gene_direct\n"
                "Chr01\ttest\tgene\t11\t19\t.\t+\t.\tID=gene_standard\n"
                "Chr01\ttest\tmRNA\t11\t19\t.\t+\t.\tID=tx_standard;Parent=gene_standard\n"
                "Chr01\ttest\tCDS\t11\t19\t.\t+\t0\tParent=tx_standard\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                PrimaryAnnotationError,
                "CDS Parent IDs do not identify declared transcripts: gene_direct",
            ):
                standardize_primary_annotation(
                    genome,
                    gff,
                    output,
                    "sample",
                    gene_as_transcript=True,
                    gffread="none",
                )
            self.assertFalse(output.exists())

    def test_gene_as_transcript_option_does_not_rewrite_a_standard_graph(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome, gff = write_single_fixture(root)
            output = root / "output"
            standardize_primary_annotation(
                genome,
                gff,
                output,
                "sample",
                gene_as_transcript=True,
                gffread="none",
            )
            self.assertEqual(list(read_fasta(output / "sample.protein.faa")), ["tx1"])
            primary_gff = (output / "sample.primary.gff3").read_text(encoding="utf-8")
            self.assertEqual(
                [
                    line.split("\t")[2]
                    for line in primary_gff.splitlines()
                    if line and not line.startswith("#")
                ],
                ["gene", "mRNA", "CDS"],
            )
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["policy"]["gene_as_transcript_requested"])
            self.assertEqual(manifest["policy"]["annotation_graph_mode"], "declared_transcripts")
            self.assertEqual(manifest["counts"]["source_transcript_rows"], 1)
            self.assertEqual(manifest["counts"]["synthesized_gene_as_transcripts"], 0)

    def test_gene_as_transcript_rejects_multiple_or_undeclared_gene_parents(self) -> None:
        for label, parent_value, expected in (
            ("multiple", "gene1,gene2", "must have exactly one declared gene Parent"),
            ("undeclared", "missing_gene", "does not identify a declared gene-level feature"),
        ):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                genome = root / "gene_only.fa"
                gff = root / "gene_only.gff3"
                output = root / "output"
                genome.write_text(">Chr01\nATGAAATAA\n", encoding="utf-8")
                gff.write_text(
                    "##gff-version 3\n"
                    "Chr01\ttest\tgene\t1\t9\t.\t+\t.\tID=gene1\n"
                    "Chr01\ttest\tgene\t1\t9\t.\t+\t.\tID=gene2\n"
                    f"Chr01\ttest\tCDS\t1\t9\t.\t+\t0\tParent={parent_value}\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(PrimaryAnnotationError, expected):
                    standardize_primary_annotation(
                        genome,
                        gff,
                        output,
                        "sample",
                        gene_as_transcript=True,
                        gffread="none",
                    )
                self.assertFalse(output.exists())

    def test_gene_as_transcript_rejects_cds_attached_to_pseudogene(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome = root / "pseudogene.fa"
            gff = root / "pseudogene.gff3"
            output = root / "output"
            genome.write_text(">Chr01\nATGAAATAA\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "Chr01\ttest\tpseudogene\t1\t9\t.\t+\t.\tID=dead1\n"
                "Chr01\ttest\tCDS\t1\t9\t.\t+\t0\tParent=dead1\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                PrimaryAnnotationError,
                "declared as a pseudogene",
            ):
                standardize_primary_annotation(
                    genome,
                    gff,
                    output,
                    "sample",
                    gene_as_transcript=True,
                    gffread="none",
                )
            self.assertFalse(output.exists())

    def test_gene_as_transcript_rejects_transcript_level_canonical_rules(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome = root / "gene_only.fa"
            gff = root / "gene_only.gff3"
            output = root / "output"
            genome.write_text(">Chr01\nATGAAATAA\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "Chr01\ttest\tgene\t1\t9\t.\t+\t.\tID=gene1;tag=canonical\n"
                "Chr01\ttest\tCDS\t1\t9\t.\t+\t0\tParent=gene1\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                PrimaryAnnotationError,
                "Canonical transcript-tag rules cannot be applied",
            ):
                standardize_primary_annotation(
                    genome,
                    gff,
                    output,
                    "sample",
                    canonical_rules=(parse_canonical_rule("tag=canonical"),),
                    gene_as_transcript=True,
                    gffread="none",
                )
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
