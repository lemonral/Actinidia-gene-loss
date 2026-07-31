"""Focused tests for exact publisher-primary protein remapping."""

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

from geneloss_repro.publisher_protein_remap import (
    PublisherProteinRemapError,
    remap_publisher_primary_proteins,
)


def write_selected(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text(
        "".join(f">{identifier}\n{sequence}\n" for identifier, sequence in records),
        encoding="utf-8",
    )


def write_publisher(
    path: Path,
    records: list[tuple[str, str, str, str]],
    *,
    compressed: bool = False,
) -> None:
    text = "".join(
        f">{protein_id}\tmRNA={mrna}\tOriID={transcript_id}\tNote=publisher record\n"
        f"{sequence}\n"
        for protein_id, transcript_id, mrna, sequence in records
    )
    if compressed:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def write_gff(path: Path, rows: list[tuple[str, str, str]]) -> None:
    lines = ["##gff-version 3"]
    for index, (transcript_id, mrna, protein_id) in enumerate(rows, start=1):
        start = index * 100
        lines.extend(
            [
                f"Chr1\ttest\tgene\t{start}\t{start + 29}\t.\t+\t.\tID=g{index}",
                f"Chr1\ttest\tmRNA\t{start}\t{start + 29}\t.\t+\t.\t"
                f"ID={transcript_id};Accession={mrna};Parent=g{index}",
                f"Chr1\ttest\tCDS\t{start}\t{start + 29}\t.\t+\t0\t"
                f"ID={transcript_id}.CDS1;Parent={transcript_id};"
                f"Protein_Accession={protein_id}",
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_first_token_publisher(
    path: Path, records: list[tuple[str, str]]
) -> None:
    path.write_text(
        "".join(
            f">{protein_id} publisher description\n{sequence}\n"
            for protein_id, sequence in records
        ),
        encoding="utf-8",
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            identifier = raw[1:].split()[0]
            records[identifier] = ""
        elif identifier is not None:
            records[identifier] += raw
    return records


class PublisherPrimaryRemapTests(unittest.TestCase):
    def test_cli_help_works_without_pythonpath(self) -> None:
        script = (
            Path(__file__).parents[2]
            / "scripts"
            / "qc"
            / "remap_publisher_primary_proteins.py"
        )
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
        self.assertIn("--selected-primary-proteins", completed.stdout)
        self.assertIn("--protein-accession-attribute", completed.stdout)
        self.assertIn("--publisher-header-mode", completed.stdout)

    def test_ncbi_first_token_mode_uses_gff_mapping_and_exact_source_closure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("rna-t1", "MA")])
            gff.write_text(
                "##gff-version 3\n"
                "Chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=gene-g1\n"
                "Chr1\ttest\tmRNA\t1\t30\t.\t+\t.\t"
                "ID=rna-t1;orig_transcript_id=gnl|WGS|t1;Parent=gene-g1\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\t"
                "ID=cds-p1;Parent=rna-t1;protein_id=KAI000001.1\n"
                "Chr1\ttest\tgene\t101\t130\t.\t+\t.\tID=gene-g2\n"
                "Chr1\ttest\tmRNA\t101\t130\t.\t+\t.\t"
                "ID=rna-t2;orig_transcript_id=gnl|WGS|t2;Parent=gene-g2\n"
                "Chr1\ttest\tCDS\t101\t130\t.\t+\t0\t"
                "ID=cds-p2;Parent=rna-t2;protein_id=KAI000002.1\n",
                encoding="utf-8",
            )
            publisher.write_text(
                ">KAI000001.1 hypothetical protein one [Example species]\nMA\n"
                ">KAI000002.1 hypothetical protein two [Example species]\nMP\n",
                encoding="utf-8",
            )

            result = remap_publisher_primary_proteins(
                selected,
                gff,
                publisher,
                output,
                "ncbi_example",
                transcript_accession_attribute="orig_transcript_id",
                protein_accession_attribute="protein_id",
                publisher_header_mode="first_token",
            )

            self.assertEqual(read_fasta(result.output_protein_path), {"rna-t1": "MA"})
            summary = read_tsv(output / "ncbi_example.publisher_primary.summary.tsv")[0]
            self.assertEqual(summary["publisher_header_mode"], "first_token")
            self.assertEqual(summary["publisher_header_mapping_check"], "not_applicable")
            self.assertEqual(summary["exact_source_accession_closure"], "true")
            mapping = read_tsv(output / "ncbi_example.publisher_primary.mapping.tsv")[0]
            self.assertEqual(mapping["publisher_protein_id"], "KAI000001.1")
            self.assertEqual(mapping["publisher_header_transcript_id"], "")
            self.assertEqual(mapping["publisher_header_mRNA_accession"], "")

    def test_actinidiabase_normal_graph_explicitly_uses_self_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("Aruf1g000001.t1", "MA")])
            gff.write_text(
                "##gff-version 3\n"
                "Chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=Aruf1g000001\n"
                "Chr1\ttest\tmRNA\t1\t30\t.\t+\t.\t"
                "ID=Aruf1g000001.t1;Parent=Aruf1g000001\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\t"
                "ID=Aruf1g000001.t1:CDS;Parent=Aruf1g000001.t1\n",
                encoding="utf-8",
            )
            write_first_token_publisher(
                publisher, [("Aruf1g000001.t1", "MA")]
            )

            result = remap_publisher_primary_proteins(
                selected,
                gff,
                publisher,
                output,
                "act_rufa_actinidiabase_v1",
                transcript_accession_source="transcript_id",
                protein_accession_source="cds_parent",
                publisher_header_mode="first_token",
            )

            self.assertEqual(
                read_fasta(result.output_protein_path),
                {"Aruf1g000001.t1": "MA"},
            )
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(
                manifest["policy"]["annotation_graph_mode"],
                "declared_transcripts",
            )
            self.assertEqual(
                manifest["schema"]["transcript_accession_source"],
                "transcript_id",
            )
            self.assertEqual(
                manifest["schema"]["protein_accession_source"], "cds_parent"
            )

    def test_gene_as_transcript_parent_accession_passes_and_audits_all_models(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("Aru1.0ch01g00010.1", "MA")])
            gff.write_text(
                "##gff-version 3\n"
                "Chr1\ttest\tgene\t1\t60\t.\t+\t.\tID=Aru1.0ch01g00010.1\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\t"
                "ID=Aru1.0ch01g00010.1:CDS:1;Parent=Aru1.0ch01g00010.1\n"
                "Chr1\ttest\tCDS\t31\t60\t.\t+\t0\t"
                "ID=Aru1.0ch01g00010.1:CDS:2;Parent=Aru1.0ch01g00010.1\n"
                "Chr1\ttest\tgene\t101\t130\t.\t+\t.\tID=Aru1.0ch01g00020.1\n"
                "Chr1\ttest\tCDS\t101\t130\t.\t+\t0\t"
                "ID=Aru1.0ch01g00020.1:CDS:1;Parent=Aru1.0ch01g00020.1\n"
                "Chr1\ttest\tpseudogene\t201\t230\t.\t+\t.\tID=Aru_dead1\n",
                encoding="utf-8",
            )
            write_first_token_publisher(
                publisher,
                [
                    ("Aru1.0ch01g00010.1", "MA"),
                    ("Aru1.0ch01g00020.1", "MP"),
                ],
            )

            result = remap_publisher_primary_proteins(
                selected,
                gff,
                publisher,
                output,
                "act_rufa_aru_r1",
                gene_as_transcript=True,
                transcript_accession_source="transcript_id",
                protein_accession_source="cds_parent",
                publisher_header_mode="first_token",
            )

            self.assertEqual(
                read_fasta(result.output_protein_path),
                {"Aru1.0ch01g00010.1": "MA"},
            )
            summary = read_tsv(
                output / "act_rufa_aru_r1.publisher_primary.summary.tsv"
            )[0]
            self.assertEqual(summary["annotation_graph_mode"], "gene_as_transcript")
            self.assertEqual(summary["gene_as_transcript_requested"], "true")
            self.assertEqual(summary["gff_gene_record_count"], "3")
            self.assertEqual(summary["gff_source_transcript_row_count"], "0")
            self.assertEqual(summary["gff_coding_transcript_count"], "2")
            self.assertEqual(summary["gff_noncoding_model_count"], "1")
            inventory = read_tsv(
                output / "act_rufa_aru_r1.publisher_protein.source_inventory.tsv"
            )
            self.assertEqual(
                [row["disposition"] for row in inventory],
                ["SELECTED_AND_REMAPPED", "EXCLUDED_NONPRIMARY"],
            )
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertTrue(manifest["policy"]["gene_as_transcript_requested"])
            self.assertEqual(
                manifest["counts"]["GFF3_source_transcript_rows"], 0
            )

    def test_gene_as_transcript_declared_protein_attribute_passes_via_cli(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("gene-Acr_01g0000010", "MA")])
            gff.write_text(
                "##gff-version 3\n"
                "Chr1\ttest\tgene\t1\t60\t.\t+\t.\tID=gene-Acr_01g0000010\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\t"
                "ID=cds-GFY80192.1;Parent=gene-Acr_01g0000010;protein_id=GFY80192.1\n"
                "Chr1\ttest\tCDS\t31\t60\t.\t+\t0\t"
                "ID=cds-GFY80192.1;Parent=gene-Acr_01g0000010;protein_id=GFY80192.1\n",
                encoding="utf-8",
            )
            write_first_token_publisher(publisher, [("GFY80192.1", "MA")])
            script = (
                Path(__file__).parents[2]
                / "scripts"
                / "qc"
                / "remap_publisher_primary_proteins.py"
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--selected-primary-proteins",
                    str(selected),
                    "--gff",
                    str(gff),
                    "--publisher-proteins",
                    str(publisher),
                    "--sample-id",
                    "act_rufa_fuchu",
                    "--output-dir",
                    str(output),
                    "--gene-as-transcript",
                    "--transcript-accession-source",
                    "transcript_id",
                    "--protein-accession-attribute",
                    "protein_id",
                    "--publisher-header-mode",
                    "first_token",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(
                read_fasta(
                    output
                    / "act_rufa_fuchu.publisher_primary.remapped.protein.faa"
                ),
                {"gene-Acr_01g0000010": "MA"},
            )

    def test_gene_to_cds_graph_requires_explicit_opt_in(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("g1", "MA")])
            gff.write_text(
                "##gff-version 3\n"
                "Chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=g1;protein_id=P1\n",
                encoding="utf-8",
            )
            write_first_token_publisher(publisher, [("P1", "MA")])
            with self.assertRaisesRegex(
                PublisherProteinRemapError, "explicit --gene-as-transcript"
            ):
                remap_publisher_primary_proteins(
                    selected,
                    gff,
                    publisher,
                    output,
                    "act_test",
                    transcript_accession_source="transcript_id",
                    protein_accession_attribute="protein_id",
                    publisher_header_mode="first_token",
                )
            self.assertFalse(output.exists())

    def test_gene_as_transcript_rejects_mixed_graph(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("g1", "MA")])
            gff.write_text(
                "##gff-version 3\n"
                "Chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
                "Chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=t1;Parent=g1\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=g1;protein_id=P1\n",
                encoding="utf-8",
            )
            write_first_token_publisher(publisher, [("P1", "MA")])
            with self.assertRaisesRegex(
                PublisherProteinRemapError, "requires zero accepted transcript rows"
            ):
                remap_publisher_primary_proteins(
                    selected,
                    gff,
                    publisher,
                    output,
                    "act_test",
                    gene_as_transcript=True,
                    transcript_accession_source="transcript_id",
                    protein_accession_attribute="protein_id",
                    publisher_header_mode="first_token",
                )
            self.assertFalse(output.exists())

    def test_gene_as_transcript_rejects_structural_and_accession_ambiguity(self) -> None:
        cases = {
            "multiple_parent": (
                "Chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
                "Chr1\ttest\tgene\t31\t60\t.\t+\t.\tID=g2\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=g1,g2;protein_id=P1\n",
                "must have exactly one declared gene parent",
            ),
            "undeclared_parent": (
                "Chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=missing;protein_id=P1\n",
                "is not a declared gene-level ID",
            ),
            "pseudogene_cds": (
                "Chr1\ttest\tpseudogene\t1\t30\t.\t+\t.\tID=g1\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=g1;protein_id=P1\n",
                "attached to declared pseudogene",
            ),
            "missing_accession": (
                "Chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=g1\n",
                "must declare exactly one 'protein_id' accession",
            ),
            "multiple_accessions": (
                "Chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=g1;protein_id=P1,P2\n",
                "must declare exactly one 'protein_id' accession",
            ),
            "multipart_conflict": (
                "Chr1\ttest\tgene\t1\t60\t.\t+\t.\tID=g1\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=g1;protein_id=P1\n"
                "Chr1\ttest\tCDS\t31\t60\t.\t+\t0\tParent=g1;protein_id=P2\n",
                "transcripts map to multiple protein accessions",
            ),
            "nonselected_duplicate_mapping": (
                "Chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n"
                "Chr1\ttest\tCDS\t1\t30\t.\t+\t0\tParent=g1;protein_id=P1\n"
                "Chr1\ttest\tgene\t31\t60\t.\t+\t.\tID=g2\n"
                "Chr1\ttest\tCDS\t31\t60\t.\t+\t0\tParent=g2;protein_id=PX\n"
                "Chr1\ttest\tgene\t61\t90\t.\t+\t.\tID=g3\n"
                "Chr1\ttest\tCDS\t61\t90\t.\t+\t0\tParent=g3;protein_id=PX\n",
                "protein accessions map to multiple transcripts",
            ),
        }
        for label, (body, expected) in cases.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                selected = root / "selected.faa"
                gff = root / "annotation.gff3"
                publisher = root / "publisher.faa"
                output = root / "remap"
                write_selected(selected, [("g1", "MA")])
                gff.write_text("##gff-version 3\n" + body, encoding="utf-8")
                write_first_token_publisher(
                    publisher, [("P1", "MA"), ("P2", "MP"), ("PX", "MK")]
                )
                with self.assertRaisesRegex(PublisherProteinRemapError, expected):
                    remap_publisher_primary_proteins(
                        selected,
                        gff,
                        publisher,
                        output,
                        "act_test",
                        gene_as_transcript=True,
                        transcript_accession_source="transcript_id",
                        protein_accession_attribute="protein_id",
                        publisher_header_mode="first_token",
                    )
                self.assertFalse(output.exists())

    def test_pass_bundle_has_exact_closure_and_preserves_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.data"
            output = root / "remap"
            write_selected(selected, [("t1", "MDT"), ("t3", "MA")])
            write_gff(
                gff,
                [
                    ("t1", "mrna1", "protein1"),
                    ("t2", "mrna2", "protein2"),
                    ("t3", "mrna3", "protein3"),
                ],
            )
            # Gzip is detected by content, not suffix.
            write_publisher(
                publisher,
                [
                    ("protein1", "t1", "mrna1", "MXX*"),
                    ("protein2", "t2", "mrna2", "MPE"),
                    ("protein3", "t3", "mrna3", "MA"),
                ],
                compressed=True,
            )

            result = remap_publisher_primary_proteins(
                selected, gff, publisher, output, "act_deliciosa_adm_A"
            )
            self.assertEqual(result.source_publisher_record_count, 3)
            self.assertEqual(result.selected_primary_record_count, 2)
            self.assertEqual(result.excluded_nonprimary_record_count, 1)
            self.assertEqual(
                read_fasta(result.output_protein_path), {"t1": "MXX*", "t3": "MA"}
            )

            summary = read_tsv(
                output / "act_deliciosa_adm_A.publisher_primary.summary.tsv"
            )[0]
            self.assertEqual(summary["publication_gate"], "PASS")
            self.assertEqual(summary["exact_source_accession_closure"], "true")
            self.assertEqual(summary["exact_selected_primary_closure"], "true")
            self.assertEqual(summary["sequence_preservation"], "true")
            self.assertEqual(summary["excluded_nonprimary_record_count"], "1")

            mapping = read_tsv(
                output / "act_deliciosa_adm_A.publisher_primary.mapping.tsv"
            )
            self.assertEqual([row["selected_transcript_id"] for row in mapping], ["t1", "t3"])
            self.assertTrue(
                all(
                    row["source_sequence_sha256"] == row["output_sequence_sha256"]
                    for row in mapping
                )
            )
            inventory = read_tsv(
                output / "act_deliciosa_adm_A.publisher_protein.source_inventory.tsv"
            )
            self.assertEqual(
                [row["disposition"] for row in inventory],
                ["SELECTED_AND_REMAPPED", "EXCLUDED_NONPRIMARY", "SELECTED_AND_REMAPPED"],
            )

            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["publication_gate"], "PASS")
            self.assertEqual(manifest["counts"]["selected_primary_transcripts"], 2)
            checksums = read_tsv(output / "checksums.tsv")
            self.assertEqual(
                {row["file"] for row in checksums},
                {path.name for path in output.iterdir() if path.name != "checksums.tsv"},
            )
            for row in checksums:
                payload = (output / row["file"]).read_bytes()
                self.assertEqual(len(payload), int(row["bytes"]))
                self.assertEqual(hashlib.sha256(payload).hexdigest(), row["sha256"])
            text_outputs = "".join(
                path.read_text(encoding="utf-8")
                for path in output.iterdir()
                if path.is_file()
            )
            self.assertNotIn(str(root), text_outputs)

    def test_missing_selected_mapping_fails_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("absent", "MA")])
            write_gff(gff, [("t1", "mrna1", "protein1")])
            write_publisher(publisher, [("protein1", "t1", "mrna1", "MA")])
            with self.assertRaisesRegex(
                PublisherProteinRemapError, "Selected primary transcript IDs lack exact"
            ):
                remap_publisher_primary_proteins(
                    selected, gff, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())
            self.assertFalse(list(root.glob(".remap.staging.*")))

    def test_publisher_extra_or_gff_missing_accession_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("t1", "MA")])
            write_gff(gff, [("t1", "mrna1", "protein1")])
            write_publisher(
                publisher,
                [
                    ("protein1", "t1", "mrna1", "MA"),
                    ("extra", "extra_t", "extra_m", "MP"),
                ],
            )
            with self.assertRaisesRegex(
                PublisherProteinRemapError,
                "accession sets are not exactly equal.*publisher_only=1",
            ):
                remap_publisher_primary_proteins(
                    selected, gff, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())

    def test_gff_protein_without_publisher_record_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("t1", "MA")])
            write_gff(
                gff,
                [("t1", "mrna1", "protein1"), ("t2", "mrna2", "protein2")],
            )
            write_publisher(publisher, [("protein1", "t1", "mrna1", "MA")])
            with self.assertRaisesRegex(
                PublisherProteinRemapError,
                "accession sets are not exactly equal: gff_only=1",
            ):
                remap_publisher_primary_proteins(
                    selected, gff, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())

    def test_ambiguous_gff_transcript_to_protein_map_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("t1", "MA")])
            write_gff(gff, [("t1", "mrna1", "protein1")])
            with gff.open("a", encoding="utf-8") as handle:
                handle.write(
                    "Chr1\ttest\tCDS\t110\t119\t.\t+\t0\t"
                    "ID=t1.CDS2;Parent=t1;Protein_Accession=protein2\n"
                )
            write_publisher(
                publisher,
                [
                    ("protein1", "t1", "mrna1", "MA"),
                    ("protein2", "t1", "mrna1", "MP"),
                ],
            )
            with self.assertRaisesRegex(
                PublisherProteinRemapError, "transcripts map to multiple protein accessions"
            ):
                remap_publisher_primary_proteins(
                    selected, gff, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())

    def test_one_gff_protein_accession_cannot_map_two_transcripts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("t1", "MA")])
            write_gff(
                gff,
                [("t1", "mrna1", "protein1"), ("t2", "mrna2", "protein1")],
            )
            write_publisher(publisher, [("protein1", "t1", "mrna1", "MA")])
            with self.assertRaisesRegex(
                PublisherProteinRemapError, "protein accessions map to multiple transcripts"
            ):
                remap_publisher_primary_proteins(
                    selected, gff, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())

    def test_publisher_header_must_match_gff_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("t1", "MA")])
            write_gff(gff, [("t1", "mrna1", "protein1")])
            write_publisher(publisher, [("protein1", "wrong", "mrna1", "MA")])
            with self.assertRaisesRegex(
                PublisherProteinRemapError, "header mapping disagrees with GFF3"
            ):
                remap_publisher_primary_proteins(
                    selected, gff, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())

    def test_duplicate_selected_identifier_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            selected.write_text(">t1\nMA\n>t1 duplicate\nMA\n", encoding="utf-8")
            write_gff(gff, [("t1", "mrna1", "protein1")])
            write_publisher(publisher, [("protein1", "t1", "mrna1", "MA")])
            with self.assertRaisesRegex(PublisherProteinRemapError, "repeats identifier"):
                remap_publisher_primary_proteins(
                    selected, gff, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())

    def test_existing_output_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected.faa"
            gff = root / "annotation.gff3"
            publisher = root / "publisher.faa"
            output = root / "remap"
            write_selected(selected, [("t1", "MA")])
            write_gff(gff, [("t1", "mrna1", "protein1")])
            write_publisher(publisher, [("protein1", "t1", "mrna1", "MA")])
            output.mkdir()
            sentinel = output / "keep.txt"
            sentinel.write_text("keep\n", encoding="utf-8")
            with self.assertRaisesRegex(PublisherProteinRemapError, "already exists"):
                remap_publisher_primary_proteins(
                    selected, gff, publisher, output, "act_test"
                )
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep\n")


if __name__ == "__main__":
    unittest.main()
