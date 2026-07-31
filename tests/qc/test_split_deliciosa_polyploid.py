"""Focused tests for the fail-closed A. deliciosa A--F splitter."""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[2] / "scripts" / "qc" / "split_deliciosa_polyploid.py"
SPEC = importlib.util.spec_from_file_location("split_deliciosa_polyploid", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(f">{header}\n{sequence}\n" for header, sequence in records), encoding="utf-8")


class Fixture:
    def __init__(self, root: Path, scheme: str) -> None:
        self.root = root
        self.scheme = scheme
        self.bundle = "qinmei_test" if scheme == "qinmei_suffix" else "adm_test"
        self.mapping = root / "mapping.tsv"
        self.manifest = root / "resolved.tsv"
        self.output = root / "out"
        self.genome = root / "genome.fa"
        self.gff = root / "genes.gff3"
        self.cds = root / "cds.fa"
        self.protein = root / "protein.fa"
        self._write_mapping()
        self._write_assets()
        self._write_manifest()

    def _write_mapping(self) -> None:
        with self.mapping.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MODULE.MAPPING_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for index, label in enumerate("ABCDEF", start=1):
                writer.writerow(
                    {
                        "bundle_id": self.bundle,
                        "partition_scheme": self.scheme,
                        "partition_token": str(index),
                        "partition_label": label,
                        "assembly_unit_id": f"{self.bundle}_{label}",
                        "biological_species": "Actinidia deliciosa",
                        "individual_id": "test_individual",
                        "ploidy": "6x",
                        "expected_chromosome_count": "1",
                        "allow_unplaced_annotations": "true" if self.scheme == "qinmei_suffix" else "false",
                    }
                )

    def _write_assets(self) -> None:
        genome_records: list[tuple[str, str]] = []
        gff_lines = ["##gff-version 3\n"]
        cds_records: list[tuple[str, str]] = []
        protein_records: list[tuple[str, str]] = []
        for index, label in enumerate("ABCDEF", start=1):
            if self.scheme == "qinmei_suffix":
                seqid = f"chr1_{index}"
                feature_id = f"Ad_1_{index}g1"
                genome_header = seqid
                annotation_header = feature_id
            else:
                seqid = f"GWH_TEST_{index}"
                feature_id = f"Achdmh{index}c01g1"
                transcript_accession = f"GWHT_TEST_{index}"
                protein_accession = f"GWHP_TEST_{index}"
                genome_header = f"{seqid} Chromosome {label}1 OriSeqID=Chr01h{index}"
                cds_header = (
                    f"{transcript_accession} Protein={protein_accession} Position={seqid}:1-6 "
                    f"OriID={feature_id}.t1 OriSeqID=Chr01h{index}"
                )
                protein_header = (
                    f"{protein_accession} mRNA={transcript_accession} Position={seqid}:1-6 "
                    f"OriID={feature_id}.t1 OriSeqID=Chr01h{index}"
                )
            genome_records.append((genome_header, "ACGTAC"))
            if self.scheme == "qinmei_suffix":
                gff_lines.append(f"{seqid}\ttest\tmRNA\t1\t6\t.\t+\t.\tID={feature_id}\n")
                cds_records.append((annotation_header, "ATGGCC"))
                protein_records.append((annotation_header, "MA"))
            else:
                gff_lines.append(
                    f"{seqid}\ttest\tmRNA\t1\t6\t.\t+\t.\t"
                    f"ID={feature_id}.t1;Accession={transcript_accession};"
                    f"Protein_Accession={protein_accession}\n"
                )
                cds_records.append((cds_header, "ATGGCC"))
                protein_records.append((protein_header, "MA"))
        if self.scheme == "qinmei_suffix":
            gff_lines.append("scf1\ttest\tmRNA\t1\t3\t.\t+\t.\tID=Ad_scf1g1\n")
            cds_records.append(("Ad_scf1g1", "ATG"))
            protein_records.append(("Ad_scf1g1", "M"))
        write_fasta(self.genome, genome_records)
        self.gff.write_text("".join(gff_lines), encoding="utf-8")
        write_fasta(self.cds, cds_records)
        write_fasta(self.protein, protein_records)

    def _write_manifest(self) -> None:
        columns = list(MODULE.RESOLVED_COLUMNS)
        row = {
            "assembly_unit_id": self.bundle,
            "accession": "TEST0001",
            "genome": self.genome.name,
            "gff": self.gff.name,
            "cds": self.cds.name,
            "protein": self.protein.name,
            "genome_local_sha256": digest(self.genome),
            "gff_local_sha256": digest(self.gff),
            "cds_local_sha256": digest(self.cds),
            "protein_local_sha256": digest(self.protein),
        }
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(row)

    def run(self):
        return MODULE.run_partition(
            resolved_manifest=self.manifest,
            bundle_id=self.bundle,
            mapping=self.mapping,
            output_dir=self.output,
        )


class SplitDeliciosaTests(unittest.TestCase):
    def test_qinmei_split_passes_and_retains_unplaced_records(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary), "qinmei_suffix")
            result = fixture.run()
            self.assertEqual(result.status, "PASS")
            self.assertTrue((fixture.output / "resolved_assembly_units.tsv").is_file())
            with (fixture.output / "resolved_assembly_units.tsv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["haplotype_or_subgenome"] for row in rows], list("ABCDEF"))
            unassigned = fixture.output / "audit" / "unassigned" / f"{fixture.bundle}.protein.fa.gz"
            with gzip.open(unassigned, "rt", encoding="utf-8") as handle:
                self.assertIn(">Ad_scf1g1", handle.read())
            metadata = (fixture.output / "audit" / "run_metadata.json").read_text(encoding="utf-8")
            self.assertIn('"status": "PASS"', metadata)

    def test_adm_split_passes_with_agreeing_header_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary), "adm_header")
            result = fixture.run()
            self.assertEqual(result.status, "PASS")
            unit = f"{fixture.bundle}_F"
            genome = fixture.output / "assembly_units" / unit / f"{unit}.genome.fa.gz"
            with gzip.open(genome, "rt", encoding="utf-8") as handle:
                self.assertIn("Chromosome F1", handle.read())

    def test_checksum_mismatch_publishes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary), "qinmei_suffix")
            fixture.protein.write_text(">changed\nM\n", encoding="utf-8")
            with self.assertRaises(MODULE.PartitionInputError):
                fixture.run()
            self.assertFalse(fixture.output.exists())

    def test_unmatched_annotation_blocks_but_retains_audit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary), "qinmei_suffix")
            with fixture.protein.open("a", encoding="utf-8") as handle:
                handle.write(">mystery_gene\nM\n")
            fixture._write_manifest()
            result = fixture.run()
            self.assertEqual(result.status, "BLOCKED")
            self.assertFalse((fixture.output / "resolved_assembly_units.tsv").exists())
            unassigned = fixture.output / "audit" / "unassigned" / f"{fixture.bundle}.protein.fa.gz"
            with gzip.open(unassigned, "rt", encoding="utf-8") as handle:
                self.assertIn("mystery_gene", handle.read())
            issues = (fixture.output / "audit" / "validation_issues.tsv").read_text(encoding="utf-8")
            self.assertIn("unmatched_protein_record", issues)

    def test_exact_cds_protein_join_mismatch_blocks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary), "adm_header")
            text = fixture.protein.read_text(encoding="utf-8")
            fixture.protein.write_text(
                text.replace("mRNA=GWHT_TEST_1", "mRNA=GWHT_WRONG_1"), encoding="utf-8"
            )
            fixture._write_manifest()
            result = fixture.run()
            self.assertEqual(result.status, "BLOCKED")
            joins = (fixture.output / "audit" / "id_join_exceptions.tsv").read_text(encoding="utf-8")
            self.assertIn("adm_nonreciprocal_cds_protein_link", joins)

    def test_conflicting_adm_header_blocks_and_is_retained(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = Fixture(Path(temporary), "adm_header")
            text = fixture.genome.read_text(encoding="utf-8")
            fixture.genome.write_text(text.replace("Chromosome A1 OriSeqID=Chr01h1", "Chromosome A1 OriSeqID=Chr01h2"), encoding="utf-8")
            fixture._write_manifest()
            result = fixture.run()
            self.assertEqual(result.status, "BLOCKED")
            fasta_audit = (fixture.output / "audit" / "fasta_records.tsv").read_text(encoding="utf-8")
            self.assertIn("conflicting", fasta_audit)


if __name__ == "__main__":
    unittest.main()
