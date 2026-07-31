"""Tests for the completed-batch to primary NLR summary adapter."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "nlr" / "prepare_primary_nlr_summary_inputs.py"
JAR = "1" * 64
MOTIFS = "2" * 64
STORE = "3" * 64


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PreparePrimaryNlrSummaryInputsTest(unittest.TestCase):
    def make_fixture(self, root: Path) -> dict[str, Path]:
        bundle = root / "bundle"
        bundle.mkdir()
        fasta_text = {
            "ref.fa": ">g1\nAAA\n>g2\nCCC\n>g3\nGGG\n>g4\nTTT\n",
            "u1.fa": ">Chr01\nAAAA\n",
            "u2.fa": ">Chr01\nCCCC\n",
        }
        for name, text in fasta_text.items():
            (bundle / name).write_text(text, encoding="utf-8")
        selected_rows = [
            {
                "sample_id": "clem_scandens_reference", "species": "Clematoclethra scandens",
                "ploidy": "n/a", "analysis_role": "reference_callable",
                "input_scope": "reference_transcript_cds", "relative_fasta": "ref.fa",
                "expected_fasta_records": 4,
            },
            {
                "sample_id": "u1", "species": "Actinidia arguta", "ploidy": "4x",
                "analysis_role": "target_repertoire", "input_scope": "whole_genome",
                "relative_fasta": "u1.fa", "expected_fasta_records": 1,
            },
            {
                "sample_id": "u2", "species": "Actinidia deliciosa", "ploidy": "6x",
                "analysis_role": "target_repertoire", "input_scope": "whole_genome",
                "relative_fasta": "u2.fa", "expected_fasta_records": 1,
            },
        ]
        runner = bundle / "nlr_annotator_inputs.tsv"
        selected_fields = [
            "sample_id", "species", "ploidy", "analysis_role", "input_scope",
            "relative_fasta", "expected_fasta_records",
        ]
        write_tsv(runner, selected_fields, selected_rows)
        output_rows = []
        sample_by_name = {
            "clem_scandens_reference": ("ref.fa", 4), "u1": ("u1.fa", 1), "u2": ("u2.fa", 1)
        }
        for sample, (name, records) in sample_by_name.items():
            path = bundle / name
            output_rows.append(
                {
                    "sample_id": sample, "basename": name, "bytes": path.stat().st_size,
                    "sha256": digest(path), "fasta_records": records, "total_bases": "",
                }
            )
        output_rows.append(
            {
                "sample_id": "bundle", "basename": runner.name, "bytes": runner.stat().st_size,
                "sha256": digest(runner),
            }
        )
        manifest = {
            "schema_version": 1, "status": "PASS", "workflow": "plain_fasta_nlr_input_bundle",
            "reference_count": 1, "target_count": 2, "outputs": output_rows,
        }
        (bundle / "run_manifest.json").write_text(
            json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_tsv(
            bundle / "checksums.tsv", ["file", "bytes", "sha256"],
            [
                {"file": row["basename"], "bytes": row["bytes"], "sha256": row["sha256"]}
                for row in output_rows
            ],
        )

        nlr_root = root / "nlr"
        nlr_root.mkdir()
        write_tsv(nlr_root / "selected_inputs.tsv", selected_fields, selected_rows)
        write_tsv(
            nlr_root / "batch_metadata.tsv", ["key", "value"],
            [
                {"key": "selected_inputs", "value": 3},
                {"key": "nlr_annotator_jar_sha256", "value": JAR},
                {"key": "motifs_sha256", "value": MOTIFS},
                {"key": "store_sha256", "value": STORE},
                {"key": "completion_status", "value": "complete"},
            ],
        )
        calls = {
            "clem_scandens_reference": [
                "g1\tL1\tNBARC\t1\t2\t+\tcomplete\n",
                "g2\tL2\tTIR-NBARC\t1\t2\t+\tcomplete\n",
                "g3\tL3\tCC-NBARC\t1\t2\t+\tcomplete\n",
            ],
            "u1": [
                "Chr01\tL1\tNBARC\t1\t2\t+\tcomplete\n",
                "Chr01\tL2\tTIR-NBARC\t4\t5\t+\tcomplete\n",
            ],
            "u2": ["Chr01\tL1\tNBARC\t1\t2\t+\tcomplete\n"],
        }
        selected = {row["sample_id"]: row for row in selected_rows}
        for sample, lines in calls.items():
            directory = nlr_root / sample
            directory.mkdir()
            (directory / "nlr_calls.txt").write_text("".join(lines), encoding="utf-8")
            (directory / "nlr_loci.gff").write_text("##gff-version 3\n", encoding="utf-8")
            (directory / "stdout.log").write_text("ok\n", encoding="utf-8")
            (directory / "stderr.log").write_text("", encoding="utf-8")
            parsed = [line.rstrip("\n").split("\t") for line in lines]
            sample_row = selected[sample]
            fasta_name = sample_row["relative_fasta"]
            metadata = {
                "sample_id": sample,
                "species": sample_row["species"],
                "ploidy": sample_row["ploidy"],
                "analysis_role": sample_row["analysis_role"],
                "input_scope": sample_row["input_scope"],
                "input_fasta": str(bundle / fasta_name),
                "input_fasta_sha256": digest(bundle / fasta_name),
                "input_fasta_records": str(sample_row["expected_fasta_records"]),
                "nlr_output_rows": str(len(parsed)),
                "nlr_output_sequence_ids": str(len({row[0] for row in parsed})),
                "nlr_output_locus_ids": str(len({row[1] for row in parsed})),
                "configured_nlr_worker_threads": "2",
                "jvm_processor_cap": "2",
                "completion_status": "complete",
                "nlr_annotator_jar_sha256": JAR,
                "motifs_sha256": MOTIFS,
                "store_sha256": STORE,
            }
            write_tsv(
                directory / "run_metadata.tsv", ["key", "value"],
                [{"key": key, "value": value} for key, value in metadata.items()],
            )
            write_tsv(
                directory / "output_checksums.tsv", ["path", "sha256"],
                [
                    {"path": name, "sha256": digest(directory / name)}
                    for name in ("nlr_calls.txt", "nlr_loci.gff", "stdout.log", "stderr.log")
                ],
            )

        metadata = root / "units.tsv"
        write_tsv(
            metadata,
            [
                "assembly_unit_id", "biological_species", "haplotype_or_subgenome",
                "assembly_scope", "include", "analysis_cohort",
            ],
            [
                {
                    "assembly_unit_id": "u1", "biological_species": "Actinidia arguta",
                    "haplotype_or_subgenome": "A", "assembly_scope": "29_chromosomes",
                    "include": "true", "analysis_cohort": "primary_test",
                },
                {
                    "assembly_unit_id": "u2", "biological_species": "Actinidia deliciosa",
                    "haplotype_or_subgenome": "A", "assembly_scope": "29_chromosomes",
                    "include": "true", "analysis_cohort": "primary_test",
                },
            ],
        )
        shared = root / "shared.tsv"
        write_tsv(shared, ["reference_gene_id"], [{"reference_gene_id": "g3"}])
        losses = root / "primary_complete_loss_matrix.tsv"
        states = {
            "u1": {
                "g1": ("deleted", "true"), "g2": ("uncertain", "true"),
                "g3": ("deleted", "true"), "g4": ("retained", "true"),
            },
            "u2": {
                "g1": ("retained", "true"), "g2": ("uncertain", "false"),
                "g3": ("deleted", "true"), "g4": ("retained", "true"),
            },
        }
        loss_rows = []
        for unit, genes in states.items():
            for gene, (classification, callable_value) in genes.items():
                loss_rows.append(
                    {
                        "reference_gene_id": gene, "assembly_unit_id": unit,
                        "classification": classification, "callable": callable_value,
                        "evidence_source": "fixture", "primary_search_state": "fixture",
                    }
                )
        write_tsv(losses, list((
            "reference_gene_id", "assembly_unit_id", "classification", "callable",
            "evidence_source", "primary_search_state",
        )), loss_rows)
        return {
            "bundle": bundle, "nlr": nlr_root, "metadata": metadata,
            "shared": shared, "losses": losses,
        }

    def run_script(self, root: Path, fixture: dict[str, Path], name: str = "output"):
        output = root / name
        completed = subprocess.run(
            [
                sys.executable, str(SCRIPT), "--nlr-root", str(fixture["nlr"]),
                "--input-bundle", str(fixture["bundle"]), "--unit-metadata",
                str(fixture["metadata"]), "--loss-matrix", str(fixture["losses"]),
                "--shared-positive-genes", str(fixture["shared"]), "--output-dir", str(output),
                "--expected-units", "2", "--expected-reference-genes", "4",
                "--expected-shared-positive", "1", "--expected-worker-threads", "2",
                "--expected-jar-sha256", JAR, "--expected-motifs-sha256", MOTIFS,
                "--expected-store-sha256", STORE,
            ],
            cwd=ROOT, capture_output=True, text=True, check=False,
        )
        return completed, output

    def test_builds_nonshared_resolved_catalog_and_repertoire_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = self.make_fixture(root)
            completed, output = self.run_script(root, fixture)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            validation = json.loads((output / "validation.json").read_text(encoding="utf-8"))
            self.assertEqual(validation["status"], "PASS_PRIMARY_NLR_SUMMARY_INPUTS")
            self.assertEqual(validation["reference_nlr_gene_count"], 3)
            self.assertEqual(validation["shared_reference_nlr_gene_count_excluded"], 1)
            self.assertEqual(validation["nonshared_reference_nlr_gene_count"], 2)
            self.assertEqual(validation["callable_reference_nlr_denominator_sum"], 2)
            self.assertEqual(validation["positive_reference_nlr_loss_call_count"], 1)

            denominators = read_tsv(output / "callable_reference_nlr_denominators.tsv")
            self.assertEqual(
                {(row["assembly_unit_id"], row["reference_nlr_id"]) for row in denominators},
                {("u1", "g1"), ("u2", "g1")},
            )
            positives = read_tsv(output / "positive_reference_nlr_loss_calls.tsv")
            self.assertEqual(
                [(row["assembly_unit_id"], row["reference_nlr_id"]) for row in positives],
                [("u1", "g1")],
            )
            repertoires = read_tsv(output / "repertoire_counts.tsv")
            self.assertEqual([row["total_nlr_count"] for row in repertoires], ["2", "1"])
            universe = read_tsv(output / "reference_nlr_universe.tsv")
            shared = next(row for row in universe if row["reference_nlr_id"] == "g3")
            self.assertEqual(shared["included_in_nonshared_analysis"], "false")
            self.assertEqual(shared["exclusion_reason"], "shared_positive_complete")
            output_checksums = read_tsv(output / "output_checksums.tsv")
            self.assertEqual(len(output_checksums), len(list(output.iterdir())) - 1)

    def test_sample_output_checksum_mismatch_fails_without_publishing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = self.make_fixture(root)
            with (fixture["nlr"] / "u1" / "nlr_calls.txt").open("a", encoding="utf-8") as handle:
                handle.write("Chr01\tL9\tNBARC\t1\t2\t+\tcomplete\n")
            completed, output = self.run_script(root, fixture)
            self.assertEqual(completed.returncode, 2)
            self.assertIn("checksum mismatch", completed.stderr)
            self.assertFalse(output.exists())

    def test_duplicate_loss_grid_pair_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = self.make_fixture(root)
            with fixture["losses"].open("a", encoding="utf-8") as handle:
                handle.write("g1\tu1\tdeleted\ttrue\tfixture\tfixture\n")
            completed, output = self.run_script(root, fixture)
            self.assertEqual(completed.returncode, 2)
            self.assertIn("Duplicate primary loss pair", completed.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
