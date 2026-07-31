"""Focused tests for paired assembly-scope reconciliation."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).parents[2]
    / "scripts"
    / "qc"
    / "audit_paired_assembly_scope.py"
)
SPEC = importlib.util.spec_from_file_location("audit_paired_assembly_scope", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text(
        "".join(f">{record_id}\n{sequence}\n" for record_id, sequence in records),
        encoding="utf-8",
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ScopeFixture:
    """A two-unit cohort spanning authoritative and candidate evidence."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.manifest = root / "paired_scope.tsv"
        self.output = root / "audit"
        self.official_analysis = root / "official.analysis.fa"
        self.official_deposit = root / "official.full.fa"
        self.official_report = root / "official.assembly_report.txt"
        self.candidate_analysis = root / "candidate.analysis.fa"
        self.candidate_deposit = root / "candidate.full.fa"
        self._write_inputs()
        self._write_manifest()

    def _write_inputs(self) -> None:
        write_fasta(self.official_analysis, [("chr1", "ACGT")])
        write_fasta(
            self.official_deposit,
            [("chr1", "ACGT"), ("unplaced_1", "TTAA")],
        )
        self.official_report.write_text(
            "# Sequence-Name\tSequence-Role\tSequence-Length\n"
            "chr1\tassembled-molecule\t4\n"
            "unplaced_1\tunplaced-scaffold\t4\n",
            encoding="utf-8",
        )
        write_fasta(self.candidate_analysis, [("legacy_chr", "GGCC")])
        write_fasta(
            self.candidate_deposit,
            [("renamed_chr", "GGCC"), ("candidate_extra", "AAAA")],
        )

    def _write_manifest(self) -> None:
        rows = [
            {
                "assembly_unit_id": "act_test_official",
                "biological_species": "Actinidia testensis",
                "pair_class": "official_same_accession_full_deposit",
                "analysis_accession": "GCA_TEST.1",
                "deposit_accession": "GCA_TEST.1",
                "analysis_fasta": self.official_analysis.name,
                "deposit_fasta": self.official_deposit.name,
                "assembly_report": self.official_report.name,
                "analysis_scope": "chromosome_only",
                "deposit_scope": "full_deposit",
                "evidence_note": "Synthetic official-scope fixture.",
            },
            {
                "assembly_unit_id": "act_test_candidate",
                "biological_species": "Actinidia testensis",
                "pair_class": "combined_haplotype_legacy_candidate",
                "analysis_accession": "LEGACY_TEST",
                "deposit_accession": "CANDIDATE_TEST",
                "analysis_fasta": self.candidate_analysis.name,
                "deposit_fasta": self.candidate_deposit.name,
                "assembly_report": "",
                "analysis_scope": "legacy_partition",
                "deposit_scope": "combined_candidate",
                "evidence_note": "Synthetic candidate-scope fixture.",
            },
        ]
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=MODULE.MANIFEST_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)


class PairedAssemblyScopeAuditTest(unittest.TestCase):
    def test_manifest_defines_cohort_and_optional_exact_count_is_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = ScopeFixture(Path(temporary))
            rows = MODULE.load_manifest(fixture.manifest)
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0].assembly_unit_id, "act_test_official")

            with self.assertRaisesRegex(
                MODULE.ScopeAuditError,
                "exactly 3 assembly-unit rows; observed 2",
            ):
                MODULE.load_manifest(fixture.manifest, expected_count=3)
            with self.assertRaisesRegex(
                MODULE.ScopeAuditError,
                "expected-assembly-unit-count must be positive",
            ):
                MODULE.load_manifest(fixture.manifest, expected_count=0)

    def test_full_audit_keeps_authoritative_and_candidate_scope_distinct(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fixture = ScopeFixture(Path(temporary))
            rows = MODULE.load_manifest(fixture.manifest)
            metadata = MODULE.publish(rows, fixture.manifest, fixture.output)

            self.assertEqual(metadata["manifest"]["assembly_unit_rows"], 2)
            summary = {
                row["assembly_unit_id"]: row
                for row in read_tsv(fixture.output / "paired_scope_summary.tsv")
            }
            official = summary["act_test_official"]
            self.assertEqual(official["analysis_bp_reconciled_percent"], "100.000000")
            self.assertEqual(official["official_unplaced_record_count"], "1")
            self.assertEqual(
                official["record_reconciliation_status"],
                "official_analysis_subset_reconciled",
            )

            candidate = summary["act_test_candidate"]
            self.assertEqual(candidate["sequence_hash_only_match_count"], "1")
            self.assertEqual(
                candidate["deposit_only_interpretation"],
                "candidate_only_scope_unknown",
            )
            deposit_only = {
                (row["assembly_unit_id"], row["deposit_record_id"]): row
                for row in read_tsv(fixture.output / "deposit_only_records.tsv")
            }
            self.assertEqual(
                deposit_only[("act_test_official", "unplaced_1")][
                    "scope_interpretation"
                ],
                "official_unplaced_scaffold",
            )
            self.assertEqual(
                deposit_only[("act_test_candidate", "candidate_extra")][
                    "scope_interpretation"
                ],
                "candidate_only_scope_unknown",
            )

            input_rows = read_tsv(fixture.output / "input_checksums.tsv")
            self.assertTrue(
                all("consumer_assembly_units" in row for row in input_rows)
            )
            persisted = json.loads(
                (fixture.output / "run_metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(persisted["manifest"]["assembly_unit_rows"], 2)

    def test_empty_manifest_is_rejected_when_count_is_derived(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            manifest = Path(temporary) / "empty.tsv"
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=MODULE.MANIFEST_COLUMNS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
            with self.assertRaisesRegex(
                MODULE.ScopeAuditError, "contains no assembly-unit rows"
            ):
                MODULE.load_manifest(manifest)


if __name__ == "__main__":
    unittest.main()
