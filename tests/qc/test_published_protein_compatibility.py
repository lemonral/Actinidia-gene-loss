"""Focused tests for the publisher-protein compatibility publication gate."""

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

from geneloss_repro.publisher_protein_qc import (
    PublishedProteinCompatibilityError,
    audit_published_protein_compatibility,
)


def write_fasta(path: Path, records: list[tuple[str, str]], *, compressed: bool = False) -> None:
    text = "".join(f">{identifier} description\n{sequence}\n" for identifier, sequence in records)
    if compressed:
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PublishedProteinCompatibilityTests(unittest.TestCase):
    def test_cli_help_works_outside_repository_without_pythonpath(self) -> None:
        script = (
            Path(__file__).parents[2]
            / "scripts"
            / "qc"
            / "audit_published_protein_compatibility.py"
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
        self.assertIn("--publisher-proteins", completed.stdout)

    def test_pass_bundle_audits_only_permitted_normalizations_and_is_path_free(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            derived = root / "derived.faa"
            publisher = root / "publisher.data"
            output = root / "audit"
            write_fasta(
                derived,
                [
                    ("p1", "MK"),
                    ("p2", "MDT"),
                    ("p3", "MP"),
                    ("p4", "MKT*"),
                    ("p5", "MX"),
                ],
            )
            # Gzip is detected from content even though the suffix is not .gz.
            write_fasta(
                publisher,
                [
                    ("p1", "MK"),
                    ("p2", "MXX"),
                    ("p3", "MP*"),
                    ("p4", "MKX"),
                    ("p5", "MX"),
                ],
                compressed=True,
            )

            result = audit_published_protein_compatibility(
                derived, publisher, output, "act_test_hap1"
            )
            self.assertEqual(result.record_count, 5)
            self.assertEqual(result.exact_record_count, 2)
            self.assertEqual(result.normalized_exact_record_count, 3)
            self.assertEqual(result.terminal_stop_normalized_record_count, 2)
            self.assertEqual(result.publisher_x_wildcard_record_count, 2)
            self.assertEqual(result.publisher_x_wildcard_position_count, 3)

            rows = {
                row["protein_id"]: row
                for row in read_tsv(
                    output
                    / "act_test_hap1.published_protein_compatibility.records.tsv"
                )
            }
            self.assertEqual(rows["p1"]["status"], "PASS_EXACT")
            self.assertEqual(rows["p2"]["status"], "PASS_PUBLISHER_X_WILDCARD")
            self.assertEqual(rows["p2"]["publisher_X_wildcard_count"], "2")
            self.assertEqual(rows["p2"]["publisher_X_wildcard_positions_1based"], "2,3")
            self.assertEqual(rows["p3"]["status"], "PASS_TERMINAL_STOP_NORMALIZED")
            self.assertEqual(
                rows["p4"]["status"],
                "PASS_TERMINAL_STOP_AND_PUBLISHER_X_WILDCARD",
            )
            # Equal X/X is an exact match, not a use of the publisher wildcard.
            self.assertEqual(rows["p5"]["publisher_X_wildcard_count"], "0")

            summary = read_tsv(
                output / "act_test_hap1.published_protein_compatibility.summary.tsv"
            )[0]
            self.assertEqual(summary["status"], "PASS")
            self.assertEqual(summary["publication_gate"], "PASS")
            self.assertEqual(summary["exact_ID_set"], "true")
            self.assertEqual(summary["normalized_exact_record_count"], "3")
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["publication_gate"], "PASS")

            checksums = read_tsv(output / "checksums.tsv")
            self.assertEqual(
                {row["file"] for row in checksums},
                {
                    "act_test_hap1.published_protein_compatibility.records.tsv",
                    "act_test_hap1.published_protein_compatibility.summary.tsv",
                    "run_manifest.json",
                },
            )
            for row in checksums:
                payload = (output / row["file"]).read_bytes()
                self.assertEqual(len(payload), int(row["bytes"]))
                self.assertEqual(hashlib.sha256(payload).hexdigest(), row["sha256"])
            published_text = "".join(
                path.read_text(encoding="utf-8")
                for path in output.iterdir()
                if path.is_file()
            )
            self.assertNotIn(str(root), published_text)

    def test_exact_identifier_set_is_required_and_failure_is_atomic(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            derived = root / "derived.faa"
            publisher = root / "publisher.faa"
            output = root / "audit"
            write_fasta(derived, [("shared", "MK"), ("derived_only", "MA")])
            write_fasta(publisher, [("shared", "MK"), ("publisher_only", "MP")])
            with self.assertRaisesRegex(
                PublishedProteinCompatibilityError,
                "ID sets are not exactly equal.*derived_only=1.*publisher_only=1",
            ):
                audit_published_protein_compatibility(
                    derived, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())
            self.assertFalse(list(root.glob(".audit.staging.*")))

    def test_non_x_residue_difference_rejects_complete_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            derived = root / "derived.faa"
            publisher = root / "publisher.faa"
            output = root / "audit"
            write_fasta(derived, [("p1", "MKT")])
            write_fasta(publisher, [("p1", "MAT")])
            with self.assertRaisesRegex(
                PublishedProteinCompatibilityError, "FAIL_RESIDUE_MISMATCH"
            ):
                audit_published_protein_compatibility(
                    derived, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())

    def test_x_wildcard_is_one_way_and_requires_canonical_derived_residue(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            derived = root / "derived.faa"
            publisher = root / "publisher.faa"
            output = root / "audit"
            write_fasta(derived, [("p1", "MX"), ("p2", "MU")])
            write_fasta(publisher, [("p1", "MA"), ("p2", "MX")])
            with self.assertRaisesRegex(
                PublishedProteinCompatibilityError, "failed for 2 records"
            ):
                audit_published_protein_compatibility(
                    derived, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())

    def test_internal_or_repeated_terminal_stop_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            derived = root / "derived.faa"
            publisher = root / "publisher.faa"
            output = root / "audit"
            write_fasta(derived, [("p1", "MK")])
            write_fasta(publisher, [("p1", "MK**")])
            with self.assertRaisesRegex(
                PublishedProteinCompatibilityError, "internal or repeated stop"
            ):
                audit_published_protein_compatibility(
                    derived, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())

    def test_duplicate_identifier_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            derived = root / "derived.faa"
            publisher = root / "publisher.faa"
            output = root / "audit"
            derived.write_text(">p1\nMK\n>p1 duplicate\nMK\n", encoding="utf-8")
            write_fasta(publisher, [("p1", "MK")])
            with self.assertRaisesRegex(PublishedProteinCompatibilityError, "repeats ID"):
                audit_published_protein_compatibility(
                    derived, publisher, output, "act_test"
                )
            self.assertFalse(output.exists())

    def test_existing_output_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            derived = root / "derived.faa"
            publisher = root / "publisher.faa"
            output = root / "audit"
            write_fasta(derived, [("p1", "MK")])
            write_fasta(publisher, [("p1", "MK")])
            output.mkdir()
            sentinel = output / "keep.txt"
            sentinel.write_text("keep\n", encoding="utf-8")
            with self.assertRaisesRegex(
                PublishedProteinCompatibilityError, "already exists"
            ):
                audit_published_protein_compatibility(
                    derived, publisher, output, "act_test"
                )
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep\n")


if __name__ == "__main__":
    unittest.main()
