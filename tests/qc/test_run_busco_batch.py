"""Focused tests for BUSCO gzip staging and restart provenance."""

from __future__ import annotations

import csv
import dataclasses
import gzip
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT_DIR = ROOT / "scripts" / "qc"
SCRIPT = SCRIPT_DIR / "run_busco_batch.py"
sys.path.insert(0, str(SCRIPT_DIR))
SPEC = importlib.util.spec_from_file_location("run_busco_batch_gzip_test", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


MANIFEST_COLUMNS = (
    "sample",
    "current_or_alternative",
    "accession",
    "genome",
    "gff",
    "protein",
    "source_url",
)


class BatchFixture:
    def __init__(self, root: Path, payload: bytes, *, compressed: bool, mode: str = "genome"):
        self.root = root
        self.mode = mode
        self.sample = "act_test_hap1"
        self.output = root / "output"
        self.lineage = root / "embryophyta_odb10"
        self.lineage.mkdir()
        self.busco = root / "fake_busco"
        self.busco.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        self.busco.chmod(0o755)
        # Deliberately omit a .gz suffix: detection must use magic bytes.
        self.input = root / ("input.fasta.payload" if compressed else "input.fasta")
        self.write_payload(payload, compressed=compressed)
        self.gff = root / "input.gff3"
        self.gff.write_text("##gff-version 3\n", encoding="utf-8")
        self.manifest = root / "manifest.tsv"
        row = {
            "sample": self.sample,
            "current_or_alternative": "candidate",
            "accession": "TEST_1",
            "genome": str(self.input),
            "gff": str(self.gff),
            "protein": str(self.input),
            "source_url": "https://example.org/test",
        }
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=MANIFEST_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(row)

    def write_payload(self, payload: bytes, *, compressed: bool = True) -> None:
        data = gzip.compress(payload, mtime=0) if compressed else payload
        self.input.write_bytes(data)

    def args(self, *extra: str):
        return MODULE.build_parser().parse_args(
            [
                "--manifest",
                str(self.manifest),
                "--output-dir",
                str(self.output),
                "--mode",
                self.mode,
                "--lineage",
                str(self.lineage),
                "--busco",
                str(self.busco),
                "--jobs",
                "1",
                "--cpus-per-job",
                "1",
                *extra,
            ]
        )

    def write_summary(self, input_path: Path) -> Path:
        run_dir = self.output / "runs" / self.sample
        run_dir.mkdir(parents=True, exist_ok=True)
        summary = run_dir / f"short_summary.specific.embryophyta_odb10.{self.sample}.txt"
        summary.write_text(
            "\n".join(
                (
                    "# BUSCO version is: 5.8.2",
                    "# The lineage dataset is: embryophyta_odb10 (Creation date: 2024-01-01)",
                    f"# BUSCO was run in mode: {self.mode}",
                    "# Summarized benchmarking in BUSCO notation for file " + str(input_path),
                    "C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:10",
                    "",
                )
            ),
            encoding="utf-8",
        )
        return summary


class BuscoGzipStagingTest(unittest.TestCase):
    def test_gzip_stage_is_content_addressed_checksum_validated_and_reused(self) -> None:
        payload = b">chr1\nACGTNN\n>chr2\nTTAA\n"
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = BatchFixture(Path(temporary_directory), payload, compressed=True)
            jobs, immediate, runnable = MODULE.build_jobs(fixture.args("--validate-only"))

            self.assertEqual(len(jobs), 1)
            self.assertEqual(len(immediate), 1)
            self.assertEqual(runnable, [])
            job = jobs[0]
            self.assertEqual(job.input_path, fixture.input.resolve())
            self.assertNotEqual(job.effective_input_path, job.input_path)
            self.assertEqual(job.effective_input_path.read_bytes(), payload)
            self.assertEqual(job.command[job.command.index("--in") + 1], str(job.effective_input_path))
            self.assertEqual(job.command.count("--opt-out-run-stats"), 1)
            self.assertIn("--opt-out-run-stats", immediate[0]["command"])
            self.assertEqual(immediate[0]["input_path"], str(fixture.input.resolve()))
            self.assertIn("gzip FASTA was staged", immediate[0]["message"])

            expected_sha = hashlib.sha256(payload).hexdigest()
            self.assertEqual(job.staged_input.staged_sha256, expected_sha)
            self.assertTrue(job.effective_input_path.name.endswith(f".{expected_sha}.fasta"))
            provenance = json.loads(job.staged_input.provenance_path.read_text(encoding="utf-8"))
            self.assertEqual(provenance["source_sha256"], hashlib.sha256(fixture.input.read_bytes()).hexdigest())
            self.assertEqual(provenance["staged_sha256"], expected_sha)

            staged_mtime = job.effective_input_path.stat().st_mtime_ns
            provenance_mtime = job.staged_input.provenance_path.stat().st_mtime_ns
            jobs_again, _, _ = MODULE.build_jobs(fixture.args("--validate-only"))
            self.assertEqual(jobs_again[0].effective_input_path, job.effective_input_path)
            self.assertEqual(job.effective_input_path.stat().st_mtime_ns, staged_mtime)
            self.assertEqual(job.staged_input.provenance_path.stat().st_mtime_ns, provenance_mtime)
            self.assertEqual(
                list(job.effective_input_path.parent.glob(".*.decompress.tmp.*")),
                [],
            )

    def test_truncated_gzip_fails_without_partial_stage_or_provenance(self) -> None:
        payload = b">chr1\nACGT\n"
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = BatchFixture(Path(temporary_directory), payload, compressed=True)
            fixture.input.write_bytes(fixture.input.read_bytes()[:-5])
            with self.assertRaisesRegex(MODULE.BatchInputError, "decompress gzip FASTA"):
                MODULE.build_jobs(fixture.args("--validate-only"))

            staging_dir = fixture.output / "staged_inputs" / fixture.mode
            self.assertFalse((staging_dir / f"{fixture.sample}.provenance.json").exists())
            self.assertEqual(list(staging_dir.glob(".*.decompress.tmp.*")), [])
            self.assertEqual(list(staging_dir.glob("*.fasta")), [])

    def test_corrupt_reused_stage_fails_closed_instead_of_being_rebuilt(self) -> None:
        payload = b">chr1\nACGT\n"
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = BatchFixture(Path(temporary_directory), payload, compressed=True)
            jobs, _, _ = MODULE.build_jobs(fixture.args("--validate-only"))
            staged = jobs[0].effective_input_path
            staged.write_bytes(b">chr1\nCORRUPTED\n")

            with self.assertRaisesRegex(MODULE.BatchInputError, "failed provenance checksum"):
                MODULE.build_jobs(fixture.args("--validate-only"))

    def test_malformed_stage_provenance_fails_closed(self) -> None:
        payload = b">chr1\nACGT\n"
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = BatchFixture(Path(temporary_directory), payload, compressed=True)
            jobs, _, _ = MODULE.build_jobs(fixture.args("--validate-only"))
            provenance_path = jobs[0].staged_input.provenance_path
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            provenance["unreviewed_override"] = True
            provenance_path.write_text(json.dumps(provenance), encoding="utf-8")

            with self.assertRaisesRegex(MODULE.BatchInputError, "unexpected field set"):
                MODULE.build_jobs(fixture.args("--validate-only"))

    def test_partial_gzip_run_restarts_only_with_exact_binding(self) -> None:
        payload = b">chr1\nACGT\n"
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = BatchFixture(Path(temporary_directory), payload, compressed=True)
            jobs, _, runnable = MODULE.build_jobs(fixture.args())
            self.assertEqual(len(runnable), 1)
            job = jobs[0]
            MODULE.atomic_write_json(job.input_binding_path, MODULE.run_input_binding_payload(job))
            job.run_dir.mkdir(parents=True)

            restarted_jobs, _, restarted = MODULE.build_jobs(fixture.args())
            self.assertEqual(len(restarted), 1)
            self.assertIn("--restart", restarted_jobs[0].command)

            fixture.write_payload(b">chr1\nTTTT\n")
            with self.assertRaisesRegex(MODULE.BatchInputError, "cannot be safely reused or restarted"):
                MODULE.build_jobs(fixture.args())

            forced_jobs, _, forced = MODULE.build_jobs(fixture.args("--force"))
            self.assertEqual(len(forced), 1)
            self.assertIn("--force", forced_jobs[0].command)
            self.assertNotIn("--restart", forced_jobs[0].command)

    def test_queued_job_revalidates_stage_before_starting_busco(self) -> None:
        payload = b">chr1\nACGT\n"
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = BatchFixture(Path(temporary_directory), payload, compressed=True)
            marker = fixture.root / "busco_started"
            fixture.busco.write_text(
                f"#!/bin/sh\ntouch {marker}\nexit 0\n",
                encoding="utf-8",
            )
            jobs, _, _ = MODULE.build_jobs(fixture.args())
            jobs[0].effective_input_path.write_bytes(b">chr1\nCHANGED\n")

            result = MODULE.run_job(jobs[0])
            self.assertEqual(result["status"], "failed_runner")
            self.assertIn("changed after job planning", result["message"])
            self.assertFalse(marker.exists())
            self.assertFalse(jobs[0].input_binding_path.exists())

    def test_execution_fails_closed_if_run_stats_opt_out_is_removed(self) -> None:
        payload = b">chr1\nACGT\n"
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = BatchFixture(Path(temporary_directory), payload, compressed=False)
            marker = fixture.root / "busco_started"
            fixture.busco.write_text(
                f"#!/bin/sh\ntouch {marker}\nexit 0\n",
                encoding="utf-8",
            )
            jobs, _, _ = MODULE.build_jobs(fixture.args())
            command_without_opt_out = tuple(
                token for token in jobs[0].command if token != "--opt-out-run-stats"
            )
            malformed_job = dataclasses.replace(jobs[0], command=command_without_opt_out)

            result = MODULE.run_job(malformed_job)

            self.assertEqual(result["status"], "failed_runner")
            self.assertIn("requires exactly one --opt-out-run-stats", result["message"])
            self.assertIn("observed 0", result["message"])
            self.assertFalse(marker.exists())
            self.assertFalse(malformed_job.stdout_log.exists())
            self.assertFalse(malformed_job.stderr_log.exists())

    def test_privacy_contract_rejects_duplicate_opt_out_flag(self) -> None:
        with self.assertRaisesRegex(
            MODULE.BatchInputError,
            "requires exactly one --opt-out-run-stats; observed 2",
        ):
            MODULE.validate_busco_privacy_command(
                ("busco", "--opt-out-run-stats", "--opt-out-run-stats")
            )

    def test_complete_gzip_run_requires_binding_and_exact_summary_input(self) -> None:
        payload = b">p1\nMPEPTIDE\n"
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = BatchFixture(
                Path(temporary_directory), payload, compressed=True, mode="proteins"
            )
            jobs, _, _ = MODULE.build_jobs(fixture.args())
            job = jobs[0]
            fixture.write_summary(job.effective_input_path)

            with self.assertRaisesRegex(MODULE.BatchInputError, "run-input binding is missing"):
                MODULE.build_jobs(fixture.args())

            MODULE.atomic_write_json(job.input_binding_path, MODULE.run_input_binding_payload(job))
            skipped_jobs, immediate, runnable = MODULE.build_jobs(fixture.args())
            self.assertEqual(runnable, [])
            self.assertEqual(immediate[0]["status"], "skipped_complete")
            self.assertEqual(skipped_jobs[0].input_path, fixture.input.resolve())

            fixture.write_summary(fixture.input.resolve())
            with self.assertRaisesRegex(MODULE.BatchInputError, "does not match staged input"):
                MODULE.build_jobs(fixture.args())

    def test_plain_fasta_keeps_original_path_restart_and_skip_behavior(self) -> None:
        payload = b">p1\nMPEPTIDE\n"
        with tempfile.TemporaryDirectory() as temporary_directory:
            fixture = BatchFixture(
                Path(temporary_directory), payload, compressed=False, mode="proteins"
            )
            run_dir = fixture.output / "runs" / fixture.sample
            run_dir.mkdir(parents=True)
            jobs, _, runnable = MODULE.build_jobs(fixture.args())
            self.assertEqual(len(runnable), 1)
            self.assertIsNone(jobs[0].staged_input)
            self.assertEqual(jobs[0].effective_input_path, fixture.input.resolve())
            self.assertEqual(jobs[0].command[jobs[0].command.index("--in") + 1], str(fixture.input.resolve()))
            self.assertIn("--restart", jobs[0].command)
            self.assertFalse((fixture.output / "staged_inputs").exists())

            fixture.write_summary(fixture.input.resolve())
            _, immediate, runnable = MODULE.build_jobs(fixture.args())
            self.assertEqual(runnable, [])
            self.assertEqual(immediate[0]["status"], "skipped_complete")


if __name__ == "__main__":
    unittest.main()
