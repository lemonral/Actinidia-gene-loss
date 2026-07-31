"""Fail-closed tests for download destination handling."""

from __future__ import annotations

import importlib.util
from unittest import mock
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[2] / "scripts" / "download" / "fetch_manifest.py"
SPEC = importlib.util.spec_from_file_location("fetch_manifest", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class DownloadPathTest(unittest.TestCase):
    def test_relative_destination_stays_below_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory).resolve()
            expected = root / "downloads" / "Actinidia_eriantha" / "genome.fa.gz"
            self.assertEqual(
                MODULE.destination_path(root, "downloads/Actinidia_eriantha/genome.fa.gz"),
                expected,
            )

    def test_parent_escape_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory).resolve()
            with self.assertRaises(MODULE.ManifestError):
                MODULE.destination_path(root, "../escape.bin")

    def test_absolute_destination_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory).resolve()
            with self.assertRaises(MODULE.ManifestError):
                MODULE.destination_path(root, "/tmp/escape.bin")

    @mock.patch.object(MODULE.subprocess, "run")
    @mock.patch.object(MODULE.shutil, "which", return_value="/usr/bin/aria2c")
    def test_segmented_download_passes_proxy_to_aria2(self, _which, run) -> None:
        MODULE.aria2_download(
            "https://ftp.ncbi.nlm.nih.gov/example.fa.gz",
            Path("/tmp/example.fa.gz.part"),
            connections=4,
            proxy="http://127.0.0.1:7890",
            environment={"PATH": "/usr/bin"},
        )
        command = run.call_args.args[0]
        self.assertIn("--all-proxy", command)
        self.assertEqual(command[command.index("--all-proxy") + 1], "http://127.0.0.1:7890")
        self.assertEqual(command[-1], "https://ftp.ncbi.nlm.nih.gov/example.fa.gz")

    @mock.patch.object(MODULE.subprocess, "run")
    @mock.patch.object(MODULE.shutil, "which", return_value="/usr/bin/aria2c")
    def test_safe_aria2_resume_can_use_one_connection(self, _which, run) -> None:
        MODULE.aria2_download(
            "https://download.example.org/example.fa.gz",
            Path("/tmp/example.fa.gz.part"),
            connections=1,
            proxy=None,
            environment={"PATH": "/usr/bin"},
        )
        command = run.call_args.args[0]
        self.assertEqual(command[command.index("--max-connection-per-server") + 1], "1")
        self.assertEqual(command[command.index("--split") + 1], "1")
        self.assertNotIn("--all-proxy", command)

    @mock.patch.object(MODULE.subprocess, "run")
    def test_curl_retry_budget_and_resume_are_explicit(self, run) -> None:
        MODULE.curl_download(
            "https://download.example.org/example.fa.gz",
            Path("/tmp/example.fa.gz.part"),
            proxy=None,
            direct_domains=[],
            retries=50,
            environment={"PATH": "/usr/bin"},
        )
        command = run.call_args.args[0]
        self.assertEqual(command[command.index("--retry") + 1], "50")
        self.assertEqual(command[command.index("--continue-at") + 1], "-")
        self.assertIn("--retry-all-errors", command)


if __name__ == "__main__":
    unittest.main()
