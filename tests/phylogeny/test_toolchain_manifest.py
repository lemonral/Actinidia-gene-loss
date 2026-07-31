from __future__ import annotations

import csv
import json
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "config" / "phylogeny" / "toolchain.tsv"


class ToolchainManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with MANIFEST.open(newline="", encoding="utf-8") as handle:
            cls.rows = list(csv.DictReader(handle, delimiter="\t"))
        cls.by_id = {row["tool_id"]: row for row in cls.rows}

    def test_exact_production_tool_set_and_versions(self) -> None:
        expected = {
            "orthofinder": "3.1.5",
            "iqtree2": "2.4.0",
            "astral_pro3": "1.25.3.8",
            "paml_mcmctree": "4.10.10",
            "cafe5": "5.1.0",
            "mafft": "7.526",
            "diamond": "2.1.10",
            "blast_plus": "2.16.0+",
            "r": "4.4.0",
        }
        self.assertEqual(set(self.by_id), set(expected))
        self.assertEqual(len(self.rows), len(self.by_id), "tool_id values must be unique")
        self.assertEqual(
            {tool_id: row["version"] for tool_id, row in self.by_id.items()},
            expected,
        )

    def test_release_and_license_links_are_official_https_urls(self) -> None:
        allowed_hosts = {
            "github.com",
            "mafft.cbrc.jp",
            "ftp.ncbi.nlm.nih.gov",
            "blast.ncbi.nlm.nih.gov",
            "cran.r-project.org",
            "www.r-project.org",
        }
        for row in self.rows:
            for field in ("artifact_url", "release_url", "license_url"):
                value = row[field]
                self.assertTrue(value.startswith("https://"), f"{row['tool_id']} {field}")
                host = value.split("/", 3)[2]
                self.assertIn(host, allowed_hosts, f"{row['tool_id']} {field}")
            self.assertTrue(row["license_expression"])

    def test_every_row_has_a_fail_closed_identity_lock(self) -> None:
        sha256 = re.compile(r"^[0-9a-f]{64}$")
        md5 = re.compile(r"^[0-9a-f]{32}$")
        for row in self.rows:
            artifact_sha = row["artifact_sha256"]
            commit = row["source_commit"]
            vendor_md5 = row["vendor_md5"]
            has_archive_sha = bool(sha256.fullmatch(artifact_sha))
            has_git_commit = bool(re.fullmatch(r"[0-9a-f]{40}", commit))
            has_vendor_md5 = bool(md5.fullmatch(vendor_md5))
            self.assertTrue(
                has_archive_sha or has_git_commit or has_vendor_md5,
                f"{row['tool_id']} lacks a checksum or immutable commit",
            )
            if artifact_sha.startswith("PENDING_"):
                self.assertEqual(row["tool_id"], "blast_plus")
                self.assertTrue(has_vendor_md5)
                self.assertIn("fail closed", row["production_action"].lower())

    def test_probes_distinguish_exact_version_from_identity_smoke(self) -> None:
        representative_probe_output = {
            "orthofinder": "OrthoFinder version 3.1.5",
            "iqtree2": "IQ-TREE multicore version 2.4.0",
            "astral_pro3": "Version: v1.25.3.8",
            "paml_mcmctree": "MCMCTREE in paml version 4.10.10, 27 Jan 2026",
            "cafe5": "CAFE: Computational Analysis of gene Family Evolution\nusage: cafe5",
            "mafft": "v7.526 (2024/Apr/26)",
            "diamond": "diamond version 2.1.10",
            "blast_plus": "blastp: 2.16.0+",
            "r": "R version 4.4.0 (2024-04-24)",
        }
        smoke_only = {"cafe5"}
        for row in self.rows:
            parsed_args = json.loads(row["probe_args_json"])
            self.assertIsInstance(parsed_args, list)
            self.assertTrue(all(isinstance(argument, str) for argument in parsed_args))
            allowed_exit_codes = json.loads(row["allowed_probe_exit_codes_json"])
            self.assertTrue(allowed_exit_codes)
            self.assertTrue(all(type(code) is int for code in allowed_exit_codes))
            self.assertEqual(allowed_exit_codes, sorted(set(allowed_exit_codes)))
            self.assertTrue(all(0 <= code <= 255 for code in allowed_exit_codes))
            self.assertEqual(
                row["allowed_probe_exit_codes_json"],
                json.dumps(allowed_exit_codes, separators=(",", ":")),
            )
            pattern = re.compile(row["probe_regex"])
            output = representative_probe_output[row["tool_id"]]
            self.assertIsNotNone(
                pattern.search(output),
                row["tool_id"],
            )
            expected_level = (
                "program_identity_smoke_only"
                if row["tool_id"] in smoke_only
                else "exact_version_banner"
            )
            self.assertEqual(row["probe_evidence_level"], expected_level)
            if row["tool_id"] in smoke_only:
                self.assertNotIn(row["version"], output)
            else:
                self.assertIn(row["version"], output)
            self.assertTrue(re.fullmatch(r"[A-Za-z0-9_.+-]+", row["executable"]))

        self.assertEqual(
            {
                row["tool_id"]: json.loads(row["allowed_probe_exit_codes_json"])
                for row in self.rows
            },
            {
                tool_id: ([255] if tool_id == "paml_mcmctree" else [0])
                for tool_id in self.by_id
            },
        )

    def test_server_status_is_explicit_and_not_claimed_ready(self) -> None:
        allowed = {
            "exact_version_observed_checksum_pending",
            "different_version_observed",
            "present_version_unverified",
            "program_identity_observed_version_provenance_pending",
            "not_observed",
        }
        for row in self.rows:
            self.assertIn(row["server_inventory_status"], allowed)
            self.assertNotIn("ready", row["server_inventory_status"])
            self.assertTrue(row["server_inventory_note"])
            self.assertTrue(row["production_action"])

        exact_but_unclosed = {
            "orthofinder",
            "mafft",
            "diamond",
            "blast_plus",
            "paml_mcmctree",
            "r",
        }
        self.assertEqual(
            {
                row["tool_id"]
                for row in self.rows
                if row["server_inventory_status"]
                == "exact_version_observed_checksum_pending"
            },
            exact_but_unclosed,
        )
        self.assertEqual(
            self.by_id["cafe5"]["server_inventory_status"],
            "program_identity_observed_version_provenance_pending",
        )
        self.assertIn("Do not infer 5.1.0", self.by_id["cafe5"]["production_action"])
        for tool_id in ("iqtree2", "astral_pro3"):
            self.assertEqual(self.by_id[tool_id]["server_inventory_status"], "not_observed")

    def test_checksum_provenance_does_not_overstate_local_hashes(self) -> None:
        upstream_published = {
            row["tool_id"]
            for row in self.rows
            if "upstream-published" in row["checksum_provenance"]
        }
        self.assertEqual(upstream_published, {"orthofinder", "paml_mcmctree"})
        for row in self.rows:
            provenance = row["checksum_provenance"]
            if "locally streamed" in provenance:
                self.assertNotIn("upstream-published", provenance)

    def test_no_private_paths_or_proxy_details_are_embedded(self) -> None:
        text = MANIFEST.read_text(encoding="utf-8")
        forbidden = ("/home/", "/Users/", "mihomo", "proxy", "subscription")
        for token in forbidden:
            self.assertNotIn(token, text)


if __name__ == "__main__":
    unittest.main()
