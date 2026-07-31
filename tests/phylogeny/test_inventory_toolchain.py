from __future__ import annotations

import csv
import hashlib
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "phylogeny" / "inventory_toolchain.py"
TOOLCHAIN = ROOT / "config" / "phylogeny" / "toolchain.tsv"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ToolchainInventoryTests(unittest.TestCase):
    @staticmethod
    def write_executable(path: Path, version: str) -> None:
        path.write_text(
            "#!/usr/bin/env python3\n"
            "import os\n"
            "if 'HTTPS_PROXY' in os.environ or 'PRIVATE_API_TOKEN' in os.environ:\n"
            "    raise SystemExit(7)\n"
            "from pathlib import Path\n"
            f"print('MAFFT {version}; private_path=' + str(Path(__file__).resolve()))\n",
            encoding="utf-8",
        )
        path.chmod(0o755)

    @staticmethod
    def write_registry(
        path: Path,
        executable: Path,
        *,
        tool_id: str = "mafft",
        executable_id: str | None = None,
        duplicate: bool = False,
    ) -> None:
        executable_id = executable_id or tool_id
        row = f"{tool_id}\t{executable_id}\t{executable}\n"
        path.write_text(
            "tool_id\texecutable_id\texecutable_path\n"
            + row
            + (row if duplicate else ""),
            encoding="utf-8",
        )

    @staticmethod
    def write_raw_probe_executable(
        path: Path, raw_output: bytes, *, exit_code: int = 0
    ) -> None:
        path.write_text(
            "#!/usr/bin/env python3\n"
            "import os\n"
            "import sys\n"
            "if os.environ.get('TERM') != 'dumb':\n"
            "    raise SystemExit(9)\n"
            f"sys.stdout.buffer.write(bytes.fromhex('{raw_output.hex()}'))\n"
            "sys.stdout.buffer.flush()\n"
            f"raise SystemExit({exit_code})\n",
            encoding="utf-8",
        )
        path.chmod(0o755)

    @staticmethod
    def run_inventory(
        *arguments: str,
        environment: dict[str, str] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SCRIPT), *arguments],
            text=True,
            capture_output=True,
            check=False,
            env=environment,
        )

    def test_exact_version_emits_path_free_hash_only_audit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "private-server-mafft"
            registry = root / "private-registry.tsv"
            output = root / "safe-audit.tsv"
            self.write_executable(executable, "v7.526")
            self.write_registry(registry, executable)

            environment = os.environ.copy()
            environment["HTTPS_PROXY"] = "http://credential-bearing.invalid"
            environment["PRIVATE_API_TOKEN"] = "must-not-reach-probe"
            completed = self.run_inventory(
                "--toolchain",
                str(TOOLCHAIN),
                "--registry",
                str(registry),
                "--output",
                str(output),
                "--strict",
                environment=environment,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["tool_id"], "mafft")
            self.assertEqual(row["executable_id"], "mafft")
            self.assertEqual(row["manifest_version"], "7.526")
            self.assertEqual(row["executable_sha256"], sha256(executable))
            self.assertEqual(row["toolchain_manifest_sha256"], sha256(TOOLCHAIN))
            self.assertEqual(row["registry_sha256"], sha256(registry))
            self.assertEqual(row["probe_evidence_level"], "exact_version_banner")
            self.assertEqual(row["post_probe_executable_sha256"], sha256(executable))
            self.assertEqual(row["executable_identity_stable"], "true")
            self.assertEqual(row["matched_probe_token"], "7.526")
            self.assertEqual(row["probe_match"], "true")
            self.assertEqual(row["version_match"], "true")
            self.assertEqual(row["inventory_status"], "PASS_EXACT_VERSION_MATCH")

            emitted = output.read_text(encoding="utf-8")
            self.assertNotIn(str(root), emitted)
            self.assertNotIn("private-server-mafft", emitted)
            self.assertNotIn("private_path", emitted)
            self.assertNotIn("must-not-reach-probe", emitted)

    def test_strict_mode_returns_one_for_version_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "mafft"
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            self.write_executable(executable, "v7.525")
            self.write_registry(registry, executable)

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
                "--strict",
            )
            self.assertEqual(completed.returncode, 1, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["version_match"], "false")
            self.assertEqual(row["matched_probe_token"], "NA")
            self.assertEqual(row["probe_match"], "false")
            self.assertEqual(row["inventory_status"], "FAIL_VERSION_MISMATCH")

    def test_exact_version_matches_strictly_sgr_stripped_view(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "orthofinder"
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            raw_output = (
                b"OrthoFinder version \x1b[1;32m3.1.5\x1b[0m\r\n"
            )
            normalized_output = b"OrthoFinder version 3.1.5\n"
            self.write_raw_probe_executable(executable, raw_output)
            self.write_registry(registry, executable, tool_id="orthofinder")

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
                "--strict",
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["probe_output_bytes"], str(len(raw_output)))
            self.assertEqual(
                row["probe_output_sha256"], hashlib.sha256(raw_output).hexdigest()
            )
            self.assertEqual(
                row["normalized_probe_output_sha256"],
                hashlib.sha256(normalized_output).hexdigest(),
            )
            self.assertEqual(row["stripped_sgr_sequence_count"], "2")
            self.assertEqual(row["matched_probe_token"], "3.1.5")
            self.assertEqual(row["version_match"], "true")
            self.assertEqual(row["inventory_status"], "PASS_EXACT_VERSION_MATCH")

    def test_non_sgr_escape_sequence_is_rejected_without_weakening_regex(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "orthofinder"
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            raw_output = b"OrthoFinder version 3.1.5\x1b[2J\n"
            self.write_raw_probe_executable(executable, raw_output)
            self.write_registry(registry, executable, tool_id="orthofinder")

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
                "--strict",
            )
            self.assertEqual(completed.returncode, 1, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["probe_output_bytes"], str(len(raw_output)))
            self.assertEqual(
                row["probe_output_sha256"], hashlib.sha256(raw_output).hexdigest()
            )
            self.assertEqual(row["normalized_probe_output_sha256"], "NA")
            self.assertEqual(row["probe_match"], "false")
            self.assertEqual(row["version_match"], "false")
            self.assertEqual(
                row["inventory_status"], "FAIL_PROBE_OUTPUT_CONTROL_SEQUENCE"
            )

    def test_declared_mcmctree_exit_255_passes_only_with_exact_version(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "mcmctree"
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            raw_output = b"MCMCTREE in paml version 4.10.10, 27 Jan 2026\n"
            self.write_raw_probe_executable(executable, raw_output, exit_code=255)
            self.write_registry(registry, executable, tool_id="paml_mcmctree")

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
                "--strict",
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["allowed_probe_exit_codes_json"], "[255]")
            self.assertEqual(row["probe_exit_code"], "255")
            self.assertEqual(row["matched_probe_token"], "4.10.10")
            self.assertEqual(row["version_match"], "true")
            self.assertEqual(row["inventory_status"], "PASS_EXACT_VERSION_MATCH")

    def test_declared_exit_does_not_allow_wrong_mcmctree_version(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "mcmctree"
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            raw_output = b"MCMCTREE in paml version 4.10.9, 1 Jan 2025\n"
            self.write_raw_probe_executable(executable, raw_output, exit_code=255)
            self.write_registry(registry, executable, tool_id="paml_mcmctree")

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
                "--strict",
            )
            self.assertEqual(completed.returncode, 1, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["probe_exit_code"], "255")
            self.assertEqual(row["version_match"], "false")
            self.assertEqual(row["inventory_status"], "FAIL_VERSION_MISMATCH")

    def test_undeclared_nonzero_exit_fails_despite_exact_version(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "mafft"
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            raw_output = b"MAFFT v7.526\n"
            self.write_raw_probe_executable(executable, raw_output, exit_code=3)
            self.write_registry(registry, executable)

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
                "--strict",
            )
            self.assertEqual(completed.returncode, 1, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["allowed_probe_exit_codes_json"], "[0]")
            self.assertEqual(row["probe_exit_code"], "3")
            self.assertEqual(row["matched_probe_token"], "7.526")
            self.assertEqual(row["version_match"], "false")
            self.assertEqual(row["inventory_status"], "FAIL_PROBE_EXIT")

    def test_cafe_help_is_identity_smoke_not_invented_version_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "cafe5"
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            executable.write_text(
                "#!/usr/bin/env python3\n"
                "print('CAFE: Computational Analysis of gene Family Evolution')\n"
                "print('usage: cafe5 [options]')\n",
                encoding="utf-8",
            )
            executable.chmod(0o755)
            self.write_registry(registry, executable, tool_id="cafe5")

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
                "--strict",
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["manifest_version"], "5.1.0")
            self.assertEqual(
                row["probe_evidence_level"], "program_identity_smoke_only"
            )
            self.assertEqual(row["matched_probe_token"], "CAFE")
            self.assertEqual(row["probe_match"], "true")
            self.assertEqual(row["version_match"], "not_tested")
            self.assertEqual(row["inventory_status"], "PASS_PROGRAM_IDENTITY_SMOKE")

    def test_executable_mutation_during_probe_fails_identity_check(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "mafft"
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            executable.write_text(
                "#!/usr/bin/env python3\n"
                "from pathlib import Path\n"
                "path = Path(__file__)\n"
                "print('MAFFT v7.526')\n"
                "path.write_text(path.read_text() + '# changed\\n')\n",
                encoding="utf-8",
            )
            executable.chmod(0o755)
            self.write_registry(registry, executable)

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
                "--strict",
            )
            self.assertEqual(completed.returncode, 1, completed.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertNotEqual(
                row["executable_sha256"], row["post_probe_executable_sha256"]
            )
            self.assertEqual(row["executable_identity_stable"], "false")
            self.assertEqual(row["inventory_status"], "FAIL_EXECUTABLE_CHANGED")

    def test_duplicate_registry_key_is_rejected_before_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "mafft"
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            self.write_executable(executable, "v7.526")
            self.write_registry(registry, executable, duplicate=True)

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("duplicates tool_id and executable_id", completed.stderr)
            self.assertFalse(output.exists())

    def test_noncanonical_allowed_exit_code_declaration_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "mafft"
            registry = root / "registry.tsv"
            toolchain = root / "toolchain.tsv"
            output = root / "audit.tsv"
            self.write_executable(executable, "v7.526")
            self.write_registry(registry, executable)
            with TOOLCHAIN.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
                fieldnames = list(reader.fieldnames or ())
            for row in rows:
                if row["tool_id"] == "mafft":
                    row["allowed_probe_exit_codes_json"] = "[0,0]"
            with toolchain.open("x", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fieldnames,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)

            completed = self.run_inventory(
                "--toolchain",
                str(toolchain),
                "--registry",
                str(registry),
                "--output",
                str(output),
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("canonical, nonempty, sorted, unique", completed.stderr)
            self.assertFalse(output.exists())

    def test_relative_executable_path_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            registry.write_text(
                "tool_id\texecutable_id\texecutable_path\n"
                "mafft\tmafft\trelative/bin/mafft\n",
                encoding="utf-8",
            )

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("must be absolute", completed.stderr)
            self.assertFalse(output.exists())

    def test_existing_output_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            executable = root / "mafft"
            registry = root / "registry.tsv"
            output = root / "audit.tsv"
            self.write_executable(executable, "v7.526")
            self.write_registry(registry, executable)
            output.write_text("keep me\n", encoding="utf-8")

            completed = self.run_inventory(
                "--registry",
                str(registry),
                "--output",
                str(output),
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("refusing to overwrite", completed.stderr)
            self.assertEqual(output.read_text(encoding="utf-8"), "keep me\n")


if __name__ == "__main__":
    unittest.main()
