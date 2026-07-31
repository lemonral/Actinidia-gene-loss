"""Tests for the fail-closed sequential MCMCTree runner."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "phylogeny" / "run_mcmctree_timetree_bundle.py"
SPEC = importlib.util.spec_from_file_location("run_mcmctree_timetree_bundle", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class RunMCMCTreeBundleTests(unittest.TestCase):
    def make_bundle(
        self, root: Path, executable: Path, name: str = "bundle", explicit_in_bv: bool = False
    ) -> Path:
        bundle = root / name
        bundle.mkdir()
        (bundle / "codon_positions.phy").write_text(" 2 1\nA A\nB A\n", encoding="ascii")
        (bundle / "calibrated_topology.trees").write_text(" 2 1\n\n(A,B);\n", encoding="ascii")
        controls = {
            "prior.ctl": (0, "prior.out", "prior.mcmc.txt"),
            "hessian.ctl": (3, "hessian.out", "hessian.mcmc.txt"),
            "posterior_chain1.ctl": (2, "posterior_chain1.out", "posterior_chain1.mcmc.txt"),
            "posterior_chain2.ctl": (2, "posterior_chain2.out", "posterior_chain2.mcmc.txt"),
        }
        for name, (usedata, outfile, mcmcfile) in controls.items():
            usedata_text = f"{usedata} in.BV" if usedata == 2 and explicit_in_bv else str(usedata)
            (bundle / name).write_text(
                f"usedata = {usedata_text}\noutfile = {outfile}\nmcmcfile = {mcmcfile}\n",
                encoding="ascii",
            )
        manifest = {
            "workflow": "mcmctree_timetree_secondary_bundle",
            "status": "PASS_PREPARED",
            "calibration_claim": "TimeTree secondary-calibrated; not fossil-calibrated",
            "active_constraint_count": 4,
            "mcmctree": {"version": "4.10.10", "sha256": digest(executable)},
        }
        (bundle / "run_manifest.json").write_text(json.dumps(manifest) + "\n", encoding="utf-8")
        rows = [
            f"{path.name}\t{digest(path)}\n"
            for path in sorted(bundle.iterdir())
            if path.name != "checksums.tsv"
        ]
        (bundle / "checksums.tsv").write_text("file\tsha256\n" + "".join(rows), encoding="utf-8")
        return bundle

    def make_executable(self, root: Path) -> Path:
        executable = root / "mcmctree"
        executable.write_text(
            """#!/usr/bin/env python3
import pathlib, re, sys
if len(sys.argv) == 1:
    print('MCMCTREE in paml version 4.10.10')
    raise SystemExit(255)
text = pathlib.Path(sys.argv[1]).read_text()
def value(name):
    return re.search(r'^\\s*' + name + r'\\s*=\\s*(\\S+)', text, re.M).group(1)
pathlib.Path(value('outfile')).write_text('output\\n')
mode = int(value('usedata'))
if mode == 3:
    pathlib.Path('out.BV').write_text('hessian\\n')
else:
    pathlib.Path(value('mcmcfile')).write_text('sample\\n')
""",
            encoding="utf-8",
        )
        executable.chmod(0o755)
        baseml = root / "baseml"
        baseml.write_text(
            "#!/usr/bin/env python3\nprint('BASEML in paml version 4.10.10')\n",
            encoding="utf-8",
        )
        baseml.chmod(0o755)
        return executable

    def test_runner_completes_all_stages_and_binds_hessian(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = self.make_executable(root)
            bundle = self.make_bundle(root, executable)
            output = root / "output"
            original = os.sys.argv
            try:
                os.sys.argv = [
                    str(SCRIPT), "--bundle", str(bundle), "--mcmctree", str(executable),
                    "--output-dir", str(output),
                ]
                self.assertEqual(MODULE.main(), 0)
            finally:
                os.sys.argv = original
            state = json.loads((output / "state.json").read_text())
            self.assertEqual(state["status"], "PASS_RUN_COMPLETE")
            self.assertEqual(len(state["completed"]), 4)
            self.assertEqual(state["baseml"]["version"], "4.10.10")
            hessian = digest(output / "runs" / "hessian" / "out.BV")
            self.assertEqual(state["hessian_sha256"], hessian)
            self.assertEqual(digest(output / "runs" / "posterior_chain1" / "in.BV"), hessian)
            self.assertEqual(digest(output / "runs" / "posterior_chain2" / "in.BV"), hessian)

    def test_checksum_tampering_is_rejected_before_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = self.make_executable(root)
            bundle = self.make_bundle(root, executable)
            (bundle / "prior.ctl").write_text("usedata = 9\n", encoding="ascii")
            output = root / "output"
            original = os.sys.argv
            try:
                os.sys.argv = [
                    str(SCRIPT), "--bundle", str(bundle), "--mcmctree", str(executable),
                    "--output-dir", str(output),
                ]
                self.assertEqual(MODULE.main(), 2)
            finally:
                os.sys.argv = original
            self.assertFalse(output.exists())

    def test_registered_missing_baseml_failure_resumes_without_prior_rerun(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = self.make_executable(root)
            bundle = self.make_bundle(root, executable)
            output = root / "output"
            original = os.sys.argv
            try:
                os.sys.argv = [
                    str(SCRIPT), "--bundle", str(bundle), "--mcmctree", str(executable),
                    "--output-dir", str(output),
                ]
                self.assertEqual(MODULE.main(), 0)
            finally:
                os.sys.argv = original
            prior = output / "runs" / "prior" / "prior.out"
            prior_sha256 = digest(prior)
            state_path = output / "state.json"
            state = json.loads(state_path.read_text())
            state["status"] = "ERROR"
            state["active_stage"] = "hessian"
            state["error"] = f"missing, empty, or symlink file: {output}/runs/hessian/out.BV"
            state["completed"] = state["completed"][:1]
            state.pop("baseml")
            state_path.write_text(json.dumps(state) + "\n", encoding="utf-8")
            shutil.rmtree(output / "runs" / "hessian")
            shutil.rmtree(output / "runs" / "posterior_chain1")
            shutil.rmtree(output / "runs" / "posterior_chain2")
            failed = output / "runs" / "hessian"
            failed.mkdir()
            (failed / "out.BV").write_bytes(b"")
            (failed / "console.stderr").write_text("sh: 1: baseml: not found\n", encoding="utf-8")
            try:
                os.sys.argv = [
                    str(SCRIPT), "--bundle", str(bundle), "--mcmctree", str(executable),
                    "--output-dir", str(output), "--resume-failed-missing-baseml",
                ]
                self.assertEqual(MODULE.main(), 0)
            finally:
                os.sys.argv = original
            resumed = json.loads(state_path.read_text())
            self.assertEqual(resumed["status"], "PASS_RUN_COMPLETE")
            self.assertEqual(len(resumed["completed"]), 4)
            self.assertEqual(digest(prior), prior_sha256)
            self.assertTrue((output / "runs" / "hessian.initial_missing_baseml").is_dir())
            self.assertTrue(resumed["recoveries"][0]["prior_reused_without_rerun"])

    def test_corrected_bundle_reuses_exact_prior_and_hessian(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = self.make_executable(root)
            old_bundle = self.make_bundle(root, executable, "bundle_v5")
            source = root / "source"
            original = os.sys.argv
            try:
                os.sys.argv = [
                    str(SCRIPT), "--bundle", str(old_bundle), "--mcmctree", str(executable),
                    "--output-dir", str(source),
                ]
                self.assertEqual(MODULE.main(), 0)
            finally:
                os.sys.argv = original
            source_state_path = source / "state.json"
            source_state = json.loads(source_state_path.read_text())
            source_state["status"] = "ERROR"
            source_state["active_stage"] = "posterior_chain1"
            source_state["error"] = "posterior_chain1: MCMCTree exited 1"
            source_state["completed"] = source_state["completed"][:2]
            source_state_path.write_text(json.dumps(source_state) + "\n", encoding="utf-8")
            (source / "runs" / "posterior_chain1" / "posterior_chain1.mcmc.txt").unlink()
            (source / "runs" / "posterior_chain1" / "console.stderr").write_text(
                "error: file name empty.\n", encoding="utf-8"
            )
            shutil.rmtree(source / "runs" / "posterior_chain2")
            corrected = self.make_bundle(root, executable, "bundle_v6", explicit_in_bv=True)
            output = root / "output"
            prior_sha256 = digest(source / "runs" / "prior" / "prior.out")
            hessian_sha256 = digest(source / "runs" / "hessian" / "out.BV")
            try:
                os.sys.argv = [
                    str(SCRIPT), "--bundle", str(corrected), "--mcmctree", str(executable),
                    "--output-dir", str(output), "--reuse-prior-hessian-from", str(source),
                ]
                self.assertEqual(MODULE.main(), 0)
            finally:
                os.sys.argv = original
            state = json.loads((output / "state.json").read_text())
            self.assertEqual(state["status"], "PASS_RUN_COMPLETE")
            self.assertEqual(len(state["completed"]), 4)
            self.assertEqual(digest(output / "runs" / "prior" / "prior.out"), prior_sha256)
            self.assertEqual(digest(output / "runs" / "hessian" / "out.BV"), hessian_sha256)
            self.assertEqual(state["reused_stages"]["stages"], ["prior", "hessian"])


if __name__ == "__main__":
    unittest.main()
