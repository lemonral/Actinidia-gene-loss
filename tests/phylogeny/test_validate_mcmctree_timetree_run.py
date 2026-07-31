"""Tests for two-chain MCMCTree validation and pooled ultrametric publication."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
from Bio import Phylo


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "phylogeny" / "validate_mcmctree_timetree_run.py"
SPEC = importlib.util.spec_from_file_location("validate_mcmctree_timetree_run", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ValidateMCMCTreeRunTests(unittest.TestCase):
    def make_fixture(self, root: Path, *, divergent: bool = False) -> tuple[Path, Path]:
        executable = root / "mcmctree"
        baseml = root / "baseml"
        executable.write_text("mcmctree\n", encoding="ascii")
        baseml.write_text("baseml\n", encoding="ascii")
        bundle = root / "bundle"
        bundle.mkdir()
        (bundle / "codon_positions.phy").write_text(" 3 1\nA A\nB A\nC A\n", encoding="ascii")
        (bundle / "calibrated_topology.trees").write_text(
            " 3 1\n\n((A,B) 'B(0.3,0.7)',C) 'B(0.9,1.3)';\n", encoding="ascii"
        )
        controls = {
            "prior.ctl": (0, "prior.out", "prior.mcmc.txt", 100),
            "hessian.ctl": (3, "hessian.out", "hessian.mcmc.txt", 100),
            "posterior_chain1.ctl": (2, "posterior_chain1.out", "posterior_chain1.mcmc.txt", 1000),
            "posterior_chain2.ctl": (2, "posterior_chain2.out", "posterior_chain2.mcmc.txt", 1000),
        }
        for name, (mode, outfile, mcmcfile, nsample) in controls.items():
            (bundle / name).write_text(
                f"usedata = {mode}\noutfile = {outfile}\nmcmcfile = {mcmcfile}\n"
                f"sampfreq = 10\nnsample = {nsample}\n",
                encoding="ascii",
            )
        manifest = {
            "schema_version": 1,
            "workflow": "mcmctree_timetree_secondary_bundle",
            "status": "PASS_PREPARED",
            "calibration_claim": MODULE.CALIBRATION_CLAIM,
            "tip_count": 3,
            "tip_order": ["A", "B", "C"],
            "time_unit_ma": 100.0,
            "active_constraints": [
                {
                    "constraint_id": "root",
                    "node_label": "root",
                    "descendant_a": "A",
                    "descendant_b": "C",
                    "minimum_ma": 90.0,
                    "maximum_ma": 130.0,
                },
                {
                    "constraint_id": "ab",
                    "node_label": "ab",
                    "descendant_a": "A",
                    "descendant_b": "B",
                    "minimum_ma": 30.0,
                    "maximum_ma": 70.0,
                },
            ],
        }
        (bundle / "run_manifest.json").write_text(json.dumps(manifest) + "\n", encoding="utf-8")
        checksum_rows = [
            (path.name, digest(path)) for path in sorted(bundle.iterdir()) if path.is_file()
        ]
        (bundle / "checksums.tsv").write_text(
            "file\tsha256\n" + "".join(f"{name}\t{sha}\n" for name, sha in checksum_rows),
            encoding="utf-8",
        )

        run = root / "run"
        (run / "runs").mkdir(parents=True)
        completed = []
        rng = np.random.default_rng(20260720)
        numbered = "((1_A,2_B) 5,3_C) 4;"
        for stage, control_name in MODULE.STAGES:
            stage_dir = run / "runs" / stage
            stage_dir.mkdir()
            for name in (*MODULE.SHARED_INPUTS, control_name):
                (stage_dir / name).write_bytes((bundle / name).read_bytes())
            (stage_dir / "console.stdout").write_text("ok\n", encoding="ascii")
            (stage_dir / "console.stderr").write_bytes(b"")
            outputs = {}
            if stage == "hessian":
                (stage_dir / "hessian.out").write_text("hessian\n", encoding="ascii")
                (stage_dir / "out.BV").write_text("frozen\n", encoding="ascii")
                output_names = ["hessian.out", "out.BV"]
            elif stage.startswith("posterior"):
                (stage_dir / "in.BV").write_text("frozen\n", encoding="ascii")
                count = 1001
                root_age = rng.normal(1.1 + (0.12 if divergent and stage.endswith("2") else 0), 0.02, count)
                child_age = rng.normal(0.5, 0.01, count)
                mu = rng.normal(0.1, 0.01, count)
                likelihood = rng.normal(-50, 1, count)
                generations = np.concatenate(([1], np.arange(1, 1001) * 10))
                table = stage_dir / f"{stage}.mcmc.txt"
                with table.open("w", encoding="ascii") as handle:
                    handle.write("Gen\tt_n4\tt_n5\tmu1\tlnL\n")
                    for values in zip(generations, root_age, child_age, mu, likelihood, strict=True):
                        handle.write("\t".join(str(value) for value in values) + "\n")
                report = stage_dir / f"{stage}.out"
                report.write_text(
                    "Species tree for FigTree.  Branch lengths = posterior mean times; 95% CIs = labels\n"
                    + numbered + "\n\n"
                    + f"t_n4          {root_age.mean():.4f} ( 1.0, 1.2)\n"
                    + f"t_n5          {child_age.mean():.4f} ( 0.4, 0.6)\n",
                    encoding="ascii",
                )
                output_names = [report.name, table.name]
            else:
                (stage_dir / "prior.out").write_text("prior\n", encoding="ascii")
                (stage_dir / "prior.mcmc.txt").write_text("prior samples\n", encoding="ascii")
                output_names = ["prior.out", "prior.mcmc.txt"]
            for name in output_names:
                path = stage_dir / name
                outputs[name] = {"bytes": path.stat().st_size, "sha256": digest(path)}
            completed.append(
                {
                    "stage": stage,
                    "returncode": 0,
                    "control_sha256": digest(bundle / control_name),
                    "stdout_sha256": digest(stage_dir / "console.stdout"),
                    "stderr_sha256": digest(stage_dir / "console.stderr"),
                    "outputs": outputs,
                }
            )
        state = {
            "schema_version": 1,
            "workflow": "sequential_mcmctree_timetree_secondary",
            "status": "PASS_RUN_COMPLETE",
            "calibration_claim": MODULE.CALIBRATION_CLAIM,
            "prepared_bundle": {
                "path": str(bundle.resolve()),
                "manifest_sha256": digest(bundle / "run_manifest.json"),
                "checksums_sha256": digest(bundle / "checksums.tsv"),
            },
            "mcmctree": {"path": str(executable), "sha256": digest(executable), "version": "4.10.10"},
            "baseml": {"path": str(baseml), "sha256": digest(baseml), "version": "4.10.10"},
            "hessian_sha256": digest(run / "runs" / "hessian" / "out.BV"),
            "completed": completed,
        }
        (run / "state.json").write_text(json.dumps(state) + "\n", encoding="utf-8")
        return bundle, run

    def invoke(self, bundle: Path, run: Path, output: Path) -> int:
        original = os.sys.argv
        try:
            os.sys.argv = [
                str(SCRIPT), "--bundle", str(bundle), "--run-dir", str(run),
                "--output-dir", str(output), "--minimum-effective-sample-size", "100",
                "--maximum-split-rhat", "1.01",
            ]
            return MODULE.main()
        finally:
            os.sys.argv = original

    def test_valid_chains_publish_exact_ultrametric_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle, run = self.make_fixture(root)
            output = root / "validated"
            self.assertEqual(self.invoke(bundle, run, output), 0)
            validation = json.loads((output / "validation.json").read_text())
            self.assertEqual(validation["status"], "PASS_MCMCTREE_VALIDATED_ULTRAMETRIC")
            self.assertEqual(validation["chain_count"], 2)
            tree = Phylo.read(output / "dated_tree.mean_ma.tre", "newick")
            distances = [tree.distance(tree.root, tip) for tip in tree.get_terminals()]
            self.assertLess(max(distances) - min(distances), 1e-6)
            self.assertEqual({tip.name for tip in tree.get_terminals()}, {"A", "B", "C"})

    def test_chain_tampering_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle, run = self.make_fixture(root)
            chain = run / "runs" / "posterior_chain1" / "posterior_chain1.mcmc.txt"
            chain.write_text(chain.read_text() + "tamper\n", encoding="ascii")
            output = root / "validated"
            self.assertEqual(self.invoke(bundle, run, output), 2)
            self.assertFalse(output.exists())

    def test_divergent_chains_fail_split_rhat(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle, run = self.make_fixture(root, divergent=True)
            output = root / "validated"
            self.assertEqual(self.invoke(bundle, run, output), 2)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
