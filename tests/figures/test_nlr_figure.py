"""Tests for the cohort-separated publication NLR figure."""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
import os
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
except ImportError:
    matplotlib = None


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPOSITORY_ROOT / "scripts" / "figures" / "make_nlr_figure.py"


HEADER = (
    "analysis_cohort\tcohort_role\tassembly_unit_id\tbiological_species\t"
    "haplotype_or_subgenome\tassembly_scope\ttotal_nlr_count\t"
    "positive_reference_nlr_loss_count\tcallable_reference_nlr_denominator\t"
    "positive_reference_nlr_loss_percentage\tpercentage_status\n"
)


class NlrFigureTest(unittest.TestCase):
    def run_figure(self, summary: Path, output: Path):
        environment = os.environ.copy()
        source_root = str(REPOSITORY_ROOT / "src")
        environment["PYTHONPATH"] = (
            source_root
            if not environment.get("PYTHONPATH")
            else source_root + os.pathsep + environment["PYTHONPATH"]
        )
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--unit-summary",
                str(summary),
                "--output-dir",
                str(output),
                "--basename",
                "R6_nlr_primary",
                "--dpi",
                "90",
            ],
            cwd=REPOSITORY_ROOT,
            env=environment,
            check=False,
            capture_output=True,
            text=True,
        )

    @unittest.skipUnless(matplotlib is not None, "optional matplotlib is not installed")
    def test_bundle_contains_counts_percentage_and_mathtext_taxa(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            summary = root / "nlr_unit_summary.tsv"
            summary.write_text(
                HEADER
                + "primary20\tprimary\tdeliciosa_a\tActinidia deliciosa\tA\tchromosome_partition\t100\t2\t10\t20.0\tdefined\n"
                + "primary20\tprimary\teriantha_hap1\tActinidia eriantha\tHAP1\t29_anchored_pseudochromosomes_only\t80\t0\t0\t\tundefined_zero_denominator\n",
                encoding="utf-8",
            )
            output = root / "R6_primary"
            completed = self.run_figure(summary, output)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue((output / "R6_nlr_primary.png").read_bytes().startswith(b"\x89PNG"))
            self.assertTrue((output / "R6_nlr_primary.pdf").read_bytes().startswith(b"%PDF"))
            plot_data = (output / "R6_nlr_primary.plot_data.tsv").read_text(encoding="utf-8")
            self.assertIn(r"$\mathit{Actinidia\ deliciosa}$ $\mathrm{A}$", plot_data)
            self.assertIn(r"$\mathit{Actinidia\ eriantha}$ $\mathrm{HAP1}$", plot_data)
            self.assertIn("positive_reference_nlr_loss_count", plot_data)
            self.assertIn("positive_reference_nlr_loss_percentage", plot_data)
            caption = (output / "R6_nlr_primary.caption.txt").read_text(encoding="utf-8")
            self.assertIn("primary and A. rufa sensitivity denominators are never pooled", caption)

    def test_mixed_primary_and_sensitivity_cohorts_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            summary = root / "mixed.tsv"
            summary.write_text(
                HEADER
                + "primary20\tprimary\tdeliciosa_a\tActinidia deliciosa\tA\tchromosome_partition\t100\t2\t10\t20.0\tdefined\n"
                + "rufa_sensitivity21\ta_rufa_sensitivity\trufa_chr\tActinidia rufa\tunphased\tchromosome_only_subset\t50\t1\t10\t10.0\tdefined\n",
                encoding="utf-8",
            )
            output = root / "must_not_exist"
            completed = self.run_figure(summary, output)
            self.assertEqual(completed.returncode, 2, completed.stderr)
            self.assertIn("separate bundles", completed.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
