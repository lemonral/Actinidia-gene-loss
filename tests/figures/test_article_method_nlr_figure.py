"""Regression test for the article-method NLR publication figure."""

from __future__ import annotations

import json
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "figures" / "render_article_method_nlr.py"


@unittest.skipUnless(
    importlib.util.find_spec("matplotlib") is not None
    and importlib.util.find_spec("PIL") is not None,
    "publication rendering dependencies are unavailable",
)
class ArticleMethodNlrFigureTest(unittest.TestCase):
    def test_renders_complete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            summary = root / "summary.tsv"
            fields = [
                "analysis_cohort",
                "cohort_role",
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "assembly_scope",
                "total_nlr_count",
                "article_retained_reference_nlr_count",
                "article_decayed_reference_nlr_loss_count",
                "article_deleted_reference_nlr_loss_count",
                "no_qualifying_translated_hit_count",
                "frameshift_supported_count",
                "inframe_stop_supported_count",
                "frameshift_and_stop_supported_count",
                "truncation_or_partial_alignment_candidate_count",
                "residual_sequence_mechanism_unresolved_count",
                "positive_reference_nlr_loss_count",
                "callable_reference_nlr_denominator",
                "positive_reference_nlr_loss_percentage",
                "percentage_status",
                "reference_nlr_universe_id",
            ]
            lines = ["\t".join(fields)]
            for index in range(23):
                row = {
                    "analysis_cohort": "article",
                    "cohort_role": "primary",
                    "assembly_unit_id": f"unit_{index:02d}",
                    "biological_species": "Actinidia alpha",
                    "haplotype_or_subgenome": "unphased",
                    "assembly_scope": "scope",
                    "total_nlr_count": "100",
                    "article_retained_reference_nlr_count": "8",
                    "article_decayed_reference_nlr_loss_count": "1",
                    "article_deleted_reference_nlr_loss_count": "1",
                    "no_qualifying_translated_hit_count": "1",
                    "frameshift_supported_count": "1",
                    "inframe_stop_supported_count": "0",
                    "frameshift_and_stop_supported_count": "0",
                    "truncation_or_partial_alignment_candidate_count": "0",
                    "residual_sequence_mechanism_unresolved_count": "0",
                    "positive_reference_nlr_loss_count": "2",
                    "callable_reference_nlr_denominator": "10",
                    "positive_reference_nlr_loss_percentage": "20.000000",
                    "percentage_status": "defined",
                    "reference_nlr_universe_id": "universe",
                }
                lines.append("\t".join(row[field] for field in fields))
            summary.write_text("\n".join(lines) + "\n", encoding="utf-8")
            manifest = root / "manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "status": "PASS_ARTICLE_METHOD_NLR_SUMMARY",
                        "article_shared_reference_nlrs_excluded": 2,
                        "article_nonshared_reference_nlrs": 10,
                    }
                )
                + "\n",
                encoding="utf-8",
            )
            output = root / "figure"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--unit-summary",
                    str(summary),
                    "--run-manifest",
                    str(manifest),
                    "--output-dir",
                    str(output),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            validation = json.loads(
                (output / "nlr_repertoire_and_loss_types.validation.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(validation["status"], "PASS_ARTICLE_METHOD_NLR_FIGURE")
            self.assertTrue(
                (output / "nlr_repertoire_and_loss_types.png").is_file()
            )


if __name__ == "__main__":
    unittest.main()
