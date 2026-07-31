"""Tests for callable copy-opportunity and loss-mode summaries."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "gene_loss" / "summarize_callable_copy_opportunities.py"


class CallableCopyOpportunitySummaryTests(unittest.TestCase):
    def test_complete_partial_and_copy_rates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            units = root / "units.tsv"
            columns = (
                "reference_gene_id", "assembly_unit_id", "biological_species", "aggregation_rule",
                "classification", "callable", "evidence_state", "positive_call", "confident_negative",
            )
            rows = [
                ("g1", "uA", "Species one", "deleted", "true", "positive", "true", "false"),
                ("g1", "uB", "Species one", "deleted", "true", "positive", "true", "false"),
                ("g2", "uA", "Species one", "deleted", "true", "positive", "true", "false"),
                ("g2", "uB", "Species one", "retained", "true", "confident_negative", "false", "true"),
                ("g1", "u2", "Species two", "deleted", "true", "positive", "true", "false"),
                ("g2", "u2", "Species two", "uncertain", "false", "uncertain", "false", "false"),
            ]
            with units.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for gene, unit, species, classification, callable_value, state, positive, negative in rows:
                    writer.writerow(
                        {
                            "reference_gene_id": gene, "assembly_unit_id": unit,
                            "biological_species": species, "aggregation_rule": "all_units_positive",
                            "classification": classification, "callable": callable_value,
                            "evidence_state": state, "positive_call": positive,
                            "confident_negative": negative,
                        }
                    )
            species = root / "species.tsv"
            species.write_text(
                "reference_gene_id\tbiological_species\tspecies_gene_status\tassembly_unit_count\n"
                "g1\tSpecies one\tpositive_complete\t2\n"
                "g2\tSpecies one\tpositive_partial\t2\n"
                "g1\tSpecies two\tpositive_complete\t1\n"
                "g2\tSpecies two\tuncertain\t1\n",
                encoding="utf-8",
            )
            shared = root / "shared.tsv"
            shared.write_text("reference_gene_id\ng1\n", encoding="utf-8")
            output = root / "output"
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--unit-calls", str(units),
                    "--species-matrix", str(species), "--shared-genes", str(shared),
                    "--output-dir", str(output),
                ],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with (output / "species_loss_mode_and_copy_opportunity_summary.tsv").open(newline="") as handle:
                result = {row["biological_species"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(result["Species one"]["complete_loss_gene_count"], "1")
            self.assertEqual(result["Species one"]["partial_homeolog_loss_gene_count"], "1")
            self.assertEqual(result["Species one"]["callable_copy_opportunity_loss_rate"], "0.75")
            self.assertEqual(result["Species one"]["nonshared_callable_copy_opportunity_loss_rate"], "0.5")
            self.assertEqual(result["Species two"]["nonshared_callable_copy_opportunity_loss_rate"], "NA")
            report = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["shared_positive_complete_gene_count"], 1)


if __name__ == "__main__":
    unittest.main()
