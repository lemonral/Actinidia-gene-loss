from __future__ import annotations

import argparse
import csv
import tempfile
import unittest
from pathlib import Path

from scripts.qc.build_target_asset_registry import RegistryError, run


class TargetAssetRegistryTests(unittest.TestCase):
    def test_exact_three_role_registry_and_no_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assets = {}
            for role in ("genome", "gff", "protein"):
                path = root / f"unit.{role}"
                path.write_text(f"{role}\n", encoding="utf-8")
                assets[role] = path
            output = root / "target_assets.tsv"
            arguments = argparse.Namespace(
                assembly_unit_id="act_test",
                target_scope_id="act_test.scope_v1",
                genome=assets["genome"],
                gff=assets["gff"],
                protein=assets["protein"],
                output=output,
            )
            self.assertEqual(run(arguments), output.resolve())
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["asset_role"] for row in rows], ["genome", "gff", "protein"])
            self.assertEqual({row["status"] for row in rows}, {"verified"})
            with self.assertRaisesRegex(RegistryError, "Refusing to overwrite"):
                run(arguments)


if __name__ == "__main__":
    unittest.main()
