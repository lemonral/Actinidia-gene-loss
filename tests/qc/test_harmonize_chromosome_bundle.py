from __future__ import annotations

import hashlib
import importlib.util
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "qc" / "harmonize_chromosome_bundle.py"
SPEC = importlib.util.spec_from_file_location("harmonize_chromosome_bundle", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class HarmonizeBundleWrapperTests(unittest.TestCase):
    def test_checksum_bundle_closes_exact_file_set(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            payload = root / "payload.txt"
            payload.write_text("verified\n", encoding="utf-8")
            digest = hashlib.sha256(payload.read_bytes()).hexdigest()
            (root / "checksums.tsv").write_text(
                f"file\tbytes\tsha256\npayload.txt\t{payload.stat().st_size}\t{digest}\n",
                encoding="utf-8",
            )
            bindings = MODULE.validate_checksum_bundle(root)
            self.assertEqual(bindings["payload.txt"]["sha256"], digest)
            (root / "untracked.txt").write_text("x", encoding="utf-8")
            with self.assertRaisesRegex(MODULE.PublicationError, "file closure"):
                MODULE.validate_checksum_bundle(root)

    def test_nonfinite_orientation_fraction_is_rejected(self) -> None:
        with self.assertRaisesRegex(MODULE.PublicationError, "finite"):
            MODULE.parse_exact_decimal("NaN", "dominant_fraction")
        self.assertEqual(MODULE.parse_exact_decimal("0.80", "x"), Decimal("0.80"))


if __name__ == "__main__":
    unittest.main()
