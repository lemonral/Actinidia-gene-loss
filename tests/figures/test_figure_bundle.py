"""Tests for atomic, privacy-safe publication figure bundles."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # matplotlib is an optional project dependency
    matplotlib = None
    plt = None

from geneloss_repro.figure_bundle import write_figure_bundle


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FigureBundleTest(unittest.TestCase):
    class StubFigure:
        """Exercise all non-rendering bundle logic without optional matplotlib."""

        def savefig(self, path, *, format, **kwargs):
            signature = b"\x89PNG\r\n\x1a\n" if format == "png" else b"%PDF-1.4\n"
            Path(path).write_bytes(signature + b"test figure\n")

    def make_figure(self):
        if plt is None:
            raise RuntimeError("matplotlib is unavailable")
        figure, axis = plt.subplots()
        axis.plot([1, 2], [3, 4])
        return figure

    @unittest.skipUnless(matplotlib is not None, "optional matplotlib is not installed")
    def test_json_bundle_is_complete_public_and_path_free(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "analysis_summary.tsv"
            source.write_text("sample\tcount\nA\t3\n", encoding="utf-8")
            figure = self.make_figure()
            try:
                bundle = write_figure_bundle(
                    figure=figure,
                    output_dir=root / "publication" / "R6",
                    basename="R6_nlr",
                    plot_rows=[
                        {"sample": "A. deliciosa A", "total": 42, "percentage": 12.5},
                        {"sample": "A. deliciosa B", "total": 40, "percentage": 10.0},
                    ],
                    plot_columns=["sample", "total", "percentage"],
                    caption="Figure R6. NLR loss counts and percentages.",
                    validation={"status": "pass", "checks": ["denominator", "checksum"]},
                    input_paths=[source],
                    dpi=100,
                )
            finally:
                plt.close(figure)

            expected_names = {
                "R6_nlr.png",
                "R6_nlr.pdf",
                "R6_nlr.plot_data.tsv",
                "R6_nlr.caption.txt",
                "R6_nlr.validation.json",
                "R6_nlr.manifest.json",
            }
            self.assertEqual({path.name for path in bundle.directory.iterdir()}, expected_names)
            self.assertTrue(bundle.png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
            self.assertTrue(bundle.pdf.read_bytes().startswith(b"%PDF"))
            self.assertEqual(
                bundle.plot_data.read_text(encoding="utf-8"),
                "sample\ttotal\tpercentage\n"
                "A. deliciosa A\t42\t12.5\n"
                "A. deliciosa B\t40\t10.0\n",
            )
            self.assertEqual(
                bundle.caption.read_text(encoding="utf-8"),
                "Figure R6. NLR loss counts and percentages.\n",
            )
            self.assertEqual(json.loads(bundle.validation.read_text())["status"], "pass")

            manifest_text = bundle.manifest.read_text(encoding="utf-8")
            manifest = json.loads(manifest_text)
            self.assertNotIn(str(root), manifest_text)
            self.assertNotIn("command", manifest_text.lower())
            self.assertEqual(
                manifest["inputs"],
                [{"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}],
            )
            for entry in manifest["outputs"]:
                path = bundle.directory / entry["basename"]
                self.assertEqual(entry["sha256"], sha256(path))
                self.assertEqual(entry["bytes"], path.stat().st_size)
                self.assertEqual(Path(entry["basename"]).name, entry["basename"])
            self.assertEqual(bundle.directory.stat().st_mode & 0o777, 0o755)
            for path in bundle.directory.iterdir():
                self.assertEqual(path.stat().st_mode & 0o777, 0o644)

    def test_tsv_validation_can_publish_into_existing_empty_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "R7"
            output.mkdir()
            bundle = write_figure_bundle(
                figure=self.StubFigure(),
                output_dir=output,
                basename="R7_spatial",
                plot_rows=[{"bin": 1, "losses": 3}],
                plot_columns=["bin", "losses"],
                caption="Figure R7. Spatial distribution of gene loss.",
                validation=[{"check": "coverage", "status": "pass"}],
                validation_format="tsv",
                validation_columns=["check", "status"],
                dpi=72,
            )
            self.assertEqual(
                bundle.validation.read_text(encoding="utf-8"),
                "check\tstatus\ncoverage\tpass\n",
            )
            self.assertEqual(
                json.loads(bundle.manifest.read_text())["outputs"][-1]["role"],
                "validation_tsv",
            )

    def test_nonempty_destination_and_bad_basename_fail_without_changes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            output = root / "bundle"
            output.mkdir()
            sentinel = output / "keep.txt"
            sentinel.write_text("do not replace", encoding="utf-8")
            arguments = dict(
                figure=self.StubFigure(),
                output_dir=output,
                basename="R6",
                plot_rows=[{"x": 1}],
                plot_columns=["x"],
                caption="Figure R6. Caption.",
                validation={"status": "pass"},
                dpi=72,
            )
            with self.assertRaises(FileExistsError):
                write_figure_bundle(**arguments)
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "do not replace")
            for bad in ("../R6", "R6.pdf", "R 6", "_R6", ""):
                with self.subTest(basename=bad), self.assertRaises(ValueError):
                    write_figure_bundle(
                        **{
                            **arguments,
                            "output_dir": root / bad.replace("/", "_"),
                            "basename": bad,
                        }
                    )

    def test_schema_duplicate_inputs_and_nonfinite_values_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            first = root / "one" / "same.tsv"
            second = root / "two" / "same.tsv"
            first.parent.mkdir()
            second.parent.mkdir()
            first.write_text("one", encoding="utf-8")
            second.write_text("two", encoding="utf-8")
            common = dict(
                figure=self.StubFigure(),
                output_dir=root / "bundle",
                basename="R6",
                plot_columns=["x"],
                caption="Figure R6. Caption.",
                validation={"status": "pass"},
                dpi=72,
            )
            with self.assertRaises(ValueError):
                write_figure_bundle(
                    **common,
                    plot_rows=[{"x": 1, "unexpected": 2}],
                )
            self.assertFalse((root / "bundle").exists())
            with self.assertRaises(ValueError):
                write_figure_bundle(
                    **common,
                    plot_rows=[{"x": 1}],
                    input_paths=[first, second],
                )
            self.assertFalse((root / "bundle").exists())
            with self.assertRaises(ValueError):
                write_figure_bundle(
                    **{**common, "validation": {"value": float("nan")}},
                    plot_rows=[{"x": 1}],
                )
            self.assertFalse((root / "bundle").exists())

    def test_figure_failure_does_not_publish_partial_bundle(self) -> None:
        class FailingFigure:
            calls = 0

            def savefig(self, path, **kwargs):
                self.calls += 1
                Path(path).write_bytes(b"partial")
                if self.calls == 2:
                    raise RuntimeError("simulated PDF failure")

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            output = root / "bundle"
            with self.assertRaises(RuntimeError):
                write_figure_bundle(
                    figure=FailingFigure(),
                    output_dir=output,
                    basename="R6",
                    plot_rows=[{"x": 1}],
                    plot_columns=["x"],
                    caption="Figure R6. Caption.",
                    validation={"status": "pass"},
                    dpi=72,
                )
            self.assertFalse(output.exists())
            self.assertEqual(list(root.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
