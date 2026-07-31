"""Tests for merging disjoint legacy/new expected-locus inputs."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "spatial" / "merge_deleted_locus_spatial_inputs.py"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(path: Path) -> dict[str, object]:
    return {"basename": path.name, "bytes": path.stat().st_size, "sha256": digest(path)}


def make_source(root: Path, name: str, unit: str, workflow: str) -> None:
    source = root / name
    source.mkdir()
    genome = source / f"{unit}.fa"
    gff = source / f"{unit}.gff3"
    genome.write_text(">Chr1\nAAAA\n")
    gff.write_text("##gff-version 3\nChr1\tt\tgene\t1\t4\t.\t+\t.\tID=x\n")
    calls = source / "positive_deleted_calls.tsv"
    calls.write_text(
        "assembly_unit_id\treference_gene_id\tclassification\tcallable\n"
        f"{unit}\tg1\tpositive_deleted\ttrue\n"
    )
    coordinates = source / "expected_deleted_locus_coordinates.tsv"
    coordinates.write_text(
        "assembly_unit_id\treference_gene_id\tclassification\tchromosome\texpected_locus_start_1based\texpected_locus_end_1based\tcoordinate_semantics\n"
        f"{unit}\tg1\tpositive_deleted\tChr1\t1\t4\texpected\n"
    )
    assemblies = source / "assembly_manifest.tsv"
    with assemblies.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "assembly_unit_id", "biological_species", "haplotype_or_subgenome", "assembly_scope",
                "genome", "gff", "genome_local_sha256", "gff_local_sha256",
            ),
            delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "assembly_unit_id": unit, "biological_species": "Actinidia test",
                "haplotype_or_subgenome": unit, "assembly_scope": "chromosome",
                "genome": genome.name, "gff": gff.name,
                "genome_local_sha256": digest(genome), "gff_local_sha256": digest(gff),
            }
        )
    report = {
        "status": "PASS", "workflow": workflow, "unit_count": 1,
        "outputs": {
            "positive_calls": binding(calls), "expected_locus_coordinates": binding(coordinates),
            "assembly_manifest": binding(assemblies),
        },
    }
    (source / "run_manifest.json").write_text(json.dumps(report) + "\n")
    with (source / "checksums.tsv").open("w") as handle:
        handle.write("file\tbytes\tsha256\n")
        for filename in (
            "positive_deleted_calls.tsv", "expected_deleted_locus_coordinates.tsv",
            "assembly_manifest.tsv", "run_manifest.json",
        ):
            path = source / filename
            handle.write(f"{filename}\t{path.stat().st_size}\t{digest(path)}\n")


class MergeDeletedLocusSpatialInputsTests(unittest.TestCase):
    def test_merges_disjoint_sources_and_rebinds_assets(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_source(
                root, "legacy", "u1", "legacy_conservative_deleted_expected_locus_spatial_inputs"
            )
            make_source(
                root, "new", "u2", "callable_positive_deleted_expected_locus_spatial_inputs"
            )
            sources = root / "sources.tsv"
            sources.write_text(
                "source_id\tinput_dir\texpected_unit_count\nlegacy\tlegacy\t1\nnew\tnew\t1\n"
            )
            completed = subprocess.run(
                [
                    sys.executable, str(SCRIPT), "--sources", str(sources),
                    "--data-root", str(root), "--expected-total-units", "2",
                    "--output-dir", str(root / "out"),
                ], text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads((root / "out" / "run_manifest.json").read_text())
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["unit_count"], 2)
            self.assertEqual(report["positive_deleted_count"], 2)
            with (root / "out" / "assembly_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["assembly_unit_id"] for row in rows}, {"u1", "u2"})
            for row in rows:
                self.assertTrue(((root / "out") / row["genome"]).resolve().is_file())


if __name__ == "__main__":
    unittest.main()
