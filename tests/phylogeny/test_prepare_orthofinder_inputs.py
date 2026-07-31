from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "phylogeny" / "prepare_orthofinder_inputs.py"


def load_module():
    spec = importlib.util.spec_from_file_location("prepare_orthofinder_inputs", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_prepare_binds_passed_audit_and_symlinks_proteome(tmp_path: Path) -> None:
    module = load_module()
    data = tmp_path / "data"
    data.mkdir()
    protein = data / "x.faa"
    cds = data / "x.cds"
    protein.write_text(">g\nMA\n", encoding="utf-8")
    cds.write_text(">g\nATGGCT\n", encoding="utf-8")
    manifest = tmp_path / "pairs.tsv"
    manifest.write_text(
        "terminal_id\tprotein_path\tcds_path\tuse_for_orthofinder\tuse_for_codon_tree\trole\n"
        "x\tx.faa\tx.cds\ttrue\ttrue\tingroup\n",
        encoding="utf-8",
    )
    audit = tmp_path / "audit"
    audit.mkdir()
    summary = audit / "protein_cds_pair_audit.tsv"
    fields = [
        "terminal_id", "protein_sha256", "protein_records", "codon_eligible_ids", "codon_gate"
    ]
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "terminal_id": "x",
                "protein_sha256": sha(protein),
                "protein_records": "1",
                "codon_eligible_ids": "1",
                "codon_gate": "PASS",
            }
        )
    rejected = audit / "x.rejected_ids.tsv"
    rejected.write_text("terminal_id\trecord_id\treason\tdetail\n", encoding="utf-8")
    outputs = [summary, rejected]
    (audit / "provenance.json").write_text(
        json.dumps(
            {
                "manifest_sha256": sha(manifest),
                "outputs": [
                    {"basename": path.name, "sha256": sha(path)} for path in outputs
                ],
            }
        ),
        encoding="utf-8",
    )

    output = module.prepare(manifest, data, audit, tmp_path / "prepared")

    assert (output / "proteomes" / "x.faa").is_symlink()
    binding = json.loads((output / "input_binding.json").read_text(encoding="utf-8"))
    assert binding["status"] == "PASS"
    assert binding["terminal_count"] == 1


def test_prepare_rejects_audit_manifest_mismatch(tmp_path: Path) -> None:
    module = load_module()
    manifest = tmp_path / "pairs.tsv"
    manifest.write_text(
        "terminal_id\tprotein_path\tcds_path\tuse_for_orthofinder\tuse_for_codon_tree\n"
        "x\tx\tx\ttrue\ttrue\n",
        encoding="utf-8",
    )
    audit = tmp_path / "audit"
    audit.mkdir()
    (audit / "protein_cds_pair_audit.tsv").write_text("terminal_id\n", encoding="utf-8")
    (audit / "provenance.json").write_text(
        json.dumps({"manifest_sha256": "0" * 64, "outputs": []}),
        encoding="utf-8",
    )

    with pytest.raises(module.InputError, match="does not bind"):
        module.prepare(manifest, tmp_path, audit, tmp_path / "prepared")
