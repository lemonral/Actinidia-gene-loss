from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "phylogeny" / "audit_protein_cds_pairs.py"


def load_module():
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("audit_protein_cds_pairs", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def test_audit_pair_filters_bad_cds_without_disabling_protein(tmp_path: Path) -> None:
    module = load_module()
    protein = tmp_path / "protein.fa"
    cds = tmp_path / "cds.fa"
    write(protein, ">good\nMA\n>bad\nMA\n")
    write(cds, ">good\nATGGCT\n>bad\nATGTAA\n")
    row = {
        "terminal_id": "sample",
        "protein_path": protein.name,
        "cds_path": cds.name,
        "use_for_orthofinder": "true",
        "use_for_codon_tree": "true",
    }

    summary, rejected = module.audit_pair(row, tmp_path)

    assert summary.protein_records == 2
    assert summary.codon_eligible_ids == 1
    assert summary.use_for_orthofinder == "true"
    assert summary.internal_stop_failures == 0  # terminal stops are permitted
    assert summary.translation_failures == 1
    assert {item["record_id"] for item in rejected} == {"bad"}


def test_audit_pair_reports_internal_stop(tmp_path: Path) -> None:
    module = load_module()
    protein = tmp_path / "protein.fa"
    cds = tmp_path / "cds.fa"
    write(protein, ">x\nMAA\n")
    write(cds, ">x\nATGTAGGCT\n")
    row = {
        "terminal_id": "sample",
        "protein_path": protein.name,
        "cds_path": cds.name,
        "use_for_orthofinder": "true",
        "use_for_codon_tree": "true",
    }

    summary, rejected = module.audit_pair(row, tmp_path)

    assert summary.internal_stop_failures == 1
    assert summary.codon_eligible_ids == 0
    assert summary.codon_gate == "FAIL"
    assert {item["reason"] for item in rejected} == {
        "internal_stop",
        "protein_cds_translation_mismatch",
    }
