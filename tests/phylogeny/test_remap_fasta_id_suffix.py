from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "phylogeny" / "remap_fasta_id_suffix.py"


def load_module():
    spec = importlib.util.spec_from_file_location("remap_fasta_id_suffix", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_suffix_remap_requires_exact_target_and_preserves_sequence(tmp_path: Path) -> None:
    module = load_module()
    source = tmp_path / "source.fa"
    expected = tmp_path / "expected.fa"
    output = tmp_path / "output.fa"
    source.write_text(">a description\nATGGCT\n>b\nATGAAA\n", encoding="utf-8")
    expected.write_text(">a.t\nMA\n>b.t\nMK\n", encoding="utf-8")

    provenance = module.remap(source, expected, output, ".t")

    assert output.read_text(encoding="utf-8") == (
        ">a.t description\nATGGCT\n>b.t\nATGAAA\n"
    )
    assert provenance["record_count"] == 2
    assert provenance["id_set_closure"] == "PASS"


def test_suffix_remap_fails_closed_on_nonmatching_target(tmp_path: Path) -> None:
    module = load_module()
    source = tmp_path / "source.fa"
    expected = tmp_path / "expected.fa"
    source.write_text(">a\nATG\n", encoding="utf-8")
    expected.write_text(">x.t\nM\n", encoding="utf-8")

    with pytest.raises(module.RemapError, match="do not equal expected IDs"):
        module.remap(source, expected, tmp_path / "output.fa", ".t")


def test_suffix_remap_can_explicitly_filter_source_extra_ids(tmp_path: Path) -> None:
    module = load_module()
    source = tmp_path / "source.fa"
    expected = tmp_path / "expected.fa"
    output = tmp_path / "output.fa"
    source.write_text(">a\nATG\n>extra\nAAA\n", encoding="utf-8")
    expected.write_text(">a.t\nM\n", encoding="utf-8")

    provenance = module.remap(
        source, expected, output, ".t", allow_source_extra=True
    )

    assert output.read_text(encoding="utf-8") == ">a.t\nATG\n"
    assert provenance["record_count"] == 1
    assert provenance["source_extra_id_count"] == 1
