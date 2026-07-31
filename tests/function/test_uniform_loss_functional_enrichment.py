from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/function/run_uniform_loss_functional_enrichment.py"
SPEC = importlib.util.spec_from_file_location("uniform_functional_enrichment", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_bh_adjust_is_monotone_in_p_value_order():
    values = [0.04, 0.001, 0.02, 1.0]
    adjusted = MODULE.bh_adjust(values)
    ordered = sorted(zip(values, adjusted))
    assert all(left[1] <= right[1] for left, right in zip(ordered, ordered[1:]))
    assert all(value >= raw for raw, value in zip(values, adjusted))


def test_kegg_pathway_prefixes_collapse_to_one_identifier():
    assert MODULE.PATHWAY_RE.fullmatch("ko04141").group(1) == "04141"
    assert MODULE.PATHWAY_RE.fullmatch("map04141").group(1) == "04141"


def test_go_roots_are_explicitly_excluded():
    assert MODULE.GO_ROOTS == {"GO:0003674", "GO:0005575", "GO:0008150"}
