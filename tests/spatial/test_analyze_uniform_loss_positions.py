from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/spatial/analyze_uniform_loss_positions.py"
SPEC = importlib.util.spec_from_file_location("uniform_spatial", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_normalized_end_distance():
    assert MODULE.normalized_end_distance(1, 101, "x") == 0
    assert MODULE.normalized_end_distance(51, 101, "x") == 1
    assert MODULE.normalized_end_distance(101, 101, "x") == 0


def test_logistic_recovers_positive_position_slope():
    rng = np.random.default_rng(17)
    count = 4000
    x = rng.uniform(0, 1, count)
    unit = rng.integers(0, 4, count, dtype=np.int16)
    source = (unit >= 2).astype(np.int8)
    probability = 1 / (1 + np.exp(-(-2 + 2.0 * x)))
    y = rng.binomial(1, probability).astype(float)
    gene = np.arange(count, dtype=np.int32)
    result = MODULE.fit_logistic(unit, x, source, y, gene, interaction=False)
    slope = next(row for row in result["coefficients"] if row["term"] == "normalized_end_distance")
    assert slope["estimate_log_odds"] > 1.4
    assert slope["wald_p"] < 1e-8


def test_primary_model_uses_only_observed_pseudogenized_and_retained_loci():
    classes = np.asarray(["retained", "deleted", "pseudogenized", "retained"])
    sources = np.asarray([0, 0, 1, 1], dtype=np.int8)
    specifications = {
        name: (selected, interaction)
        for name, selected, interaction in MODULE.model_specifications(classes, sources)
    }
    selected, interaction = specifications["primary_pseudogenized_unit_fe_source_interaction"]
    assert selected.tolist() == [True, False, True, True]
    assert interaction is True
    deleted_selected, _ = specifications["deleted_expected_locus_vs_retained_sensitivity"]
    assert deleted_selected.tolist() == [True, True, False, True]
