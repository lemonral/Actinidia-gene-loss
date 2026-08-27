from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "statistics" / "compare_ploidy_loss_rates.py"
INPUT = ROOT / "results" / "tables" / "ploidy_comparison" / "unit_loss_summary.tsv"


def load_module():
    specification = importlib.util.spec_from_file_location(
        "compare_ploidy_loss_rates",
        SCRIPT,
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def test_primary_exact_ploidy_comparison() -> None:
    module = load_module()
    observations = module.read_observations(INPUT, expected_genomes=23)
    polyploid = [
        row.positive_loss_rate
        for row in observations
        if row.ploidy_group == "polyploid"
    ]
    diploid = [
        row.positive_loss_rate
        for row in observations
        if row.ploidy_group == "diploid"
    ]
    result = module.exact_tests(polyploid, diploid)

    assert len(polyploid) == 11
    assert len(diploid) == 12
    assert result.mann_whitney_u_polyploid == 110
    assert result.mann_whitney_u_complement == 22
    assert result.mann_whitney_u_reported == 22
    assert result.permutation_assignments == math.comb(23, 11)
    assert math.isclose(
        result.mann_whitney_exact_p,
        0.00562245669258726,
        rel_tol=0,
        abs_tol=1e-15,
    )
    assert math.isclose(
        result.permutation_exact_p,
        0.004138074874378549,
        rel_tol=0,
        abs_tol=1e-15,
    )
