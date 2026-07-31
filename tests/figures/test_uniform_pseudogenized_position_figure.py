from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/figures/render_uniform_pseudogenized_positions.py"
SPEC = importlib.util.spec_from_file_location("uniform_position_figure", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_plot_schema_keeps_deleted_out_of_primary_panel():
    assert "pseudogenized_count" in MODULE.PLOT_COLUMNS
    assert "observed_locus_denominator" in MODULE.PLOT_COLUMNS
    assert "deleted" not in MODULE.PLOT_COLUMNS
    assert "uncertain" not in MODULE.PLOT_COLUMNS
