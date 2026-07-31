from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/figures/render_uniform_loss_functional_enrichment.py"
SPEC = importlib.util.spec_from_file_location("uniform_functional_figure", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)

CATEGORY_SCRIPT = ROOT / "scripts/figures/render_uniform_loss_functional_categories.py"
CATEGORY_SPEC = importlib.util.spec_from_file_location(
    "uniform_functional_categories", CATEGORY_SCRIPT
)
CATEGORY_MODULE = importlib.util.module_from_spec(CATEGORY_SPEC)
assert CATEGORY_SPEC.loader is not None
CATEGORY_SPEC.loader.exec_module(CATEGORY_MODULE)


def test_primary_and_strict_pseudogenized_are_separate_foregrounds():
    assert dict(MODULE.FOREGROUNDS) == {
        "pooled_nonshared_combined": "Deleted + strict pseudogenized",
        "pooled_nonshared_strict_pseudogenized": "Strict pseudogenized only",
    }


def test_figure_does_not_treat_uncertain_as_loss():
    serialized = " ".join(label for _, label in MODULE.FOREGROUNDS)
    assert "uncertain" not in serialized.lower()


def test_functional_categories_are_explicit_and_nonoverlapping():
    assert [row[0] for row in CATEGORY_MODULE.CATEGORIES] == [
        "GO_BP",
        "GO_MF",
        "GO_CC",
        "KEGG_KO",
        "KEGG_PATHWAY",
    ]
    assert CATEGORY_MODULE.functional_category(
        {"ontology": "GO", "go_namespace": "biological_process"}
    ) == "GO_BP"
    assert CATEGORY_MODULE.functional_category(
        {"ontology": "KEGG_PATHWAY", "go_namespace": ""}
    ) == "KEGG_PATHWAY"


def test_category_summary_closes_expected_pooled_counts():
    rows, validation, _ = CATEGORY_MODULE.prepare(
        ROOT / "results/tables/functional_enrichment_uniform"
    )
    assert len(rows) == 140
    assert validation["status"] == "pass"
    pooled = {
        (row["evidence_mode"], row["category_id"]): row["significant_term_count"]
        for row in rows
        if row["foreground_scope"] == "pooled"
    }
    assert pooled[("deleted_plus_strict_pseudogenized", "GO_BP")] == 12
    assert pooled[("deleted_plus_strict_pseudogenized", "GO_MF")] == 17
    assert pooled[("deleted_plus_strict_pseudogenized", "GO_CC")] == 0
    assert pooled[("all_required_units_strict_pseudogenized", "KEGG_KO")] == 21
    assert pooled[("all_required_units_strict_pseudogenized", "KEGG_PATHWAY")] == 32
    displays = {row["biological_species"]: row["display_label"] for row in rows}
    assert r"\mathrm{unphased}" not in displays["Actinidia rufa"]
    assert r"\mathrm{A}" in displays["Actinidia x zhejiangensis parental lineage A"]
