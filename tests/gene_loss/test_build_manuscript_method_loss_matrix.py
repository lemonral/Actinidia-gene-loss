from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/gene_loss/build_manuscript_method_loss_matrix.py"
SPEC = importlib.util.spec_from_file_location("manuscript_method_matrix", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def uniform(**updates: str) -> dict[str, str]:
    row = {
        "frameshift_events": "0",
        "inframe_stop_codons": "0",
        "disruption_supported": "false",
        "callable": "true",
        "classification": "uncertain",
        "evidence_reason": "local_sequence_without_disruptive_event",
    }
    row.update(updates)
    return row


def miniprot(**updates: str) -> dict[str, str]:
    row = {
        "qualifying_local_alignment": "true",
        "query_length_aa": "100",
        "query_aligned_start_0based": "0",
        "query_aligned_end_0based_exclusive": "100",
        "query_coverage": "1",
        "frameshift_events": "0",
        "inframe_stop_codons": "0",
    }
    row.update(updates)
    return row


def cause(row: dict[str, str], state: dict[str, str] | None = None):
    return MODULE.refined_decayed_cause(
        "decayed", row, state, terminal_missing_fraction=0.20
    )


def test_explicit_disruptions_are_separate_from_manuscript_classification():
    assert cause(uniform(frameshift_events="1", disruption_supported="true"))[0] == "frameshift_supported"
    assert cause(uniform(inframe_stop_codons="2", disruption_supported="true"))[0] == "stop_supported"
    assert cause(uniform(frameshift_events="1", inframe_stop_codons="2", disruption_supported="true"))[0] == "frameshift_and_stop_supported"
    assert cause(uniform(frameshift_events="1"))[0] == "frameshift_or_stop_below_strict_quality_gate"


def test_terminal_truncation_is_candidate_only():
    label, level = cause(uniform(), miniprot(query_aligned_start_0based="25"))
    assert label == "n_terminal_alignment_truncation_candidate"
    assert level == "alignment_candidate_only"
    label, level = cause(uniform(), miniprot(query_aligned_end_0based_exclusive="75"))
    assert label == "c_terminal_alignment_truncation_candidate"
    assert level == "alignment_candidate_only"
    label, _ = cause(
        uniform(),
        miniprot(query_aligned_start_0based="25", query_aligned_end_0based_exclusive="75"),
    )
    assert label == "both_terminal_alignment_truncation_candidate"


def test_original_classes_are_not_overwritten_by_refinement():
    row = uniform()
    assert MODULE.refined_decayed_cause("retained", row, None, terminal_missing_fraction=0.20)[0] == "not_applicable_retained"
    assert MODULE.refined_decayed_cause("deleted", row, None, terminal_missing_fraction=0.20)[0] == "no_qualifying_genomewide_tblastx_hit"
    assert MODULE.refined_decayed_cause(
        "deleted",
        uniform(frameshift_events="1", disruption_supported="true"),
        None,
        terminal_missing_fraction=0.20,
    )[0] == "frameshift_supported"
    assert MODULE.refined_decayed_cause("not_called_loss", row, None, terminal_missing_fraction=0.20)[0] == "not_called_outside_historical_scope"
