from __future__ import annotations

import gzip
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "downstream" / "prepare_unit_loss_evidence_scaffold.py"


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")


def test_unit_scaffold_assigns_parallel_and_repeated_events(tmp_path: Path) -> None:
    metadata = tmp_path / "units.tsv"
    write_tsv(
        metadata,
        [
            {
                "assembly_unit_id": "u1",
                "biological_species": "Species one",
                "haplotype_or_subgenome": "A",
                "include": "true",
            },
            {
                "assembly_unit_id": "u2",
                "biological_species": "Species one",
                "haplotype_or_subgenome": "B",
                "include": "true",
            },
            {
                "assembly_unit_id": "u3",
                "biological_species": "Species two",
                "haplotype_or_subgenome": "",
                "include": "true",
            },
            {
                "assembly_unit_id": "u4",
                "biological_species": "Species three",
                "haplotype_or_subgenome": "",
                "include": "true",
            },
        ],
    )
    tip_map = tmp_path / "tip_map.tsv"
    write_tsv(
        tip_map,
        [
            {"tree_tip": "T1", "biological_species": "Species one", "include": "true"},
            {"tree_tip": "T2", "biological_species": "Species two", "include": "true"},
            {"tree_tip": "T3", "biological_species": "Species three", "include": "true"},
        ],
    )
    tree = tmp_path / "tree.tre"
    tree.write_text("(T1:1,(T2:1,T3:1):1);\n", encoding="utf-8")

    states = {
        "g1": {"u1": "decayed", "u2": "decayed", "u3": "decayed", "u4": "decayed"},
        "g2": {"u1": "decayed", "u2": "retained", "u3": "retained", "u4": "retained"},
        "g3": {"u1": "decayed", "u2": "decayed", "u3": "retained", "u4": "retained"},
        "g4": {"u1": "deleted", "u2": "retained", "u3": "deleted", "u4": "retained"},
        "g5": {
            "u1": "decayed",
            "u2": "retained",
            "u3": "retained",
            "u4": "not_called_loss",
        },
    }
    causes = {
        ("g1", "u1"): "frameshift_supported",
        ("g1", "u2"): "stop_supported",
        ("g1", "u3"): "frameshift_and_stop_supported",
        ("g1", "u4"): "local_sequence_no_explicit_coding_disruption",
        ("g2", "u1"): "frameshift_supported",
        ("g3", "u1"): "stop_supported",
        ("g3", "u2"): "stop_supported",
        ("g5", "u1"): "frameshift_and_stop_supported",
    }
    rows: list[dict[str, object]] = []
    for gene, unit_states in states.items():
        for unit, state in unit_states.items():
            cause = causes.get((gene, unit), "not_applicable_retained")
            positive = state in {"decayed", "deleted"}
            if state == "deleted":
                cause = "no_qualifying_genomewide_tblastx_hit"
            elif state == "not_called_loss":
                cause = "not_called_outside_historical_scope"
            confirmed = cause in {
                "frameshift_supported",
                "stop_supported",
                "frameshift_and_stop_supported",
            }
            rows.append(
                {
                    "reference_gene_id": gene,
                    "assembly_unit_id": unit,
                    "manuscript_classification": state,
                    "manuscript_positive_loss": str(positive).lower(),
                    "uniform_classification": "pseudogenized" if confirmed else state,
                    "uniform_evidence_reason": "fixture",
                    "refined_decayed_cause": cause,
                    "refined_cause_evidence_level": (
                        "explicit_coding_disruption" if confirmed else "fixture"
                    ),
                    "query_coverage": "1" if confirmed else "",
                    "exact_alignment_identity": "0.9" if confirmed else "",
                    "alignment_score": "500" if confirmed else "",
                    "frameshift_events": (
                        "1"
                        if cause
                        in {"frameshift_supported", "frameshift_and_stop_supported"}
                        else "0" if confirmed else ""
                    ),
                    "inframe_stop_codons": (
                        "1"
                        if cause
                        in {"stop_supported", "frameshift_and_stop_supported"}
                        else "0" if confirmed else ""
                    ),
                }
            )
    matrix = tmp_path / "matrix.tsv.gz"
    with gzip.open(matrix, "wt", encoding="utf-8", newline="") as handle:
        pd.DataFrame(rows).to_csv(
            handle,
            sep="\t",
            index=False,
            lineterminator="\n",
        )

    output = tmp_path / "output"
    command = [
        sys.executable,
        str(SCRIPT),
        "--manuscript-matrix",
        str(matrix),
        "--unit-metadata",
        str(metadata),
        "--tip-map",
        str(tip_map),
        "--time-tree",
        str(tree),
        "--output-dir",
        str(output),
        "--expected-units",
        "4",
        "--expected-reference-genes",
        "5",
        "--expected-lineages",
        "3",
        "--expected-positive-rows",
        "10",
        "--expected-shared-genes",
        "1",
        "--expected-decayed-frameshift-only",
        "2",
        "--expected-decayed-stop-only",
        "3",
        "--expected-decayed-both",
        "2",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    assert completed.returncode == 0, completed.stderr

    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["status"] == "PASS_UNIT_RESOLVED_ARTICLE_LOSS_SCAFFOLD"
    assert manifest["counts"]["scaffold_terminal_nodes"] == 4
    assert manifest["counts"]["scaffold_event_rows"] == 5
    summary = pd.read_csv(output / "scaffold_pattern_summary.tsv", sep="\t")
    observed = dict(zip(summary["scaffold_pattern"], summary["reference_gene_count"]))
    assert observed == {
        "no_loss": 0,
        "single_terminal_event": 1,
        "single_internal_event": 2,
        "repeated_independent_events": 1,
        "ambiguous_not_called": 1,
    }
    mechanisms = pd.read_csv(
        output / "unit_decayed_mechanism_summary.tsv",
        sep="\t",
    )
    assert mechanisms["frameshift_only"].sum() == 2
    assert mechanisms["inframe_stop_only"].sum() == 3
    assert mechanisms["frameshift_and_inframe_stop"].sum() == 2
