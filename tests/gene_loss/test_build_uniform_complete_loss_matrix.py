from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/gene_loss/build_uniform_complete_loss_matrix.py"
SPEC = importlib.util.spec_from_file_location("uniform_matrix", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def checksum(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bind(path: Path):
    return {"basename": path.name, "bytes": path.stat().st_size, "sha256": checksum(path)}


def bundle(root: Path, unit: str, table: str, rows: list[dict[str, str]]) -> Path:
    root.mkdir(parents=True)
    table_path = root / table
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    (root / "run_manifest.json").write_text(json.dumps({"status": "PASS", "unit": unit, "metrics": {}}))
    with (root / "checksums.tsv").open("w") as handle:
        handle.write("file\tbytes\tsha256\n")
        handle.write(f"{table}\t{table_path.stat().st_size}\t{checksum(table_path)}\n")
    return root


def test_strict_disruption_gate_and_complete_grid(tmp_path: Path):
    data = tmp_path / "data"; data.mkdir()
    protein = data / "protein.faa"
    protein.write_text(">g1\nMAA\n>g2\nMAA\n>g3\nMAA\n>g4\nMAA\n>g5\nMAA\n")
    synorth = data / "synorth.tsv"
    synorth.write_text("t1\tChr1\t10\t20\tg1\tRef\t1\t9\n")
    candidate_rows = []
    state_rows = []
    for gene, raw_class, cov, ident, score, fs, st, qualifying in [
        ("g2", "deleted", "", "", "", "", "", "false"),
        ("g3", "pseudogenized", "0.9", "0.8", "200", "1", "0", "true"),
        ("g4", "pseudogenized", "0.6", "0.8", "200", "1", "0", "true"),
    ]:
        candidate_rows.append({"unit": "u", "reference_gene": gene, "callable": "true", "callability_reason": "callable", "target_chromosome": "Chr1", "target_interval_start_1based": "30", "target_interval_end_1based": "90"})
        state_rows.append({"unit": "u", "reference_gene": gene, "callable": "true", "classification": raw_class, "qualifying_local_alignment": qualifying, "alignment_target_start_1based": "40" if qualifying == "true" else "", "alignment_target_end_1based": "60" if qualifying == "true" else "", "query_coverage": cov, "exact_alignment_identity": ident, "alignment_score": score, "frameshift_events": fs, "inframe_stop_codons": st, "evidence_reason": "raw"})
    bundle(data / "cand", "u", "candidates.tsv", candidate_rows)
    bundle(data / "uniform", "u", "uniform_candidate_loss_states.tsv", state_rows)
    config = tmp_path / "config.tsv"
    config.write_text("\t".join(MODULE.CONFIG_COLUMNS) + "\n" + "u\tnew\tsynorth.tsv\t5\tcand\tuniform\n")
    output = tmp_path / "out"
    import sys
    old = sys.argv
    try:
        sys.argv = [str(SCRIPT), "--config", str(config), "--data-root", str(data), "--reference-protein", str(protein), "--output-dir", str(output)]
        assert MODULE.main() == 0
    finally:
        sys.argv = old
    with (output / "complete_unit_loss_matrix.tsv").open() as handle:
        rows = {row["reference_gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["g1"]["classification"] == "retained"
    assert rows["g2"]["classification"] == "deleted"
    assert rows["g3"]["classification"] == "pseudogenized"
    assert rows["g3"]["disruption_supported"] == "true"
    assert rows["g4"]["classification"] == "uncertain"
    assert rows["g5"]["classification"] == "uncertain"
    assert rows["g5"]["callable"] == "false"
    assert sum(row["positive_loss"] == "true" for row in rows.values()) == 2
    with (output / "strict_disruption_calls.tsv").open() as handle:
        disruptions = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["reference_gene_id"] for row in disruptions] == ["g3"]
