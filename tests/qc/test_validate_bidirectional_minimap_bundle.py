from __future__ import annotations

import hashlib
import json
from pathlib import Path

from scripts.qc.validate_bidirectional_minimap_bundle import (
    EXPECTED_ARGV,
    EXPECTED_COMPARISONS,
    ValidationError,
    validate,
)


def record(path: Path) -> dict[str, object]:
    data = path.read_bytes()
    return {"basename": path.name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def fixture(tmp_path: Path) -> tuple[Path, dict[str, Path]]:
    inputs = {role: tmp_path / f"{role}.fa" for role in ("target", "hy4a", "hy4p")}
    for role, path in inputs.items():
        path.write_text(f">{role}\nACGT\n", encoding="utf-8")
    comparisons = {}
    for name, (query, reference) in EXPECTED_COMPARISONS.items():
        paf = tmp_path / f"{name}.paf"
        stderr = tmp_path / f"{name}.stderr.log"
        paf.write_text("q\t4\t0\t4\t+\tt\t4\t0\t4\t4\t4\t60\n", encoding="utf-8")
        stderr.write_text("ok\n", encoding="utf-8")
        comparisons[name] = {
            "exit_code": 0,
            "finished_at_utc": "2026-01-01T00:00:00+00:00",
            "query_role": query,
            "reference_role": reference,
            "paf": record(paf),
            "stderr": record(stderr),
        }
    status = {
        "schema_version": 1,
        "workflow": "bidirectional_chromosome_minimap",
        "status": "completed",
        "unit": "unit",
        "started_at_utc": "2026-01-01T00:00:00+00:00",
        "finished_at_utc": "2026-01-01T01:00:00+00:00",
        "minimap2_version": "2.28-r1209",
        "fixed_argv": EXPECTED_ARGV,
        "inputs": {role: record(path) for role, path in inputs.items()},
        "comparisons": comparisons,
    }
    status_path = tmp_path / "status.json"
    status_path.write_text(json.dumps(status), encoding="utf-8")
    return status_path, inputs


def test_valid_bundle_passes(tmp_path: Path) -> None:
    status, inputs = fixture(tmp_path)
    result = validate(
        status_path=status,
        target_genome=inputs["target"],
        hy4a_genome=inputs["hy4a"],
        hy4p_genome=inputs["hy4p"],
        expected_unit="unit",
    )
    assert result["status"] == "PASS"
    assert len(result["comparisons"]) == 4


def test_changed_paf_fails(tmp_path: Path) -> None:
    status, inputs = fixture(tmp_path)
    (tmp_path / "target_to_hy4a.paf").write_text("changed\n", encoding="utf-8")
    try:
        validate(
            status_path=status,
            target_genome=inputs["target"],
            hy4a_genome=inputs["hy4a"],
            hy4p_genome=inputs["hy4p"],
            expected_unit="unit",
        )
    except ValidationError as error:
        assert "mismatch" in str(error)
    else:
        raise AssertionError("changed PAF was accepted")
