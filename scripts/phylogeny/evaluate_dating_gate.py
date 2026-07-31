#!/usr/bin/env python3
"""Evaluate whether the frozen species topology has usable dating constraints."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

from Bio import Phylo


class GateError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def regular(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink() or resolved.stat().st_size == 0:
        raise GateError(f"missing, empty, or symlink file: {resolved}")
    return resolved


def binding(path: Path) -> dict[str, object]:
    source = regular(path)
    return {"basename": source.name, "bytes": source.stat().st_size, "sha256": sha256(source)}


def read_tsv(path: Path) -> list[dict[str, str]]:
    source = regular(path)
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise GateError(f"missing TSV header: {source}")
        rows = list(reader)
    return rows


def active_status(value: str) -> bool:
    normalized = value.strip().lower()
    return normalized == "pass" or normalized == "active" or normalized.startswith("active_")


def validate_active_secondary(
    rows: list[dict[str, str]], tree: object, tips: list[str], table_path: Path
) -> list[dict[str, object]]:
    required = {
        "constraint_id",
        "node_label",
        "descendant_a",
        "descendant_b",
        "minimum_ma",
        "maximum_ma",
        "distribution",
        "timetree_version",
        "retrieved_utc",
        "query_url",
        "contributing_studies",
        "raw_artifact_relative_path",
        "raw_artifact_sha256",
        "transformation_note",
        "status",
    }
    tip_set = set(tips)
    tip_by_name = {tip.name: tip for tip in tree.get_terminals()}
    repo_root = table_path.expanduser().resolve().parents[2]
    seen_ids: set[str] = set()
    seen_mrcas: set[tuple[str, ...]] = set()
    validated: list[dict[str, object]] = []
    for row in rows:
        if not required.issubset(row):
            raise GateError("secondary-constraint table is missing required columns")
        if not active_status(row["status"]):
            continue
        if any(not row[column].strip() for column in required):
            raise GateError(f"active secondary constraint has blank fields: {row['constraint_id']}")
        constraint_id = row["constraint_id"]
        if constraint_id in seen_ids:
            raise GateError(f"duplicate secondary constraint id: {constraint_id}")
        seen_ids.add(constraint_id)

        descendant_a = row["descendant_a"]
        descendant_b = row["descendant_b"]
        if descendant_a == descendant_b or {descendant_a, descendant_b} - tip_set:
            raise GateError(f"secondary constraint descendants do not match topology: {constraint_id}")
        mrca = tree.common_ancestor(tip_by_name[descendant_a], tip_by_name[descendant_b])
        mrca_signature = tuple(sorted(tip.name for tip in mrca.get_terminals()))
        if mrca_signature in seen_mrcas:
            raise GateError(f"multiple secondary constraints target the same MRCA: {constraint_id}")
        seen_mrcas.add(mrca_signature)

        minimum = float(row["minimum_ma"])
        maximum = float(row["maximum_ma"])
        if not 0 < minimum < maximum:
            raise GateError(f"invalid secondary constraint interval: {constraint_id}")
        try:
            retrieved = dt.datetime.fromisoformat(row["retrieved_utc"].replace("Z", "+00:00"))
        except ValueError as error:
            raise GateError(f"invalid retrieval timestamp: {constraint_id}") from error
        if retrieved.tzinfo is None:
            raise GateError(f"retrieval timestamp lacks timezone: {constraint_id}")

        relative = Path(row["raw_artifact_relative_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise GateError(f"unsafe raw artifact path: {constraint_id}")
        raw_path = regular(repo_root / relative)
        if sha256(raw_path) != row["raw_artifact_sha256"]:
            raise GateError(f"raw TimeTree checksum mismatch: {constraint_id}")
        raw = json.loads(raw_path.read_text(encoding="utf-8"))
        if float(raw["precomputed_ci_low"]) != minimum or float(raw["precomputed_ci_high"]) != maximum:
            raise GateError(f"TimeTree interval does not match raw response: {constraint_id}")
        expected_url = (
            "https://timetree.org/api/pairwise/"
            f"{raw['taxon_a_id']}/{raw['taxon_b_id']}/summaryjson"
        )
        if row["query_url"] != expected_url:
            raise GateError(f"TimeTree query URL does not match raw response: {constraint_id}")
        if not row["contributing_studies"].startswith(f"{raw['all_total']}_"):
            raise GateError(f"TimeTree contributing-estimate count mismatch: {constraint_id}")
        validated.append(
            {
                "constraint_id": constraint_id,
                "node_label": row["node_label"],
                "descendant_a": descendant_a,
                "descendant_b": descendant_b,
                "minimum_ma": minimum,
                "maximum_ma": maximum,
                "mrca_tip_count": len(mrca_signature),
                "raw_artifact": binding(raw_path),
            }
        )
    return validated


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology-freeze", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--calibrations", required=True, type=Path)
    parser.add_argument("--secondary-constraints", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        output = args.output_dir.expanduser().resolve()
        if output.exists():
            raise GateError(f"refusing to overwrite output: {output}")

        freeze_path = regular(args.topology_freeze)
        freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
        topology_path = regular(args.topology)
        if freeze.get("status") != "PASS":
            raise GateError("topology freeze is not PASS")
        accepted = freeze.get("accepted_topology", {})
        if accepted.get("sha256") != sha256(topology_path):
            raise GateError("topology does not match the frozen accepted topology")

        tree = Phylo.read(str(topology_path), "newick")
        tips = sorted(tip.name for tip in tree.get_terminals())
        if len(tips) != len(set(tips)) or any(not tip for tip in tips):
            raise GateError("topology tips are missing or duplicated")

        calibration_rows = read_tsv(args.calibrations)
        required_calibration_columns = {
            "calibration_id",
            "specimen_or_material",
            "activation_gate",
            "status",
        }
        if calibration_rows and not required_calibration_columns.issubset(calibration_rows[0]):
            raise GateError("calibration table is missing required columns")
        fossil_rows = [
            row
            for row in calibration_rows
            if row["specimen_or_material"].strip().lower() != "not a fossil calibration"
        ]
        active_fossils = [row for row in fossil_rows if active_status(row["status"])]

        secondary_rows = read_tsv(args.secondary_constraints)
        validated_secondary = validate_active_secondary(
            secondary_rows, tree, tips, args.secondary_constraints
        )
        active_secondary = [row for row in secondary_rows if active_status(row["status"])]

        secondary_status = "PASS" if active_secondary else "BLOCKED_NO_DECLARED_TIMETREE_CONSTRAINT"
        if active_fossils:
            primary_status = "PASS_FOSSIL_CALIBRATED"
            production_mcmctree_allowed = True
            calibration_basis = "active_fossil_constraints"
        elif active_secondary:
            primary_status = "PASS_SECONDARY_TIMETREE_CALIBRATION"
            production_mcmctree_allowed = True
            calibration_basis = "user_authorized_secondary_timetree_constraints"
        else:
            primary_status = "BLOCKED_NO_ACTIVE_CALIBRATION"
            production_mcmctree_allowed = False
            calibration_basis = "none"

        report = {
            "schema_version": 1,
            "workflow": "species_tree_dating_activation_gate",
            "status": primary_status,
            "production_mcmctree_allowed": production_mcmctree_allowed,
            "calibration_basis": calibration_basis,
            "secondary_timetree_sensitivity_status": secondary_status,
            "topology_freeze": binding(freeze_path),
            "accepted_topology": binding(topology_path),
            "tip_count": len(tips),
            "tips": tips,
            "calibrations": binding(args.calibrations),
            "fossil_row_count": len(fossil_rows),
            "active_fossil_count": len(active_fossils),
            "active_fossil_ids": [row["calibration_id"] for row in active_fossils],
            "inactive_fossils": [
                {
                    "calibration_id": row["calibration_id"],
                    "activation_gate": row["activation_gate"],
                    "status": row["status"],
                }
                for row in fossil_rows
                if row not in active_fossils
            ],
            "secondary_constraints": binding(args.secondary_constraints),
            "active_secondary_constraint_count": len(active_secondary),
            "active_secondary_constraint_ids": [row["constraint_id"] for row in active_secondary],
            "validated_secondary_constraints": validated_secondary,
            "decision": (
                "Run MCMCTree when either an active fossil set or a complete user-authorized "
                "TimeTree secondary-calibration set is present. TimeTree-calibrated output must "
                "remain explicitly labelled as secondary-calibrated, not fossil-calibrated."
            ),
        }

        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
        try:
            report_path = staging / "dating_gate.json"
            report_path.write_text(
                json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            (staging / "checksums.tsv").write_text(
                f"file\tsha256\ndating_gate.json\t{sha256(report_path)}\n",
                encoding="utf-8",
            )
            os.rename(staging, output)
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        print(f"{primary_status}\t{output}")
        return 0
    except (OSError, UnicodeError, json.JSONDecodeError, GateError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
