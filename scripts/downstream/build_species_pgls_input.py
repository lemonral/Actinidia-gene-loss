#!/usr/bin/env python3
"""Build checksum-bound biological-species PGLS input from passed loss calls."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from geneloss_repro.io_utils import SchemaError, write_tsv  # noqa: E402
from geneloss_repro.pgls import (  # noqa: E402
    AGGREGATION_POLICY,
    ANALYSIS_LEVEL,
    INPUT_PASS_CHECKS,
    INPUT_PASS_SCHEMA,
    LOSS_SCOPE,
    MIN_PRIMARY_SPECIES,
    PLOIDY_PASS_CHECKS,
    PLOIDY_PASS_SCHEMA,
    PRIMARY_PREDICTOR,
    _capture_snapshot,
    _nonnegative_integer,
    _read_exact_tsv,
    _rename_directory_no_replace,
    _require_unchanged,
    _validate_species_loss_manifest,
)


MATRIX_REQUIRED = {
    "reference_gene_id",
    "biological_species",
    "aggregation_rule",
    "species_gene_status",
    "species_positive_by_rule",
    "assembly_unit_count",
    "callable_unit_count",
    "positive_unit_count",
    "confident_negative_unit_count",
    "uncertain_unit_count",
}
PLOIDY_REQUIRED = {
    "biological_species",
    "ploidy",
    "ploidy_source",
    "source_reference",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path) -> dict[str, str | int]:
    return {"basename": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}


def require_manifest_output(
    manifest: dict[str, object], directory: Path, basename: str
) -> tuple[Path, str]:
    outputs = manifest["outputs"]
    if not isinstance(outputs, list):
        raise SchemaError("species-loss manifest outputs must be a list")
    matches = [row for row in outputs if isinstance(row, dict) and row.get("basename") == basename]
    if len(matches) != 1:
        raise SchemaError(f"species-loss manifest must bind exactly one {basename}")
    row = matches[0]
    if set(row) != {"basename", "sha256"}:
        raise SchemaError(f"species-loss manifest binding for {basename} has wrong keys")
    digest = row["sha256"]
    if not isinstance(digest, str) or len(digest) != 64:
        raise SchemaError(f"species-loss manifest binding for {basename} has invalid SHA-256")
    path = directory / basename
    if not path.is_file() or sha256(path) != digest:
        raise SchemaError(f"species-loss manifest binding for {basename} does not match bytes")
    return path, digest


def canonical_int(value: str, *, column: str, line: int, path: Path) -> int:
    return _nonnegative_integer(value, column=column, line=line, path=path)


def read_ploidy(path: Path, expected_species: list[str]) -> tuple[list[dict[str, str]], dict[str, float]]:
    fields, rows = _read_exact_tsv(path)
    missing = sorted(PLOIDY_REQUIRED.difference(fields))
    if missing:
        raise SchemaError(f"{path.name}: missing ploidy columns: {', '.join(missing)}")
    by_species: dict[str, float] = {}
    normalized: list[dict[str, str]] = []
    for row, line in rows:
        species = row["biological_species"]
        if not species or species in by_species:
            raise SchemaError(f"{path.name}:{line}: empty or duplicate biological_species")
        ploidy = canonical_int(row["ploidy"], column="ploidy", line=line, path=path)
        if ploidy <= 0 or ploidy > 1024:
            raise SchemaError(f"{path.name}:{line}: ploidy must be an integer in [1,1024]")
        if not row["ploidy_source"] or not row["source_reference"]:
            raise SchemaError(f"{path.name}:{line}: ploidy provenance fields must be non-empty")
        by_species[species] = math.log2(ploidy)
        normalized.append(row)
    if list(by_species) != expected_species:
        raise SchemaError(f"{path.name}: ploidy ledger order/set differs from species-loss cohort")
    return normalized, by_species


def derive_counts(
    matrix_path: Path,
    shared_path: Path,
    expected_species: list[str],
    expected_gene_count: int,
) -> tuple[dict[str, tuple[int, int]], int]:
    fields, rows = _read_exact_tsv(matrix_path)
    missing = sorted(MATRIX_REQUIRED.difference(fields))
    if missing:
        raise SchemaError(f"{matrix_path.name}: missing matrix columns: {', '.join(missing)}")
    shared_fields, shared_rows = _read_exact_tsv(shared_path)
    if "reference_gene_id" not in shared_fields:
        raise SchemaError(f"{shared_path.name}: missing reference_gene_id")
    shared_ids: set[str] = set()
    for row, line in shared_rows:
        gene = row["reference_gene_id"]
        if not gene or gene in shared_ids:
            raise SchemaError(f"{shared_path.name}:{line}: empty or duplicate reference gene")
        shared_ids.add(gene)

    expected_pairs = expected_gene_count * len(expected_species)
    if len(rows) != expected_pairs:
        raise SchemaError(
            f"{matrix_path.name}: found {len(rows)} rows; expected {expected_pairs}"
        )
    by_gene: dict[str, dict[str, str]] = {}
    counts = {species: [0, 0] for species in expected_species}
    for row, line in rows:
        gene = row["reference_gene_id"]
        species = row["biological_species"]
        if not gene or species not in counts:
            raise SchemaError(f"{matrix_path.name}:{line}: invalid gene/species key")
        species_rows = by_gene.setdefault(gene, {})
        if species in species_rows:
            raise SchemaError(f"{matrix_path.name}:{line}: duplicate gene/species row")
        if row["aggregation_rule"] != "all_units_positive":
            raise SchemaError(f"{matrix_path.name}:{line}: PGLS forbids any-unit aggregation")
        unit_count = canonical_int(
            row["assembly_unit_count"], column="assembly_unit_count", line=line, path=matrix_path
        )
        callable_count = canonical_int(
            row["callable_unit_count"], column="callable_unit_count", line=line, path=matrix_path
        )
        positive_count = canonical_int(
            row["positive_unit_count"], column="positive_unit_count", line=line, path=matrix_path
        )
        negative_count = canonical_int(
            row["confident_negative_unit_count"],
            column="confident_negative_unit_count",
            line=line,
            path=matrix_path,
        )
        uncertain_count = canonical_int(
            row["uncertain_unit_count"], column="uncertain_unit_count", line=line, path=matrix_path
        )
        if unit_count <= 0 or callable_count > unit_count:
            raise SchemaError(f"{matrix_path.name}:{line}: invalid unit/callable counts")
        status = row["species_gene_status"]
        if status == "positive_complete":
            if not (
                callable_count == positive_count == unit_count
                and negative_count == uncertain_count == 0
                and row["species_positive_by_rule"] == "true"
            ):
                raise SchemaError(f"{matrix_path.name}:{line}: invalid positive_complete row")
        elif status == "not_positive":
            if not (
                callable_count == negative_count == unit_count
                and positive_count == uncertain_count == 0
                and row["species_positive_by_rule"] == "false"
            ):
                raise SchemaError(f"{matrix_path.name}:{line}: invalid not_positive row")
        elif status in {"positive_partial", "uncertain"}:
            if row["species_positive_by_rule"] != "false":
                raise SchemaError(f"{matrix_path.name}:{line}: partial/uncertain row is positive")
        else:
            raise SchemaError(f"{matrix_path.name}:{line}: invalid species_gene_status {status!r}")
        species_rows[species] = status
        if gene not in shared_ids and status in {"positive_complete", "not_positive"}:
            counts[species][1] += 1
            if status == "positive_complete":
                counts[species][0] += 1

    if len(by_gene) != expected_gene_count or any(
        list(species_rows) != expected_species for species_rows in by_gene.values()
    ):
        raise SchemaError(f"{matrix_path.name}: gene/species grid order or closure is invalid")
    expected_shared = {
        gene
        for gene, species_rows in by_gene.items()
        if all(species_rows[species] == "positive_complete" for species in expected_species)
    }
    if shared_ids != expected_shared:
        raise SchemaError("shared positive-complete gene set is not the exact species intersection")
    for species, (positive, denominator) in counts.items():
        if denominator <= 0 or positive > denominator:
            raise SchemaError(f"{species}: invalid derived PGLS numerator/denominator")
    return {species: (values[0], values[1]) for species, values in counts.items()}, len(shared_ids)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-loss-dir", required=True, type=Path)
    parser.add_argument("--ploidy-ledger", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output = args.output_dir
    if output.exists() or output.is_symlink():
        raise SystemExit(f"ERROR: output already exists: {output}")
    manifest_path = args.species_loss_dir / "species_loss_summary.json"
    snapshots = {
        "manifest": _capture_snapshot(manifest_path),
        "ploidy": _capture_snapshot(args.ploidy_ledger),
    }
    try:
        manifest = json.loads(snapshots["manifest"].payload.decode("utf-8"))
        species_rows = manifest.get("species_aggregation") if isinstance(manifest, dict) else None
        if not isinstance(species_rows, list):
            raise SchemaError("species-loss manifest lacks species_aggregation")
        species = [row.get("biological_species") for row in species_rows if isinstance(row, dict)]
        if len(species) < MIN_PRIMARY_SPECIES or any(not isinstance(item, str) for item in species):
            raise SchemaError(f"PGLS input requires at least {MIN_PRIMARY_SPECIES} valid species")
        _validate_species_loss_manifest(snapshots["manifest"], expected_species=species)
        matrix_path, expected_matrix_sha256 = require_manifest_output(
            manifest, args.species_loss_dir, "species_gene_matrix.tsv"
        )
        shared_path, expected_shared_sha256 = require_manifest_output(
            manifest, args.species_loss_dir, "shared_positive_complete_genes.tsv"
        )
        snapshots["matrix"] = _capture_snapshot(matrix_path)
        snapshots["shared"] = _capture_snapshot(shared_path)
        if snapshots["matrix"].sha256 != expected_matrix_sha256:
            raise SchemaError("species_gene_matrix.tsv changed after manifest validation")
        if snapshots["shared"].sha256 != expected_shared_sha256:
            raise SchemaError(
                "shared_positive_complete_genes.tsv changed after manifest validation"
            )
        with tempfile.TemporaryDirectory(prefix="species-pgls-input-snapshots.") as frozen:
            frozen_root = Path(frozen)
            frozen_ploidy = frozen_root / f"ploidy.{args.ploidy_ledger.name}"
            frozen_matrix = frozen_root / f"matrix.{matrix_path.name}"
            frozen_shared = frozen_root / f"shared.{shared_path.name}"
            frozen_ploidy.write_bytes(snapshots["ploidy"].payload)
            frozen_matrix.write_bytes(snapshots["matrix"].payload)
            frozen_shared.write_bytes(snapshots["shared"].payload)
            _, ploidy = read_ploidy(frozen_ploidy, species)
            counts, shared_count = derive_counts(
                frozen_matrix,
                frozen_shared,
                species,
                int(manifest["reference_gene_count"]),
            )
        if shared_count != manifest["shared_positive_complete_gene_count"]:
            raise SchemaError("shared count differs from species-loss manifest")

        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
        try:
            analysis_rows = [
                {
                    "biological_species": name,
                    "analysis_level": ANALYSIS_LEVEL,
                    "loss_scope": LOSS_SCOPE,
                    "lineage_specific_nonshared_positive_loss_count": counts[name][0],
                    "callable_denominator": counts[name][1],
                    PRIMARY_PREDICTOR: format(ploidy[name], ".12g"),
                }
                for name in species
            ]
            data_path = write_tsv(
                staging / "pgls_input.tsv",
                analysis_rows,
                [
                    "biological_species",
                    "analysis_level",
                    "loss_scope",
                    "lineage_specific_nonshared_positive_loss_count",
                    "callable_denominator",
                    PRIMARY_PREDICTOR,
                ],
            )
            ploidy_report_path = staging / "ploidy_ledger_pass.json"
            ploidy_report_path.write_text(
                json.dumps(
                    {
                        "schema_version": PLOIDY_PASS_SCHEMA,
                        "workflow": "species_ploidy_ledger_validation",
                        "workflow_version": "1.0.0",
                        "status": "PASS",
                        "analysis_level": ANALYSIS_LEVEL,
                        "predictor": PRIMARY_PREDICTOR,
                        "ploidy_ledger": snapshots["ploidy"].public_binding(),
                        "biological_species": species,
                        "checks": {key: True for key in sorted(PLOIDY_PASS_CHECKS)},
                    },
                    indent=2,
                    sort_keys=True,
                    allow_nan=False,
                )
                + "\n",
                encoding="utf-8",
            )
            input_report_path = staging / "pgls_input_pass.json"
            input_report_path.write_text(
                json.dumps(
                    {
                        "schema_version": INPUT_PASS_SCHEMA,
                        "workflow": "species_pgls_input_builder",
                        "workflow_version": "1.0.0",
                        "status": "PASS",
                        "analysis_level": ANALYSIS_LEVEL,
                        "loss_scope": LOSS_SCOPE,
                        "predictor": PRIMARY_PREDICTOR,
                        "input_data": binding(data_path),
                        "species_count": len(species),
                        "biological_species": species,
                        "aggregation_policy": AGGREGATION_POLICY,
                        "upstream_bindings": {
                            "species_loss_manifest": snapshots["manifest"].public_binding(),
                            "species_gene_matrix": snapshots["matrix"].public_binding(),
                            "shared_positive_complete_gene_set": snapshots[
                                "shared"
                            ].public_binding(),
                            "ploidy_ledger": snapshots["ploidy"].public_binding(),
                            "ploidy_ledger_pass_report": binding(ploidy_report_path),
                        },
                        "checks": {key: True for key in sorted(INPUT_PASS_CHECKS)},
                    },
                    indent=2,
                    sort_keys=True,
                    allow_nan=False,
                )
                + "\n",
                encoding="utf-8",
            )
            files = [data_path, ploidy_report_path, input_report_path]
            write_tsv(
                staging / "checksums.sha256.tsv",
                [
                    {"relative_path": path.name, "sha256": sha256(path)}
                    for path in sorted(files, key=lambda item: item.name)
                ],
                ["relative_path", "sha256"],
            )
            for snapshot in snapshots.values():
                _require_unchanged(snapshot)
            _rename_directory_no_replace(staging, output)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            raise
    except (OSError, ValueError, SchemaError, json.JSONDecodeError) as exc:
        raise SystemExit(f"ERROR: {exc}")
    print(f"Built validated PGLS input for {len(species)} biological species in {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
