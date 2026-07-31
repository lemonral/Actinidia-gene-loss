#!/usr/bin/env python3
"""Create a checksum-bound OrthoFinder proteome directory from an audited cohort."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path


SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


class InputError(ValueError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_true(value: str, context: str) -> bool:
    normalized = value.strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise InputError(f"{context}: expected true or false, found {value!r}")


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate_audit(audit_dir: Path, manifest: Path) -> dict[str, dict[str, str]]:
    provenance_path = audit_dir / "provenance.json"
    summary_path = audit_dir / "protein_cds_pair_audit.tsv"
    if not provenance_path.is_file() or not summary_path.is_file():
        raise InputError("audit directory lacks provenance.json or protein_cds_pair_audit.tsv")
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if provenance.get("manifest_sha256") != sha256_file(manifest):
        raise InputError("audit provenance does not bind the exact cohort manifest")
    for item in provenance.get("outputs", []):
        output = audit_dir / item["basename"]
        if not output.is_file() or sha256_file(output) != item["sha256"]:
            raise InputError(f"audit output checksum mismatch: {output.name}")
    rows = read_tsv(summary_path)
    by_terminal = {row["terminal_id"]: row for row in rows}
    if len(by_terminal) != len(rows):
        raise InputError("audit summary contains duplicate terminal IDs")
    return by_terminal


def prepare(manifest: Path, data_root: Path, audit_dir: Path, output_dir: Path) -> Path:
    if output_dir.exists():
        raise InputError(f"refusing to overwrite output directory: {output_dir}")
    manifest_rows = read_tsv(manifest)
    selected = [
        row
        for row in manifest_rows
        if parse_true(row["use_for_orthofinder"], f"{row['terminal_id']}: use_for_orthofinder")
    ]
    terminal_ids = [row["terminal_id"] for row in selected]
    if len(terminal_ids) != len(set(terminal_ids)):
        raise InputError("cohort manifest contains duplicate terminal IDs")
    if any(not SAFE_ID.fullmatch(identifier) for identifier in terminal_ids):
        raise InputError("terminal IDs must match [A-Za-z0-9_.-]+")

    audit = validate_audit(audit_dir, manifest)
    if set(audit) != {row["terminal_id"] for row in manifest_rows}:
        raise InputError("audit terminal set does not equal the cohort manifest")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    try:
        proteomes = temporary / "proteomes"
        proteomes.mkdir()
        output_rows: list[dict[str, object]] = []
        for row in selected:
            terminal_id = row["terminal_id"]
            audited = audit[terminal_id]
            if audited["codon_gate"] != "PASS":
                raise InputError(f"{terminal_id}: audit codon gate is not PASS")
            source_relpath = row["protein_path"]
            if Path(source_relpath).is_absolute():
                raise InputError(f"{terminal_id}: absolute source paths are prohibited")
            source = data_root / source_relpath
            if not source.is_file():
                raise InputError(f"{terminal_id}: protein source is missing")
            source_sha = sha256_file(source)
            if source_sha != audited["protein_sha256"]:
                raise InputError(f"{terminal_id}: current protein hash differs from the audit")
            link_name = f"{terminal_id}.faa"
            os.symlink(source.resolve(), proteomes / link_name)
            output_rows.append(
                {
                    "terminal_id": terminal_id,
                    "role": row.get("role", ""),
                    "protein_source_relpath": source_relpath,
                    "protein_sha256": source_sha,
                    "protein_records": audited["protein_records"],
                    "proteome_link": f"proteomes/{link_name}",
                    "audit_codon_eligible_ids": audited["codon_eligible_ids"],
                }
            )

        input_manifest = temporary / "orthofinder_input_manifest.tsv"
        write_tsv(input_manifest, list(output_rows[0]), output_rows)
        binding = {
            "schema_version": 1,
            "terminal_count": len(output_rows),
            "cohort_manifest_basename": manifest.name,
            "cohort_manifest_sha256": sha256_file(manifest),
            "audit_provenance_sha256": sha256_file(audit_dir / "provenance.json"),
            "audit_summary_sha256": sha256_file(audit_dir / "protein_cds_pair_audit.tsv"),
            "orthofinder_input_manifest_sha256": sha256_file(input_manifest),
            "path_policy": "repository and data-root relative paths only in manifests",
            "large_file_policy": "proteomes are symlinked; exact source bytes are SHA-256 bound",
            "status": "PASS",
        }
        (temporary / "input_binding.json").write_text(
            json.dumps(binding, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.rename(temporary, output_dir)
    except Exception:
        for child in sorted(temporary.rglob("*"), reverse=True):
            if child.is_symlink() or child.is_file():
                child.unlink()
            elif child.is_dir():
                child.rmdir()
        if temporary.exists():
            temporary.rmdir()
        raise
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--audit-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        output = prepare(args.manifest, args.data_root, args.audit_dir, args.output_dir)
    except InputError as error:
        raise SystemExit(str(error)) from error
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
