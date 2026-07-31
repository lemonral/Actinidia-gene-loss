#!/usr/bin/env python3
"""Build one fail-closed, consumer-compatible analysis-unit manifest.

The repository deliberately separates biological assembly declarations from
download instructions.  This program joins ``config/assemblies.tsv`` to
``config/downloads.tsv`` only after every enabled download has been reconciled
with an exact downloader JSON record and the bytes currently stored below
``--data-root``.  The resulting TSV is the sole manifest that should be passed
to assembly-QC, BUSCO, or operational gene-loss consumers.

Rows are emitted for assembly units selected by at least one of
``include_qc``, ``include_gene_loss``, or ``include_species_tree``.  By default
each selected unit must have one enabled genome, GFF, and protein asset.  CDS
is retained when available but is not required unless requested explicitly.

Canonical columns such as ``assembly_unit_id`` and ``biological_species`` are
written alongside compatibility aliases used by the current scripts:
``sample``, ``species``, ``target_haplotype``, ``current_or_alternative``,
``include_downstream``, ``genome``, ``gff``, ``protein``, and ``source_url``.
The aliases do not create additional biological replicates: ``sample`` and
``target_haplotype`` both identify the declared assembly unit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Iterable
from urllib.parse import urlparse


ASSEMBLY_REQUIRED_COLUMNS = {
    "assembly_unit_id",
    "biological_species",
    "individual_id",
    "haplotype_or_subgenome",
    "ploidy",
    "accession",
    "version",
    "source_bundle_id",
    "assembly_scope",
    "partition_rule",
    "genome_url",
    "annotation_url",
    "protein_url",
    "expected_genome_sha256",
    "expected_annotation_sha256",
    "expected_protein_sha256",
    "qc_status",
    "include_qc",
    "include_gene_loss",
    "include_species_tree",
    "exclusion_reason",
    "notes",
}
DOWNLOAD_REQUIRED_COLUMNS = {
    "asset_id",
    "assembly_unit_id",
    "asset_type",
    "url",
    "relative_path",
    "expected_bytes",
    "md5",
    "sha256",
    "download",
    "source_note",
}
INCLUDE_COLUMNS = ("include_qc", "include_gene_loss", "include_species_tree")
SUPPORTED_ROLES = ("genome", "gff", "protein", "cds")
DEFAULT_REQUIRED_ROLES = ("genome", "gff", "protein")
SAFE_UNIT_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
HEX_RE = re.compile(r"^[0-9a-fA-F]+$")

CANONICAL_COLUMNS = (
    "assembly_unit_id",
    "biological_species",
    "individual_id",
    "haplotype_or_subgenome",
    "ploidy",
    "accession",
    "version",
    "source_bundle_id",
    "assembly_scope",
    "partition_rule",
    "qc_status",
    "include_qc",
    "include_gene_loss",
    "include_species_tree",
    "exclusion_reason",
    "notes",
)
COMPATIBILITY_COLUMNS = (
    "sample",
    "species",
    "individual",
    "haplotype",
    "target_haplotype",
    "current_or_alternative",
    "include_downstream",
    "source_url",
    "legacy_loss_fasta_stem",
)
ASSET_PATH_COLUMNS = SUPPORTED_ROLES
ASSET_METADATA_SUFFIXES = (
    "asset_id",
    "source_url",
    "relative_path",
    "expected_bytes",
    "publisher_md5",
    "publisher_sha256",
    "local_sha256",
    "integrity_level",
)
OUTPUT_COLUMNS = (
    *CANONICAL_COLUMNS,
    *COMPATIBILITY_COLUMNS,
    *ASSET_PATH_COLUMNS,
    *(f"{role}_{suffix}" for role in SUPPORTED_ROLES for suffix in ASSET_METADATA_SUFFIXES),
)


class ResolutionError(RuntimeError):
    """Raised when declarations, downloads, reports, or local bytes disagree."""


def parse_boolean(value: str, *, path: Path, line_number: int, field: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"true", "false"}:
        raise ResolutionError(
            f"{path}:{line_number}: {field} must be exactly true or false"
        )
    return normalized == "true"


def validate_hex(value: str, length: int, *, context: str) -> str:
    normalized = value.strip().lower()
    if normalized and (len(normalized) != length or not HEX_RE.fullmatch(normalized)):
        raise ResolutionError(f"{context}: expected a {length}-character hexadecimal value")
    return normalized


def read_tsv(path: Path, required_columns: set[str], label: str) -> list[tuple[int, dict[str, str]]]:
    path = path.expanduser().resolve()
    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except OSError as error:
        raise ResolutionError(f"Cannot open {label} {path}: {error}") from error
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ResolutionError(f"{path}: missing a TSV header")
        reader.fieldnames = [field.strip() for field in reader.fieldnames]
        missing = required_columns.difference(reader.fieldnames)
        if missing:
            raise ResolutionError(
                f"{path}: missing {label} columns: {', '.join(sorted(missing))}"
            )
        rows: list[tuple[int, dict[str, str]]] = []
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise ResolutionError(f"{path}:{line_number}: more fields than the TSV header")
            row = {key: (value or "").strip() for key, value in raw.items()}
            if not any(row.values()):
                continue
            rows.append((line_number, row))
    if not rows:
        raise ResolutionError(f"{path}: {label} has no data rows")
    return rows


def load_assemblies(path: Path) -> list[dict[str, str]]:
    resolved = path.expanduser().resolve()
    rows = read_tsv(resolved, ASSEMBLY_REQUIRED_COLUMNS, "assembly manifest")
    assemblies: list[dict[str, str]] = []
    seen: set[str] = set()
    for line_number, row in rows:
        unit = row["assembly_unit_id"]
        if not SAFE_UNIT_ID.fullmatch(unit):
            raise ResolutionError(
                f"{resolved}:{line_number}: unsafe assembly_unit_id {unit!r}; use 1-128 "
                "ASCII letters, digits, dots, underscores, or hyphens"
            )
        if unit in seen:
            raise ResolutionError(f"{resolved}:{line_number}: duplicate assembly_unit_id {unit!r}")
        seen.add(unit)
        for field in (
            "biological_species",
            "individual_id",
            "haplotype_or_subgenome",
            "ploidy",
            "accession",
            "version",
            "source_bundle_id",
            "assembly_scope",
            "partition_rule",
            "qc_status",
        ):
            if not row[field]:
                raise ResolutionError(f"{resolved}:{line_number}: {unit} has an empty {field}")
        for field in INCLUDE_COLUMNS:
            row[field] = "true" if parse_boolean(
                row[field], path=resolved, line_number=line_number, field=field
            ) else "false"
        for field in (
            "expected_genome_sha256",
            "expected_annotation_sha256",
            "expected_protein_sha256",
        ):
            row[field] = validate_hex(
                row[field], 64, context=f"{resolved}:{line_number}: {field}"
            )
        assemblies.append(row)
    return assemblies


def is_under(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def resolve_data_path(data_root: Path, relative_path: str, *, context: str) -> Path:
    raw = Path(relative_path)
    if not relative_path.strip() or raw.is_absolute():
        raise ResolutionError(f"{context}: relative_path must be non-empty and relative")
    resolved = (data_root / raw).resolve()
    if not is_under(resolved, data_root):
        raise ResolutionError(f"{context}: relative_path escapes --data-root: {relative_path!r}")
    return resolved


def load_downloads(path: Path, data_root: Path) -> list[dict[str, object]]:
    resolved = path.expanduser().resolve()
    rows = read_tsv(resolved, DOWNLOAD_REQUIRED_COLUMNS, "download manifest")
    downloads: list[dict[str, object]] = []
    asset_ids: set[str] = set()
    relative_paths: set[str] = set()
    enabled_roles: set[tuple[str, str]] = set()
    for line_number, text_row in rows:
        asset_id = text_row["asset_id"]
        unit = text_row["assembly_unit_id"]
        role = text_row["asset_type"].lower()
        if not asset_id:
            raise ResolutionError(f"{resolved}:{line_number}: empty asset_id")
        if asset_id in asset_ids:
            raise ResolutionError(f"{resolved}:{line_number}: duplicate asset_id {asset_id!r}")
        asset_ids.add(asset_id)
        if not unit:
            raise ResolutionError(f"{resolved}:{line_number}: {asset_id} has no assembly_unit_id")
        if role not in SUPPORTED_ROLES:
            raise ResolutionError(
                f"{resolved}:{line_number}: {asset_id} has unsupported asset_type {role!r}"
            )
        parsed_url = urlparse(text_row["url"])
        if parsed_url.scheme not in {"http", "https"} or not parsed_url.hostname:
            raise ResolutionError(
                f"{resolved}:{line_number}: {asset_id} has an unsupported URL {text_row['url']!r}"
            )
        relative = text_row["relative_path"]
        if relative in relative_paths:
            raise ResolutionError(
                f"{resolved}:{line_number}: duplicate relative_path {relative!r}"
            )
        relative_paths.add(relative)
        local_path = resolve_data_path(
            data_root, relative, context=f"{resolved}:{line_number}: {asset_id}"
        )
        enabled = parse_boolean(
            text_row["download"], path=resolved, line_number=line_number, field="download"
        )
        if enabled and (unit, role) in enabled_roles:
            raise ResolutionError(
                f"{resolved}:{line_number}: duplicate enabled {role} role for {unit}"
            )
        if enabled:
            enabled_roles.add((unit, role))
        expected_bytes: int | None = None
        if text_row["expected_bytes"]:
            try:
                expected_bytes = int(text_row["expected_bytes"])
            except ValueError as error:
                raise ResolutionError(
                    f"{resolved}:{line_number}: {asset_id} expected_bytes is not an integer"
                ) from error
            if expected_bytes <= 0:
                raise ResolutionError(
                    f"{resolved}:{line_number}: {asset_id} expected_bytes must be positive"
                )
        md5 = validate_hex(text_row["md5"], 32, context=f"{resolved}:{line_number}: {asset_id} md5")
        sha256 = validate_hex(
            text_row["sha256"], 64, context=f"{resolved}:{line_number}: {asset_id} sha256"
        )
        if enabled and expected_bytes is None and not (md5 or sha256):
            raise ResolutionError(
                f"{resolved}:{line_number}: enabled asset {asset_id} has neither an expected "
                "size nor a publisher checksum"
            )
        row: dict[str, object] = dict(text_row)
        row.update(
            {
                "asset_type": role,
                "download": enabled,
                "expected_bytes": expected_bytes,
                "md5": md5,
                "sha256": sha256,
                "local_path": local_path,
                "manifest_line": line_number,
            }
        )
        downloads.append(row)
    return downloads


def load_report(path: Path) -> dict[str, dict[str, object]]:
    resolved = path.expanduser().resolve()
    try:
        payload = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ResolutionError(f"Cannot read downloader report {resolved}: {error}") from error
    if not isinstance(payload, list):
        raise ResolutionError(f"{resolved}: downloader report must contain a JSON list")
    reports: dict[str, dict[str, object]] = {}
    for index, item in enumerate(payload, start=1):
        context = f"{resolved}: record {index}"
        if not isinstance(item, dict):
            raise ResolutionError(f"{context} is not a JSON object")
        asset_id = item.get("asset_id")
        if not isinstance(asset_id, str) or not asset_id:
            error_text = item.get("error")
            if item.get("status") == "failed" and isinstance(error_text, str):
                raise ResolutionError(f"{context}: downloader recorded failure: {error_text}")
            raise ResolutionError(f"{context}: missing asset_id")
        if asset_id in reports:
            raise ResolutionError(f"{context}: duplicate report asset_id {asset_id!r}")
        reports[asset_id] = item
    return reports


def file_digests(path: Path, need_md5: bool) -> tuple[str, str]:
    sha = hashlib.sha256()
    md5 = hashlib.md5() if need_md5 else None
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            sha.update(block)
            if md5 is not None:
                md5.update(block)
    return sha.hexdigest(), md5.hexdigest() if md5 is not None else ""


def validate_enabled_assets(
    downloads: Iterable[dict[str, object]],
    report_rows: dict[str, dict[str, object]],
) -> dict[tuple[str, str], dict[str, object]]:
    enabled = [row for row in downloads if row["download"] is True]
    enabled_ids = {str(row["asset_id"]) for row in enabled}
    report_ids = set(report_rows)
    if report_ids != enabled_ids:
        missing = sorted(enabled_ids - report_ids)
        extra = sorted(report_ids - enabled_ids)
        details = []
        if missing:
            details.append("missing report rows: " + ", ".join(missing))
        if extra:
            details.append("unexpected report rows: " + ", ".join(extra))
        raise ResolutionError("Downloader report does not exactly cover enabled assets; " + "; ".join(details))

    by_role: dict[tuple[str, str], dict[str, object]] = {}
    for row in enabled:
        asset_id = str(row["asset_id"])
        unit = str(row["assembly_unit_id"])
        role = str(row["asset_type"])
        report = report_rows[asset_id]
        context = f"download report asset {asset_id}"
        if report.get("status") != "verified":
            raise ResolutionError(f"{context}: status is {report.get('status')!r}, not 'verified'")
        problems = report.get("problems")
        if not isinstance(problems, list) or problems:
            raise ResolutionError(f"{context}: problems must be an empty JSON list")
        for field, expected in (
            ("assembly_unit_id", unit),
            ("asset_type", role),
            ("relative_path", str(row["relative_path"])),
        ):
            if report.get(field) != expected:
                raise ResolutionError(
                    f"{context}: {field} {report.get(field)!r} does not match {expected!r}"
                )
        local_path = row["local_path"]
        assert isinstance(local_path, Path)
        if not local_path.is_file():
            raise ResolutionError(f"{asset_id}: local asset is missing or not a file: {local_path}")
        observed_bytes = local_path.stat().st_size
        if isinstance(report.get("bytes"), bool) or not isinstance(report.get("bytes"), int):
            raise ResolutionError(f"{context}: bytes must be a JSON integer")
        if report["bytes"] != observed_bytes:
            raise ResolutionError(
                f"{context}: recorded bytes {report['bytes']} do not match local bytes {observed_bytes}"
            )
        expected_bytes = row["expected_bytes"]
        if expected_bytes is not None and observed_bytes != expected_bytes:
            raise ResolutionError(
                f"{asset_id}: local bytes {observed_bytes} do not match expected_bytes {expected_bytes}"
            )
        report_sha = report.get("sha256")
        if not isinstance(report_sha, str):
            raise ResolutionError(f"{context}: sha256 must be a string")
        report_sha = validate_hex(report_sha, 64, context=f"{context}: sha256")
        local_sha, local_md5 = file_digests(local_path, bool(row["md5"]))
        if local_sha != report_sha:
            raise ResolutionError(f"{asset_id}: local SHA-256 does not match downloader report")
        if row["sha256"] and local_sha != row["sha256"]:
            raise ResolutionError(f"{asset_id}: local SHA-256 does not match publisher SHA-256")
        report_md5 = report.get("md5", "")
        if not isinstance(report_md5, str):
            raise ResolutionError(f"{context}: md5 must be a string")
        if row["md5"]:
            normalized_report_md5 = validate_hex(report_md5, 32, context=f"{context}: md5")
            if local_md5 != row["md5"] or normalized_report_md5 != row["md5"]:
                raise ResolutionError(f"{asset_id}: local/report MD5 does not match publisher MD5")
        declared = bool(row["md5"] or row["sha256"])
        if report.get("publisher_checksum_declared") is not declared:
            raise ResolutionError(
                f"{context}: publisher_checksum_declared disagrees with download manifest"
            )
        integrity_level = (
            "publisher_sha256" if row["sha256"]
            else "publisher_md5" if row["md5"]
            else "expected_size_plus_local_sha256"
        )
        validated = dict(row)
        validated["local_sha256"] = local_sha
        validated["integrity_level"] = integrity_level
        key = (unit, role)
        if key in by_role:
            raise ResolutionError(f"duplicate enabled {role} role for {unit}")
        by_role[key] = validated
    return by_role


def validate_assembly_asset_contract(assembly: dict[str, str], role: str, asset: dict[str, object]) -> None:
    unit = assembly["assembly_unit_id"]
    declared_url_field = {
        "genome": "genome_url",
        "gff": "annotation_url",
        "protein": "protein_url",
        "cds": "",
    }[role]
    if declared_url_field:
        declared_url = assembly[declared_url_field]
        if declared_url and declared_url != asset["url"]:
            raise ResolutionError(
                f"{unit}: {declared_url_field} disagrees with enabled {role} download URL"
            )
    declared_hash_field = {
        "genome": "expected_genome_sha256",
        "gff": "expected_annotation_sha256",
        "protein": "expected_protein_sha256",
        "cds": "",
    }[role]
    if declared_hash_field:
        declared_hash = assembly[declared_hash_field]
        if declared_hash and declared_hash != asset["local_sha256"]:
            raise ResolutionError(
                f"{unit}: {declared_hash_field} disagrees with the resolved local {role} SHA-256"
            )


def build_output_rows(
    assemblies: list[dict[str, str]],
    assets_by_role: dict[tuple[str, str], dict[str, object]],
    required_roles: Iterable[str],
) -> list[dict[str, str]]:
    required = tuple(dict.fromkeys(required_roles))
    if not required:
        raise ResolutionError("At least one --required-role is required")
    unknown = sorted(set(required).difference(SUPPORTED_ROLES))
    if unknown:
        raise ResolutionError("Unsupported required roles: " + ", ".join(unknown))
    assembly_ids = {row["assembly_unit_id"] for row in assemblies}
    unknown_asset_units = sorted({unit for unit, _role in assets_by_role}.difference(assembly_ids))
    if unknown_asset_units:
        raise ResolutionError(
            "Enabled downloads refer to unknown assembly units: " + ", ".join(unknown_asset_units)
        )

    selected = [
        row for row in assemblies
        if any(row[column] == "true" for column in INCLUDE_COLUMNS)
    ]
    if not selected:
        raise ResolutionError("No assembly units are selected by an include_* column")

    output_rows: list[dict[str, str]] = []
    for assembly in selected:
        unit = assembly["assembly_unit_id"]
        missing = [role for role in required if (unit, role) not in assets_by_role]
        if missing:
            raise ResolutionError(
                f"{unit}: missing enabled required asset roles: {', '.join(missing)}"
            )
        role_assets: dict[str, dict[str, object] | None] = {
            role: assets_by_role.get((unit, role)) for role in SUPPORTED_ROLES
        }
        for role, asset in role_assets.items():
            if asset is not None:
                validate_assembly_asset_contract(assembly, role, asset)

        row = {column: assembly[column] for column in CANONICAL_COLUMNS}
        row.update(
            {
                "sample": unit,
                "species": assembly["biological_species"],
                "individual": assembly["individual_id"],
                "haplotype": assembly["haplotype_or_subgenome"],
                "target_haplotype": unit,
                "current_or_alternative": (
                    "current" if assembly["include_gene_loss"] == "true" else "alternative"
                ),
                "include_downstream": assembly["include_gene_loss"],
                "source_url": (
                    str(role_assets["genome"]["url"])
                    if role_assets["genome"] is not None else ""
                ),
                "legacy_loss_fasta_stem": "",
            }
        )
        for role in SUPPORTED_ROLES:
            asset = role_assets[role]
            row[role] = str(asset["local_path"]) if asset is not None else ""
            values = {
                "asset_id": str(asset["asset_id"]) if asset is not None else "",
                "source_url": str(asset["url"]) if asset is not None else "",
                "relative_path": str(asset["relative_path"]) if asset is not None else "",
                "expected_bytes": (
                    str(asset["expected_bytes"])
                    if asset is not None and asset["expected_bytes"] is not None else ""
                ),
                "publisher_md5": str(asset["md5"]) if asset is not None else "",
                "publisher_sha256": str(asset["sha256"]) if asset is not None else "",
                "local_sha256": str(asset["local_sha256"]) if asset is not None else "",
                "integrity_level": str(asset["integrity_level"]) if asset is not None else "",
            }
            for suffix, value in values.items():
                row[f"{role}_{suffix}"] = value
        output_rows.append(row)
    return output_rows


def write_output(path: Path, rows: list[dict[str, str]]) -> None:
    resolved = path.expanduser().resolve()
    if resolved.exists() and not resolved.is_file():
        raise ResolutionError(f"Output path exists and is not a regular file: {resolved}")
    resolved.parent.mkdir(parents=True, exist_ok=True)
    temporary = resolved.with_name(resolved.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(resolved)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def resolve_manifests(
    *,
    assemblies_path: Path,
    downloads_path: Path,
    data_root: Path,
    report_path: Path,
    output_path: Path,
    required_roles: Iterable[str] = DEFAULT_REQUIRED_ROLES,
) -> list[dict[str, str]]:
    root = data_root.expanduser().resolve()
    if not root.is_dir():
        raise ResolutionError(f"--data-root is missing or not a directory: {root}")
    assemblies = load_assemblies(assemblies_path)
    downloads = load_downloads(downloads_path, root)
    report = load_report(report_path)
    validated_assets = validate_enabled_assets(downloads, report)
    rows = build_output_rows(assemblies, validated_assets, required_roles)
    write_output(output_path, rows)
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assemblies", required=True, type=Path)
    parser.add_argument("--downloads", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--download-report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--required-role",
        action="append",
        choices=SUPPORTED_ROLES,
        help="Required enabled role for every selected unit; repeat as needed (default: genome,gff,protein)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        rows = resolve_manifests(
            assemblies_path=args.assemblies,
            downloads_path=args.downloads,
            data_root=args.data_root,
            report_path=args.download_report,
            output_path=args.output,
            required_roles=args.required_role or DEFAULT_REQUIRED_ROLES,
        )
    except (OSError, UnicodeError, csv.Error, ResolutionError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(f"Resolved {len(rows)} analysis units into {args.output.expanduser().resolve()}")


if __name__ == "__main__":
    main()
