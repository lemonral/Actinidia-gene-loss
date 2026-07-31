#!/usr/bin/env python3
"""Prepare exact, normalized QC inputs for the primary Actinidia cohort.

This is a read-only join over already completed basic-statistics and BUSCO
tables.  It never launches BUSCO or scans genome FASTA files.  The explicit
selection manifest chooses one row from each producer table, optionally
rekeys a legacy sample to its corrected assembly-unit identifier, and binds
the analyzed primary gene set to its standardization manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile


class PrimaryQCError(RuntimeError):
    pass


SELECTION_COLUMNS = (
    "assembly_unit_id",
    "source_sample",
    "source_accession",
    "basic_stats_table",
    "genome_busco_table",
    "protein_busco_table",
    "primary_annotation_manifest",
    "biological_species",
    "individual_id",
    "haplotype_or_subgenome",
    "accession",
    "publisher_assembly_scope",
    "local_qc_scope",
    "publisher_assembly_provenance",
    "publisher_annotation_provenance",
    "publisher_protein_provenance",
    "decision_reason",
)

METADATA_COLUMNS = (
    "assembly_unit_id", "qc_sample", "biological_species", "individual_id",
    "haplotype_or_subgenome", "accession", "decision_status",
    "publisher_assembly_scope", "local_qc_scope",
    "publisher_assembly_provenance", "publisher_annotation_provenance",
    "publisher_protein_provenance", "decision_reason",
)

BASIC_COLUMNS = (
    "sample", "current_or_alternative", "accession", "source_url",
    "genome_path", "gff_path", "protein_path", "genome_sequence_count",
    "genome_total_bp", "genome_ungapped_bp", "genome_n_bp",
    "genome_n_percent", "genome_gc_bp", "genome_gc_percent",
    "genome_longest_bp", "genome_n50_bp", "genome_l50",
    "gff_feature_rows", "gff_invalid_rows", "gff_gene_count",
    "gff_mrna_count", "gff_transcript_count",
    "gff_mrna_or_transcript_count", "gff_cds_count", "gff_exon_count",
    "protein_sequence_count", "protein_empty_sequence_count",
    "protein_total_aa", "protein_longest_aa", "protein_n50_aa",
    "protein_l50", "protein_internal_stop_record_count",
    "protein_terminal_stop_record_count", "protein_internal_stop_character_count",
    "protein_nonstandard_character_record_count",
    "protein_nonstandard_character_count",
)

BUSCO_COLUMNS = (
    "sample", "busco_version", "dataset", "dataset_creation_date", "mode",
    "input_path", "C_percent", "S_percent", "D_percent", "F_percent",
    "M_percent", "n", "C_count", "S_count", "D_count", "F_count",
    "M_count", "short_summary_path",
)

LEGACY_BASIC_COLUMNS = (
    "assembly_unit_id", "legacy_sample",
    *(column for column in BASIC_COLUMNS if column not in {"sample", "genome_path", "gff_path", "protein_path"}),
    "source_basename", "source_sha256",
)
LEGACY_BUSCO_COLUMNS = (
    "assembly_unit_id", "legacy_sample",
    *(column for column in BUSCO_COLUMNS if column not in {"sample", "input_path", "short_summary_path"}),
    "source_basename", "source_sha256",
)

ANALYSIS_SCOPE_COLUMNS = (
    "assembly_unit_id", "analysis_chromosome_count",
    "publisher_annotated_gene_count", "analysis_primary_gene_count",
    "analysis_primary_transcript_count", "invalid_coding_gene_count",
    "primary_gene_set_policy", "primary_annotation_manifest_basename",
    "primary_annotation_manifest_sha256",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path: Path) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size == 0:
        raise PrimaryQCError(f"missing or empty file: {path}")
    return {"basename": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}


def resolve(root: Path, value: str, *, allow_na: bool = False) -> Path | None:
    if allow_na and value == "NA":
        return None
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise PrimaryQCError(f"unsafe data-root-relative path: {value!r}")
    path = (root / relative).absolute()
    if not path.is_relative_to(root) or not path.is_file():
        raise PrimaryQCError(f"missing or unsafe input: {value!r}")
    return path


def read_exact(path: Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != columns:
            raise PrimaryQCError(f"{path.name}: columns differ from exact producer schema")
        rows = list(reader)
    if not rows:
        raise PrimaryQCError(f"{path.name}: empty table")
    return rows


def unique_row(path: Path, columns: tuple[str, ...], sample: str) -> dict[str, str]:
    matches = [row for row in read_exact(path, columns) if row["sample"] == sample]
    if len(matches) != 1:
        raise PrimaryQCError(f"{path.name}: expected one row for {sample!r}, found {len(matches)}")
    return dict(matches[0])


def normalized_source_row(
    path: Path, columns: tuple[str, ...], sample: str, *, table_kind: str
) -> dict[str, str]:
    """Read a current producer row or a checksum-bound legacy public row."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = tuple((csv.reader(handle, delimiter="\t").__next__()))
    if header == columns:
        return unique_row(path, columns, sample)
    legacy_columns = LEGACY_BASIC_COLUMNS if table_kind == "basic" else LEGACY_BUSCO_COLUMNS
    if header != legacy_columns:
        raise PrimaryQCError(f"{path.name}: unsupported {table_kind} source schema")
    matches = [row for row in read_exact(path, legacy_columns) if row["assembly_unit_id"] == sample]
    if len(matches) != 1:
        raise PrimaryQCError(f"{path.name}: expected one legacy row for {sample!r}, found {len(matches)}")
    legacy = matches[0]
    if table_kind == "basic":
        converted = {column: legacy.get(column, "") for column in BASIC_COLUMNS}
        converted.update({
            "sample": sample,
            "genome_path": f"legacy_exact_bound/{sample}.genome",
            "gff_path": f"legacy_exact_bound/{sample}.annotation",
            "protein_path": f"legacy_exact_bound/{sample}.protein",
        })
    else:
        converted = {column: legacy.get(column, "") for column in BUSCO_COLUMNS}
        converted.update({
            "sample": sample,
            "input_path": f"legacy_exact_bound/{sample}.{table_kind}",
            "short_summary_path": f"legacy_exact_bound/{sample}.short_summary.txt",
        })
    return converted


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def primary_scope(basic: dict[str, str], manifest_path: Path | None) -> dict[str, str]:
    if manifest_path is None:
        return {
            "analysis_chromosome_count": basic["genome_sequence_count"],
            "publisher_annotated_gene_count": basic["gff_gene_count"],
            "analysis_primary_gene_count": basic["protein_sequence_count"],
            "analysis_primary_transcript_count": basic["protein_sequence_count"],
            "invalid_coding_gene_count": "NA",
            "primary_gene_set_policy": "legacy exact-bound analyzed protein set; no new extraction",
            "primary_annotation_manifest_basename": "NA",
            "primary_annotation_manifest_sha256": "NA",
        }
    try:
        report = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise PrimaryQCError(f"cannot read primary annotation manifest {manifest_path}: {error}") from error
    counts = report.get("counts")
    policy = report.get("policy")
    if report.get("status") != "PASS" or report.get("publication_gate") != "PASS":
        raise PrimaryQCError(f"primary annotation manifest is not PASS: {manifest_path}")
    if not isinstance(counts, dict) or not isinstance(policy, dict):
        raise PrimaryQCError(f"malformed primary annotation manifest: {manifest_path}")
    required = ("chromosome_sequences", "selected_genes", "selected_transcripts", "invalid_coding_genes")
    if any(not isinstance(counts.get(key), int) for key in required):
        raise PrimaryQCError(f"primary annotation counts are incomplete: {manifest_path}")
    return {
        "analysis_chromosome_count": str(counts["chromosome_sequences"]),
        "publisher_annotated_gene_count": basic["gff_gene_count"],
        "analysis_primary_gene_count": str(counts["selected_genes"]),
        "analysis_primary_transcript_count": str(counts["selected_transcripts"]),
        "invalid_coding_gene_count": str(counts["invalid_coding_genes"]),
        "primary_gene_set_policy": str(policy.get("fallback_selection", "declared primary set")),
        "primary_annotation_manifest_basename": manifest_path.name,
        "primary_annotation_manifest_sha256": sha256(manifest_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--expected-units", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    staging: Path | None = None
    try:
        root = args.data_root.resolve()
        selections = read_exact(args.selection, SELECTION_COLUMNS)
        ids = [row["assembly_unit_id"] for row in selections]
        if len(selections) != args.expected_units or len(ids) != len(set(ids)):
            raise PrimaryQCError("selection unit count or uniqueness failed")
        output = args.output_dir.absolute()
        if output.exists():
            raise PrimaryQCError(f"refusing to overwrite output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))

        metadata_rows: list[dict[str, str]] = []
        basic_rows: list[dict[str, str]] = []
        genome_rows: list[dict[str, str]] = []
        protein_rows: list[dict[str, str]] = []
        scope_rows: list[dict[str, str]] = []
        input_bindings: dict[str, dict[str, object]] = {}
        for selected in selections:
            unit = selected["assembly_unit_id"]
            source_sample = selected["source_sample"]
            basic_path = resolve(root, selected["basic_stats_table"])
            genome_path = resolve(root, selected["genome_busco_table"])
            protein_path = resolve(root, selected["protein_busco_table"])
            manifest_path = resolve(root, selected["primary_annotation_manifest"], allow_na=True)
            assert basic_path and genome_path and protein_path
            for path in (basic_path, genome_path, protein_path):
                input_bindings.setdefault(str(path.relative_to(root)), binding(path))
            if manifest_path:
                input_bindings.setdefault(str(manifest_path.relative_to(root)), binding(manifest_path))
            basic = normalized_source_row(basic_path, BASIC_COLUMNS, source_sample, table_kind="basic")
            genome = normalized_source_row(genome_path, BUSCO_COLUMNS, source_sample, table_kind="genome_busco")
            protein = normalized_source_row(protein_path, BUSCO_COLUMNS, source_sample, table_kind="protein_busco")
            if genome["mode"] != "euk_genome_min" or protein["mode"] != "proteins":
                raise PrimaryQCError(f"{unit}: unexpected BUSCO mode")
            if basic["accession"] != selected["source_accession"]:
                raise PrimaryQCError(f"{unit}: source accession differs from the selected producer row")
            # Metadata may correct a previously misidentified legacy label.
            # Preserve the old value in the checksum-bound selection row but
            # normalize the publication join to the accepted accession.
            basic["accession"] = selected["accession"]
            for row in (basic, genome, protein):
                row["sample"] = unit
            metadata_rows.append({
                "assembly_unit_id": unit, "qc_sample": unit,
                "biological_species": selected["biological_species"],
                "individual_id": selected["individual_id"],
                "haplotype_or_subgenome": selected["haplotype_or_subgenome"],
                "accession": selected["accession"], "decision_status": "current",
                "publisher_assembly_scope": selected["publisher_assembly_scope"],
                "local_qc_scope": selected["local_qc_scope"],
                "publisher_assembly_provenance": selected["publisher_assembly_provenance"],
                "publisher_annotation_provenance": selected["publisher_annotation_provenance"],
                "publisher_protein_provenance": selected["publisher_protein_provenance"],
                "decision_reason": selected["decision_reason"],
            })
            basic_rows.append(basic)
            genome_rows.append(genome)
            protein_rows.append(protein)
            scope_rows.append({"assembly_unit_id": unit, **primary_scope(basic, manifest_path)})

        order = {unit: index for index, unit in enumerate(ids)}
        metadata_rows.sort(key=lambda row: order[row["assembly_unit_id"]])
        basic_rows.sort(key=lambda row: order[row["sample"]])
        genome_rows.sort(key=lambda row: order[row["sample"]])
        protein_rows.sort(key=lambda row: order[row["sample"]])
        scope_rows.sort(key=lambda row: order[row["assembly_unit_id"]])
        files = {
            "assembly_decisions.tsv": (METADATA_COLUMNS, metadata_rows),
            "basic_stats.tsv": (BASIC_COLUMNS, basic_rows),
            "genome_busco.tsv": (BUSCO_COLUMNS, genome_rows),
            "protein_busco.tsv": (BUSCO_COLUMNS, protein_rows),
            "analysis_scope_primary_annotation.tsv": (ANALYSIS_SCOPE_COLUMNS, scope_rows),
        }
        for name, (columns, rows) in files.items():
            write_tsv(staging / name, columns, rows)
        report = {
            "schema_version": 1, "workflow": "primary_actinidia_qc_input_preparation",
            "status": "PASS", "assembly_unit_count": len(ids), "assembly_units": ids,
            "selection": binding(args.selection), "source_inputs": input_bindings,
            "outputs": {name: binding(staging / name) for name in files},
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        checksum_names = [*files, "run_manifest.json"]
        with (staging / "checksums.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("file", "bytes", "sha256"))
            for name in sorted(checksum_names):
                item = binding(staging / name)
                writer.writerow((name, item["bytes"], item["sha256"]))
        os.replace(staging, output)
        staging = None
        return 0
    except (OSError, ValueError, PrimaryQCError) as error:
        print(f"ERROR: {error}")
        return 2
    finally:
        if staging is not None:
            shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
