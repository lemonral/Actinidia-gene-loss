"""Build an explicit reference-gene × target-haplotype loss master table.

The historical workflow had separate candidate, tBLASTX and retained-anchor
files. Downstream rate analyses need one complete table, but must not silently
turn untested genes into retained genes. This module makes the non-candidate
policy an explicit command-line decision and records rate eligibility per row.
"""

from __future__ import annotations

from pathlib import Path

from .io_utils import SchemaError, read_tsv, write_tsv
from .synorth import read_reference_coords


def _single_target(rows: list[dict[str, str]], label: str) -> str:
    samples = {row["target_sample"] for row in rows if row.get("target_sample", "")}
    if len(samples) != 1:
        raise SchemaError(f"{label}: expected exactly one target_sample, found {sorted(samples)}")
    return next(iter(samples))


def _metadata_ploidy(path: str | Path | None, target_sample: str) -> str:
    if path is None:
        return ""
    rows = read_tsv(path)
    if not rows:
        raise SchemaError(f"{path}: sample metadata has no rows")
    header = set(rows[0])
    sample_column = "target_haplotype" if "target_haplotype" in header else "sample_id" if "sample_id" in header else ""
    if not sample_column or "ploidy" not in header:
        raise SchemaError(f"{path}: sample metadata needs target_haplotype or sample_id, plus ploidy")
    matches = [row for row in rows if row[sample_column] == target_sample]
    if len(matches) != 1:
        raise SchemaError(f"{path}: expected one metadata row for {target_sample!r}, found {len(matches)}")
    return matches[0]["ploidy"]


def build_loss_master(
    reference_coords_path: str | Path,
    classification_path: str | Path,
    retained_anchors_path: str | Path,
    output_path: str | Path,
    noncandidate_class: str = "unassessed",
    sample_metadata_path: str | Path | None = None,
    run_id: str = "unspecified",
) -> list[dict[str, object]]:
    """Join candidate classifications and retained anchors into one full table.

    ``noncandidate_class`` is deliberately limited to ``unassessed`` (safe
    default) or ``retained_by_synorth`` (a stated denominator policy). The
    latter should only be used when the user accepts that SynOrths noncandidate
    status is sufficient retention evidence for the intended rate definition.
    """
    if noncandidate_class not in {"unassessed", "retained_by_synorth"}:
        raise SchemaError("noncandidate_class must be unassessed or retained_by_synorth")
    if not run_id.strip():
        raise SchemaError("run_id cannot be empty")

    reference_rows = read_reference_coords(reference_coords_path)
    reference_by_id = {str(row["reference_gene"]): row for row in reference_rows}
    if len(reference_by_id) != len(reference_rows):
        raise SchemaError(f"{reference_coords_path}: duplicate reference_gene IDs")

    classifications = read_tsv(classification_path, required=["target_sample", "reference_gene", "classification"])
    if not classifications:
        raise SchemaError(f"{classification_path}: no candidate classifications")
    target_sample = _single_target(classifications, str(classification_path))
    classification_by_id = {row["reference_gene"]: row for row in classifications}
    if len(classification_by_id) != len(classifications):
        raise SchemaError(f"{classification_path}: duplicate reference_gene classifications")

    retained = read_tsv(retained_anchors_path, required=["target_sample", "reference_gene"])
    if retained:
        retained_sample = _single_target(retained, str(retained_anchors_path))
        if retained_sample != target_sample:
            raise SchemaError("classification and retained-anchor tables have different target_sample values")
    retained_ids = {row["reference_gene"] for row in retained}
    if len(retained_ids) != len(retained):
        raise SchemaError(f"{retained_anchors_path}: duplicate retained reference_gene IDs")

    unknown = (set(classification_by_id) | retained_ids) - set(reference_by_id)
    if unknown:
        preview = ", ".join(sorted(unknown)[:5])
        raise SchemaError(f"Input rows contain IDs absent from reference coordinates (e.g. {preview})")
    overlap = set(classification_by_id) & retained_ids
    if overlap:
        preview = ", ".join(sorted(overlap)[:5])
        raise SchemaError(f"A reference gene is both classified candidate and retained anchor (e.g. {preview})")

    ploidy = _metadata_ploidy(sample_metadata_path, target_sample)
    output: list[dict[str, object]] = []
    for identifier, reference in reference_by_id.items():
        if identifier in classification_by_id:
            source = "tblastx_candidate"
            classification = classification_by_id[identifier]["classification"].strip()
            eligible = classification.lower() in {"pseudogenized", "deleted"}
        elif identifier in retained_ids:
            source = "synorth_anchor"
            classification = "retained"
            eligible = True
        else:
            source = "not_in_candidate_or_anchor"
            classification = noncandidate_class
            eligible = noncandidate_class == "retained_by_synorth"
        output.append({
            "reference_gene_id": identifier,
            "target_haplotype": target_sample,
            "target_sample": target_sample,
            "ploidy": ploidy,
            "classification": classification,
            "classification_source": source,
            "rate_eligible": str(eligible).lower(),
            "reference_chromosome": reference["reference_chromosome"],
            "reference_start": reference["reference_start"],
            "reference_end": reference["reference_end"],
            "reference_strand": reference["reference_strand"],
            "run_id": run_id,
        })

    output.sort(key=lambda row: (str(row["reference_chromosome"]), int(row["reference_start"]), str(row["reference_gene_id"])))
    write_tsv(output_path, output, [
        "reference_gene_id", "target_haplotype", "target_sample", "ploidy", "classification",
        "classification_source", "rate_eligible", "reference_chromosome", "reference_start",
        "reference_end", "reference_strand",
        "run_id",
    ])
    return output
