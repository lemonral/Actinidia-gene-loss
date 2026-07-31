#!/usr/bin/env python3
"""Reproduce historical shared and lineage-restricted gene-loss summaries.

The input is the final legacy-compatible A06 list directory. Each included
sample must have ``<sample>_decayed_genes.txt`` and
``<sample>_deleted_genes.txt``. The former is normalized to the maintained
historical class name ``pseudogenized``. This script is retained only for
manuscript-era reproduction; it is not the callable-aware primary analysis.
Absence from both lists means only ``not_called_loss``; it is never silently
re-labelled as retained.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


LOSS_CLASSES = ("pseudogenized", "deleted")


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def read_metadata(path: Path, include_column: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"target_haplotype", "species", include_column}
    missing = required.difference(rows[0] if rows else {})
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
    selected = [row for row in rows if truthy(row[include_column])]
    if not selected:
        raise ValueError(f"{path}: no rows selected by {include_column}")
    samples = [row["target_haplotype"] for row in selected]
    if len(samples) != len(set(samples)):
        raise ValueError(f"{path}: duplicate target_haplotype values")
    return selected


def read_ids(path: Path) -> set[str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    values = {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}
    if not values:
        raise ValueError(f"{path}: no gene identifiers")
    return values


def read_loss_sets(loss_dir: Path, samples: list[str]) -> dict[str, dict[str, set[str]]]:
    result: dict[str, dict[str, set[str]]] = {}
    for sample in samples:
        pseudogenized = read_ids(loss_dir / f"{sample}_decayed_genes.txt")
        deleted = read_ids(loss_dir / f"{sample}_deleted_genes.txt")
        overlap = pseudogenized.intersection(deleted)
        if overlap:
            example = sorted(overlap)[0]
            raise ValueError(f"{sample}: gene occurs in both classes: {example}")
        result[sample] = {"pseudogenized": pseudogenized, "deleted": deleted}
    return result


def read_reference_ids(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    ids: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            first = line.rstrip("\n").split("\t", 1)[0].split(None, 1)[0]
            if first.lower() in {"transcript_id", "gene_id", "id"}:
                continue
            ids.add(first)
    if not ids:
        raise ValueError(f"{path}: no reference IDs")
    return ids


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--loss-dir", type=Path, required=True)
    parser.add_argument("--sample-metadata", type=Path, required=True)
    parser.add_argument("--include-column", default="include_manuscript")
    parser.add_argument("--reference-coords", type=Path)
    parser.add_argument("--exclude-sample", action="append", default=[])
    parser.add_argument("--cohort-name", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    metadata = read_metadata(args.sample_metadata, args.include_column)
    excluded = set(args.exclude_sample)
    unknown_exclusions = excluded.difference(row["target_haplotype"] for row in metadata)
    if unknown_exclusions:
        raise ValueError("excluded samples are not selected by metadata: " + ", ".join(sorted(unknown_exclusions)))
    metadata = [row for row in metadata if row["target_haplotype"] not in excluded]
    if not metadata:
        raise ValueError("cohort is empty after exclusions")

    samples = [row["target_haplotype"] for row in metadata]
    species_by_sample = {row["target_haplotype"]: row["species"] for row in metadata}
    samples_by_species: dict[str, set[str]] = defaultdict(set)
    for sample, species in species_by_sample.items():
        samples_by_species[species].add(sample)
    species = sorted(samples_by_species)

    loss_by_sample_class = read_loss_sets(args.loss_dir, samples)
    class_by_sample_gene: dict[str, dict[str, str]] = {}
    loss_by_sample: dict[str, set[str]] = {}
    all_loss_genes: set[str] = set()
    for sample in samples:
        classes: dict[str, str] = {}
        for loss_class in LOSS_CLASSES:
            for gene in loss_by_sample_class[sample][loss_class]:
                classes[gene] = loss_class
        class_by_sample_gene[sample] = classes
        loss_by_sample[sample] = set(classes)
        all_loss_genes.update(classes)

    reference_ids = read_reference_ids(args.reference_coords)
    if reference_ids is not None:
        outside_reference = all_loss_genes.difference(reference_ids)
        if outside_reference:
            raise ValueError(f"{len(outside_reference)} loss IDs are absent from reference coordinates")

    shared = set.intersection(*(loss_by_sample[sample] for sample in samples))
    prevalence_rows: list[dict[str, object]] = []
    lineage_rows: list[dict[str, object]] = []
    long_rows: list[dict[str, object]] = []
    histogram = Counter()

    for gene in sorted(all_loss_genes):
        affected_samples = [sample for sample in samples if gene in loss_by_sample[sample]]
        affected_species_any = sorted({species_by_sample[sample] for sample in affected_samples})
        affected_species_complete = sorted(
            sp for sp in species if samples_by_species[sp].issubset(affected_samples)
        )
        classes = Counter(class_by_sample_gene[sample][gene] for sample in affected_samples)
        histogram[len(affected_samples)] += 1
        row = {
            "reference_gene": gene,
            "affected_sample_count": len(affected_samples),
            "cohort_sample_count": len(samples),
            "sample_prevalence": len(affected_samples) / len(samples),
            "affected_species_any_count": len(affected_species_any),
            "affected_species_complete_count": len(affected_species_complete),
            "cohort_species_count": len(species),
            "shared_all_terminals": gene in shared,
            "lineage_restricted_any_species": len(affected_species_any) == 1,
            "pseudogenized_call_count": classes["pseudogenized"],
            "deleted_call_count": classes["deleted"],
            "affected_samples": ";".join(affected_samples),
            "affected_species_any": ";".join(affected_species_any),
            "affected_species_complete": ";".join(affected_species_complete),
            "cohort_name": args.cohort_name,
            "run_id": args.run_id,
        }
        prevalence_rows.append(row)
        if len(affected_species_any) == 1:
            lineage_rows.append(row)
        for sample in affected_samples:
            long_rows.append(
                {
                    "reference_gene": gene,
                    "target_haplotype": sample,
                    "species": species_by_sample[sample],
                    "classification": class_by_sample_gene[sample][gene],
                    "shared_all_terminals": gene in shared,
                    "analysis_set": "shared" if gene in shared else "non_shared",
                    "cohort_name": args.cohort_name,
                    "run_id": args.run_id,
                }
            )

    sample_rows = []
    for sample in samples:
        loss = loss_by_sample[sample]
        sample_rows.append(
            {
                "target_haplotype": sample,
                "species": species_by_sample[sample],
                "pseudogenized_count": len(loss_by_sample_class[sample]["pseudogenized"]),
                "deleted_count": len(loss_by_sample_class[sample]["deleted"]),
                "total_loss_count": len(loss),
                "shared_loss_count": len(shared),
                "non_shared_loss_count": len(loss.difference(shared)),
                "shared_fraction_of_sample_loss": len(shared) / len(loss) if loss else 0,
                "cohort_name": args.cohort_name,
                "run_id": args.run_id,
            }
        )

    species_rows = []
    for sp in species:
        sp_samples = sorted(samples_by_species[sp])
        any_loss = set.union(*(loss_by_sample[sample] for sample in sp_samples))
        complete_loss = set.intersection(*(loss_by_sample[sample] for sample in sp_samples))
        species_rows.append(
            {
                "species": sp,
                "haplotype_count": len(sp_samples),
                "haplotypes": ";".join(sp_samples),
                "any_haplotype_loss_count": len(any_loss),
                "all_haplotypes_loss_count": len(complete_loss),
                "partial_haplotype_loss_count": len(any_loss.difference(complete_loss)),
                "shared_all_terminals_count": len(shared),
                "cohort_name": args.cohort_name,
                "run_id": args.run_id,
            }
        )

    histogram_rows = [
        {
            "affected_sample_count": count,
            "cohort_sample_count": len(samples),
            "gene_count": histogram[count],
            "cohort_name": args.cohort_name,
            "run_id": args.run_id,
        }
        for count in range(1, len(samples) + 1)
    ]

    summary = {
        "analysis_role": "historical_manuscript_reproduction_only",
        "cohort_name": args.cohort_name,
        "run_id": args.run_id,
        "sample_count": len(samples),
        "species_count": len(species),
        "samples": samples,
        "excluded_samples": sorted(excluded),
        "reference_gene_count": len(reference_ids) if reference_ids is not None else None,
        "genes_called_lost_in_at_least_one_sample": len(all_loss_genes),
        "shared_loss_gene_count": len(shared),
        "shared_fraction_of_any_loss_union": len(shared) / len(all_loss_genes),
        "lineage_restricted_gene_count": len(lineage_rows),
        "source_loss_directory_basename": args.loss_dir.resolve().name,
        "decayed_semantics": "historical pseudogenized label; not used by primary callable analysis",
        "absence_semantics": "not_called_loss; absence is not assumed retained",
    }

    args.output_dir.mkdir(parents=True, exist_ok=False)
    columns = list(prevalence_rows[0])
    write_tsv(args.output_dir / "gene_loss_prevalence.tsv", prevalence_rows, columns)
    write_tsv(args.output_dir / "shared_loss_genes.tsv", [row for row in prevalence_rows if row["shared_all_terminals"]], columns)
    write_tsv(args.output_dir / "lineage_restricted_loss_genes.tsv", lineage_rows, columns)
    write_tsv(args.output_dir / "loss_calls_long.tsv", long_rows, list(long_rows[0]))
    write_tsv(args.output_dir / "sample_loss_summary.tsv", sample_rows, list(sample_rows[0]))
    write_tsv(args.output_dir / "species_loss_summary.tsv", species_rows, list(species_rows[0]))
    write_tsv(args.output_dir / "prevalence_histogram.tsv", histogram_rows, list(histogram_rows[0]))
    (args.output_dir / "shared_loss_gene_ids.txt").write_text("\n".join(sorted(shared)) + "\n", encoding="utf-8")
    (args.output_dir / "cohort_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
