"""Chromosome and within-chromosome spatial summaries for pseudogene fragments.

The manuscript's spatial result uses pseudogenized/degraded genes, not all
putative candidates.  This implementation starts with the final classification
table and uses each selected tBLASTX hit's target coordinates explicitly.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import statistics
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterable, Mapping, Sequence

from .gff import parse_attributes, read_gene_catalog
from .io_utils import SchemaError, bh_adjust, format_number, natural_key, parse_int, read_tsv, write_tsv


def _gammaincc(a: float, x: float) -> float:
    """Regularized upper incomplete gamma Q(a, x), dependency-free.

    Numerical Recipes' series/continued-fraction implementation is sufficient
    for chi-square p-values in the reporting tables.  It avoids making simple
    table regeneration depend on a particular SciPy installation.
    """
    if a <= 0 or x < 0:
        raise ValueError("a must be >0 and x must be >=0")
    if x == 0:
        return 1.0
    gln = math.lgamma(a)
    eps = 3e-14
    maximum = 10000
    tiny = 1e-300
    if x < a + 1:
        ap = a
        total = 1.0 / a
        delta = total
        for _ in range(maximum):
            ap += 1.0
            delta *= x / ap
            total += delta
            if abs(delta) < abs(total) * eps:
                p = total * math.exp(-x + a * math.log(x) - gln)
                return max(0.0, min(1.0, 1.0 - p))
        raise RuntimeError("gamma series failed to converge")
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / max(b, tiny)
    h = d
    for iteration in range(1, maximum + 1):
        an = -iteration * (iteration - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            return max(0.0, min(1.0, math.exp(-x + a * math.log(x) - gln) * h))
    raise RuntimeError("gamma continued fraction failed to converge")


def chi_square_gof(observed: list[int], opportunities: list[int]) -> dict[str, object]:
    """Goodness-of-fit against gene-density opportunity, with transparent QC."""
    if len(observed) != len(opportunities) or len(observed) < 2:
        return {"statistic": "", "df": "", "p_value": "", "status": "insufficient_categories"}
    total_observed = sum(observed)
    total_opportunity = sum(opportunities)
    if total_observed == 0:
        return {"statistic": "", "df": "", "p_value": "", "status": "no_pseudogenized_fragments"}
    if total_opportunity <= 0:
        return {"statistic": "", "df": "", "p_value": "", "status": "zero_gene_opportunity"}
    expected = [total_observed * opportunity / total_opportunity for opportunity in opportunities]
    usable = [(obs, exp) for obs, exp in zip(observed, expected) if exp > 0]
    if len(usable) < 2:
        return {"statistic": "", "df": "", "p_value": "", "status": "insufficient_nonzero_expected_categories"}
    statistic = sum((obs - exp) ** 2 / exp for obs, exp in usable)
    df = len(usable) - 1
    p_value = _gammaincc(df / 2, statistic / 2)
    status = "warning_expected_lt_5" if any(exp < 5 for _, exp in usable) else "ok"
    return {"statistic": statistic, "df": df, "p_value": p_value, "status": status}


def _read_chromosome_lengths(path: str | Path | None, catalog: list[dict[str, object]]) -> tuple[dict[str, int], str]:
    if path is None:
        lengths: dict[str, int] = defaultdict(int)
        for row in catalog:
            chrom = str(row["target_chromosome"])
            lengths[chrom] = max(lengths[chrom], int(row["target_end"]))
        return dict(lengths), "derived_from_max_GFF_gene_end"
    rows = read_tsv(path)
    lengths = {}
    for row in rows:
        chromosome = row.get("chromosome") or row.get("Chromosome") or row.get("sequence_id") or ""
        length = row.get("length") or row.get("Length") or ""
        if not chromosome or not length:
            raise SchemaError(f"{path}: each row needs chromosome/Chromosome and length/Length")
        lengths[chromosome] = parse_int(length, "length", str(path))
    return lengths, "explicit_chromosome_lengths_table"


def _canonical_hit_rows(classification_path: str | Path, loss_class: str) -> tuple[list[dict[str, object]], str]:
    rows = read_tsv(
        classification_path,
        required=["target_sample", "reference_gene", "classification", "best_subject_id", "best_subject_start", "best_subject_end"],
    )
    samples = {row["target_sample"] for row in rows}
    if len(samples) != 1:
        raise SchemaError(f"{classification_path}: spatial summary expects one target sample per table")
    kept: list[dict[str, object]] = []
    seen: set[str] = set()
    for row in rows:
        if row["classification"] != loss_class:
            continue
        gene = row["reference_gene"]
        if gene in seen:
            raise SchemaError(f"{classification_path}: duplicate classified reference_gene {gene!r}")
        seen.add(gene)
        if not row["best_subject_id"] or not row["best_subject_start"] or not row["best_subject_end"]:
            raise SchemaError(
                f"{classification_path}: {gene} is {loss_class} but lacks selected target hit coordinates"
            )
        start = parse_int(row["best_subject_start"], "best_subject_start", str(classification_path))
        end = parse_int(row["best_subject_end"], "best_subject_end", str(classification_path))
        kept.append({
            "target_sample": row["target_sample"], "reference_gene": gene,
            "target_chromosome": row["best_subject_id"], "target_start": min(start, end), "target_end": max(start, end),
            "target_midpoint": (min(start, end) + max(start, end)) / 2,
            "classification": loss_class,
        })
    return kept, next(iter(samples))


def _equal_width_bin(midpoint: float, chromosome_length: int, number_of_bins: int) -> int:
    # 1-based coordinates and bins.  Every coordinate 1..length maps exactly once.
    return min(number_of_bins, max(1, int((midpoint - 1) * number_of_bins // chromosome_length) + 1))


def _bin_boundaries(chromosome_length: int, number_of_bins: int) -> list[tuple[int, int, int]]:
    return [
        (index, math.floor((index - 1) * chromosome_length / number_of_bins) + 1,
         math.floor(index * chromosome_length / number_of_bins))
        for index in range(1, number_of_bins + 1)
    ]


def _legacy_nested_bins(chromosome_length: int, number_of_bins: int) -> list[tuple[int, int, int]]:
    """Reproduce the legacy nested-midpoint interval construction.

    These are overlapping intervals, not bins in the usual statistical sense.
    The mode is retained only to compare old figures; do not use it for new
    inferential tests.
    """
    midpoint = chromosome_length // 2
    width = midpoint // number_of_bins
    if width < 1:
        raise SchemaError(f"chromosome length {chromosome_length} too short for {number_of_bins} legacy bins")
    return [
        (index, max(0, midpoint - index * width), min(chromosome_length, midpoint + index * width))
        for index in range(1, number_of_bins + 1)
    ]


def spatial_summary(
    classification_path: str | Path,
    target_gff: str | Path,
    output_dir: str | Path,
    sample_id: str | None = None,
    chromosome_lengths_path: str | Path | None = None,
    number_of_bins: int = 5,
    bin_mode: str = "equal-width",
    gene_feature: str = "gene",
    loss_class: str = "pseudogenized",
) -> dict[str, Path]:
    """Write auditable chromosome/bin summaries and opportunity-based tests."""
    if number_of_bins < 2:
        raise SchemaError("number_of_bins must be >=2")
    if bin_mode not in {"equal-width", "legacy-nested-midpoint"}:
        raise SchemaError("bin_mode must be equal-width or legacy-nested-midpoint")
    hits, table_sample = _canonical_hit_rows(classification_path, loss_class)
    if sample_id is not None and sample_id != table_sample:
        raise SchemaError(f"sample_id={sample_id!r} does not match classification target_sample={table_sample!r}")
    sample = table_sample
    catalog = read_gene_catalog(target_gff, gene_feature=gene_feature)
    chromosome_lengths, length_source = _read_chromosome_lengths(chromosome_lengths_path, catalog)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    catalog_by_chromosome: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in catalog:
        catalog_by_chromosome[str(row["target_chromosome"])].append(row)
    hit_by_chromosome: dict[str, list[dict[str, object]]] = defaultdict(list)
    unmapped: list[dict[str, object]] = []
    for hit in hits:
        chrom = str(hit["target_chromosome"])
        if chrom not in catalog_by_chromosome:
            unmapped.append({**hit, "unmapped_reason": "target_chromosome_not_in_GFF_gene_catalog"})
        else:
            hit_by_chromosome[chrom].append(hit)

    chromosomes = sorted(catalog_by_chromosome, key=natural_key)
    inter_rows: list[dict[str, object]] = []
    for chromosome in chromosomes:
        total_genes = len(catalog_by_chromosome[chromosome])
        lost = len(hit_by_chromosome[chromosome])
        inter_rows.append({
            "sample_id": sample, "target_chromosome": chromosome, "total_target_genes": total_genes,
            "pseudogenized_reference_genes": lost,
            "loss_fragment_per_target_gene": lost / total_genes if total_genes else "",
            "loss_class": loss_class, "coordinate_length_source": length_source,
        })
    inter_test = chi_square_gof(
        [int(row["pseudogenized_reference_genes"]) for row in inter_rows],
        [int(row["total_target_genes"]) for row in inter_rows],
    )
    inter_test_rows = [{
        "sample_id": sample, "scope": "all_target_chromosomes", "loss_class": loss_class,
        "n_categories": len(inter_rows), "n_pseudogenized_fragments": sum(int(row["pseudogenized_reference_genes"]) for row in inter_rows),
        "chi_square": format_number(inter_test["statistic"]) if inter_test["statistic"] != "" else "",
        "df": inter_test["df"], "p_value": format_number(inter_test["p_value"]) if inter_test["p_value"] != "" else "",
        "test_status": inter_test["status"],
        "null_model": "fragment_locations_proportional_to_target_gene_count",
    }]

    intra_rows: list[dict[str, object]] = []
    intra_tests: list[dict[str, object]] = []
    for chromosome in chromosomes:
        length = chromosome_lengths.get(chromosome)
        if length is None or length <= 0:
            raise SchemaError(f"no positive chromosome length available for {chromosome!r}")
        bins = _bin_boundaries(length, number_of_bins) if bin_mode == "equal-width" else _legacy_nested_bins(length, number_of_bins)
        gene_counts = {number: 0 for number, _, _ in bins}
        fragment_counts = {number: 0 for number, _, _ in bins}
        for gene in catalog_by_chromosome[chromosome]:
            midpoint = float(gene["target_midpoint"])
            if bin_mode == "equal-width":
                gene_counts[_equal_width_bin(midpoint, length, number_of_bins)] += 1
            else:
                for number, start, end in bins:
                    if int(gene["target_start"]) >= start and int(gene["target_end"]) <= end:
                        gene_counts[number] += 1
        for hit in hit_by_chromosome[chromosome]:
            midpoint = float(hit["target_midpoint"])
            if bin_mode == "equal-width":
                fragment_counts[_equal_width_bin(midpoint, length, number_of_bins)] += 1
            else:
                for number, start, end in bins:
                    if int(hit["target_start"]) >= start and int(hit["target_end"]) <= end:
                        fragment_counts[number] += 1
        chrom_rows: list[dict[str, object]] = []
        for number, start, end in bins:
            total_genes = gene_counts[number]
            lost = fragment_counts[number]
            chrom_rows.append({
                "sample_id": sample, "target_chromosome": chromosome, "bin": number,
                "bin_start": start, "bin_end": end, "chromosome_length": length,
                "total_target_genes": total_genes, "pseudogenized_reference_genes": lost,
                "loss_fragment_per_target_gene": lost / total_genes if total_genes else "",
                "loss_class": loss_class, "bin_mode": bin_mode,
            })
        intra_rows.extend(chrom_rows)
        if bin_mode == "legacy-nested-midpoint":
            result = {"statistic": "", "df": "", "p_value": "", "status": "not_tested_overlapping_legacy_intervals"}
        else:
            result = chi_square_gof(
                [int(row["pseudogenized_reference_genes"]) for row in chrom_rows],
                [int(row["total_target_genes"]) for row in chrom_rows],
            )
        intra_tests.append({
            "sample_id": sample, "target_chromosome": chromosome, "loss_class": loss_class, "bin_mode": bin_mode,
            "n_bins": len(chrom_rows), "n_pseudogenized_fragments": sum(int(row["pseudogenized_reference_genes"]) for row in chrom_rows),
            "chi_square": result["statistic"], "df": result["df"], "p_value": result["p_value"], "test_status": result["status"],
            "null_model": "fragment_locations_proportional_to_target_gene_count_per_bin",
        })
    adjusted = bh_adjust([float(row["p_value"]) if row["p_value"] not in ("", None) else None for row in intra_tests])
    for row, p_adjusted in zip(intra_tests, adjusted):
        row["p_value_bh"] = p_adjusted if p_adjusted is not None else ""

    catalog_path = output / f"{sample}.target_gene_catalog.tsv"
    hits_path = output / f"{sample}.pseudogene_fragment_locations.tsv"
    unmapped_path = output / f"{sample}.unmapped_pseudogene_fragments.tsv"
    inter_path = output / f"{sample}.inter_chromosome.tsv"
    inter_test_path = output / f"{sample}.inter_chromosome_chi_square.tsv"
    intra_path = output / f"{sample}.intra_chromosome.tsv"
    intra_test_path = output / f"{sample}.intra_chromosome_chi_square.tsv"
    write_tsv(catalog_path, catalog, ["target_gene", "target_chromosome", "target_start", "target_end", "target_midpoint"])
    write_tsv(hits_path, hits, ["target_sample", "reference_gene", "target_chromosome", "target_start", "target_end", "target_midpoint", "classification"])
    write_tsv(unmapped_path, unmapped, ["target_sample", "reference_gene", "target_chromosome", "target_start", "target_end", "target_midpoint", "classification", "unmapped_reason"])
    write_tsv(inter_path, inter_rows, ["sample_id", "target_chromosome", "total_target_genes", "pseudogenized_reference_genes", "loss_fragment_per_target_gene", "loss_class", "coordinate_length_source"])
    write_tsv(inter_test_path, inter_test_rows, ["sample_id", "scope", "loss_class", "n_categories", "n_pseudogenized_fragments", "chi_square", "df", "p_value", "test_status", "null_model"])
    write_tsv(intra_path, intra_rows, ["sample_id", "target_chromosome", "bin", "bin_start", "bin_end", "chromosome_length", "total_target_genes", "pseudogenized_reference_genes", "loss_fragment_per_target_gene", "loss_class", "bin_mode"])
    write_tsv(intra_test_path, intra_tests, ["sample_id", "target_chromosome", "loss_class", "bin_mode", "n_bins", "n_pseudogenized_fragments", "chi_square", "df", "p_value", "test_status", "null_model", "p_value_bh"])
    return {
        "catalog": catalog_path, "fragments": hits_path, "unmapped": unmapped_path,
        "inter": inter_path, "inter_test": inter_test_path, "intra": intra_path, "intra_test": intra_test_path,
    }


# ---------------------------------------------------------------------------
# Manifest-driven multi-assembly position analysis
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _SpatialAssembly:
    assembly_unit_id: str
    biological_species: str
    haplotype_or_subgenome: str
    assembly_scope: str
    genome: Path
    gff: Path
    genome_declared_sha256: str
    gff_declared_sha256: str


@dataclass(frozen=True)
class _PositiveCall:
    assembly_unit_id: str
    reference_gene_id: str
    classification: str


@dataclass(frozen=True)
class _FeatureCoordinate:
    assembly_unit_id: str
    reference_gene_id: str
    classification: str
    chromosome: str
    start: int
    end: int

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2


@dataclass(frozen=True)
class _CentromereInterval:
    assembly_unit_id: str
    chromosome: str
    start: int
    end: int
    evidence_source: str


def _open_text_auto(path: str | Path) -> IO[str]:
    """Open plain or gzip-compressed text based on the final filename suffix."""
    source = Path(path)
    if source.suffix.lower() == ".gz":
        return gzip.open(source, "rt", encoding="utf-8-sig", newline="")
    return source.open("r", encoding="utf-8-sig", newline="")


def _read_tsv_auto(path: str | Path) -> tuple[list[dict[str, str]], list[str]]:
    """Read a strict headered TSV, including ``.tsv.gz`` inputs."""
    source = Path(path)
    try:
        handle = _open_text_auto(source)
    except OSError as exc:
        raise SchemaError(f"cannot open {source}: {exc}") from exc
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SchemaError(f"{source}: missing TSV header")
        fields = [field.strip() for field in reader.fieldnames]
        if not fields or any(not field for field in fields):
            raise SchemaError(f"{source}: TSV header contains an empty column name")
        if len(fields) != len(set(fields)):
            raise SchemaError(f"{source}: TSV header contains duplicate column names")
        reader.fieldnames = fields
        rows: list[dict[str, str]] = []
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise SchemaError(f"{source}:{line_number}: more fields than the TSV header")
            row = {str(key).strip(): (value or "").strip() for key, value in raw.items()}
            if any(row.values()):
                row["__line_number"] = str(line_number)
                rows.append(row)
    if not rows:
        raise SchemaError(f"{source}: no data rows")
    return rows, fields


def _require_columns(path: str | Path, fields: Sequence[str], required: Sequence[str]) -> None:
    missing = [column for column in required if column not in fields]
    if missing:
        raise SchemaError(
            f"{path}: missing required column(s): {', '.join(missing)}; "
            f"found: {', '.join(fields)}"
        )


def _strict_boolean(value: str, *, path: str | Path, line_number: str, field: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"true", "false"}:
        raise SchemaError(
            f"{path}:{line_number}: {field} must be exactly true or false, found {value!r}"
        )
    return normalized == "true"


def _sha256_small_input(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: str | Path, payload: Mapping[str, object]) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        delete=False,
        dir=output.parent,
        prefix=f".{output.name}.",
        suffix=".tmp",
    ) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = handle.name
    os.replace(temporary, output)
    return output


def _resolve_asset_path(value: str, manifest_path: Path, field: str, unit: str) -> Path:
    if not value:
        raise SchemaError(f"{manifest_path}: {unit}: empty {field} path")
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest_path.parent / path
    path = path.resolve()
    if not path.is_file():
        raise SchemaError(f"{manifest_path}: {unit}: {field} file does not exist: {path}")
    return path


def _load_spatial_assemblies(
    manifest_path: str | Path,
    expected_units: set[str],
    *,
    unit_column: str,
    species_column: str,
    haplotype_column: str,
    scope_column: str,
    genome_column: str,
    gff_column: str,
    include_column: str | None,
) -> dict[str, _SpatialAssembly]:
    source = Path(manifest_path).expanduser().resolve()
    rows, fields = _read_tsv_auto(source)
    required = [
        unit_column,
        species_column,
        haplotype_column,
        scope_column,
        genome_column,
        gff_column,
    ]
    if include_column:
        required.append(include_column)
    _require_columns(source, fields, required)
    by_unit: dict[str, dict[str, str]] = {}
    selected_by_manifest: set[str] = set()
    for row in rows:
        unit = row[unit_column]
        line_number = row["__line_number"]
        if not unit:
            raise SchemaError(f"{source}:{line_number}: empty {unit_column}")
        if unit in by_unit:
            raise SchemaError(f"{source}:{line_number}: duplicate assembly unit {unit!r}")
        by_unit[unit] = row
        if include_column and _strict_boolean(
            row[include_column], path=source, line_number=line_number, field=include_column
        ):
            selected_by_manifest.add(unit)

    missing = expected_units.difference(by_unit)
    if missing:
        raise SchemaError(
            f"{source}: positive-call units absent from assembly manifest: "
            f"{', '.join(sorted(missing, key=natural_key))}"
        )
    if include_column and selected_by_manifest != expected_units:
        absent_from_calls = selected_by_manifest.difference(expected_units)
        not_selected = expected_units.difference(selected_by_manifest)
        details: list[str] = []
        if absent_from_calls:
            details.append(
                "manifest-selected but absent from positive-call scope="
                + ",".join(sorted(absent_from_calls, key=natural_key))
            )
        if not_selected:
            details.append(
                "positive-call units not selected by manifest="
                + ",".join(sorted(not_selected, key=natural_key))
            )
        raise SchemaError(f"{source}: assembly-unit scope mismatch; {'; '.join(details)}")

    assemblies: dict[str, _SpatialAssembly] = {}
    for unit in sorted(expected_units, key=natural_key):
        row = by_unit[unit]
        for column in (species_column, haplotype_column, scope_column):
            if not row[column]:
                raise SchemaError(
                    f"{source}:{row['__line_number']}: {unit}: empty required metadata {column}"
                )
        assemblies[unit] = _SpatialAssembly(
            assembly_unit_id=unit,
            biological_species=row[species_column],
            haplotype_or_subgenome=row[haplotype_column],
            assembly_scope=row[scope_column],
            genome=_resolve_asset_path(row[genome_column], source, genome_column, unit),
            gff=_resolve_asset_path(row[gff_column], source, gff_column, unit),
            genome_declared_sha256=(
                row.get("genome_local_sha256", "")
                or row.get("expected_genome_sha256", "")
            ),
            gff_declared_sha256=(
                row.get("gff_local_sha256", "")
                or row.get("annotation_local_sha256", "")
                or row.get("expected_annotation_sha256", "")
            ),
        )
    for field in ("genome", "gff"):
        owners: dict[Path, list[str]] = defaultdict(list)
        for unit, assembly in assemblies.items():
            owners[getattr(assembly, field)].append(unit)
        reused = {path: units for path, units in owners.items() if len(units) > 1}
        if reused:
            example_path, units = next(iter(reused.items()))
            raise SchemaError(
                f"{source}: assembly units {', '.join(sorted(units, key=natural_key))} share the "
                f"same {field} file {example_path}; split combined polyploid bundles into one "
                "sequence/GFF partition per assembly unit before position analysis"
            )
    return assemblies


def _load_positive_calls(
    path: str | Path,
    positive_classes: set[str],
    *,
    unit_column: str,
    gene_column: str,
    classification_column: str,
) -> tuple[dict[tuple[str, str], _PositiveCall], set[str]]:
    rows, fields = _read_tsv_auto(path)
    _require_columns(path, fields, [unit_column, gene_column, classification_column])
    all_keys: set[tuple[str, str]] = set()
    call_units: set[str] = set()
    positive: dict[tuple[str, str], _PositiveCall] = {}
    for row in rows:
        line_number = row["__line_number"]
        unit = row[unit_column]
        gene = row[gene_column]
        classification = row[classification_column].lower()
        if not unit or not gene or not classification:
            raise SchemaError(
                f"{path}:{line_number}: unit, reference gene, and classification must be nonempty"
            )
        key = (unit, gene)
        if key in all_keys:
            raise SchemaError(
                f"{path}:{line_number}: duplicate final call for assembly unit/gene {unit!r}/{gene!r}"
            )
        all_keys.add(key)
        call_units.add(unit)
        if classification in positive_classes:
            positive[key] = _PositiveCall(unit, gene, classification)
    if not positive:
        raise SchemaError(
            f"{path}: no rows match positive classification(s): {', '.join(sorted(positive_classes))}"
        )
    return positive, call_units


def _load_feature_coordinates(
    path: str | Path,
    positive_calls: Mapping[tuple[str, str], _PositiveCall],
    positive_classes: set[str],
    *,
    unit_column: str,
    gene_column: str,
    chromosome_column: str,
    start_column: str,
    end_column: str,
    classification_column: str | None,
) -> dict[tuple[str, str], _FeatureCoordinate]:
    rows, fields = _read_tsv_auto(path)
    required = [unit_column, gene_column, chromosome_column, start_column, end_column]
    if classification_column:
        required.append(classification_column)
    _require_columns(path, fields, required)
    coordinates: dict[tuple[str, str], _FeatureCoordinate] = {}
    for row in rows:
        line_number = row["__line_number"]
        if classification_column:
            classification = row[classification_column].lower()
            if not classification:
                raise SchemaError(f"{path}:{line_number}: empty {classification_column}")
            if classification not in positive_classes:
                continue
        else:
            classification = ""
        unit = row[unit_column]
        gene = row[gene_column]
        chromosome = row[chromosome_column]
        if not unit or not gene or not chromosome:
            raise SchemaError(
                f"{path}:{line_number}: unit, reference gene, and chromosome must be nonempty"
            )
        start = parse_int(row[start_column], start_column, f"{path}:{line_number}")
        end = parse_int(row[end_column], end_column, f"{path}:{line_number}")
        start, end = min(start, end), max(start, end)
        if start < 1:
            raise SchemaError(f"{path}:{line_number}: feature interval begins before coordinate 1")
        key = (unit, gene)
        if key in coordinates:
            raise SchemaError(
                f"{path}:{line_number}: multiple selected feature coordinates for {unit!r}/{gene!r}; "
                "supply exactly one upstream-selected hit"
            )
        call = positive_calls.get(key)
        if call is not None and classification_column and classification != call.classification:
            raise SchemaError(
                f"{path}:{line_number}: classification {classification!r} disagrees with positive-call "
                f"classification {call.classification!r} for {unit!r}/{gene!r}"
            )
        coordinates[key] = _FeatureCoordinate(
            unit,
            gene,
            call.classification if call is not None else classification,
            chromosome,
            start,
            end,
        )

    coordinate_keys = set(coordinates)
    call_keys = set(positive_calls)
    if coordinate_keys != call_keys:
        missing = call_keys.difference(coordinate_keys)
        extra = coordinate_keys.difference(call_keys)
        details: list[str] = []
        if missing:
            examples = ", ".join(f"{unit}/{gene}" for unit, gene in sorted(missing)[:5])
            details.append(f"{len(missing)} positive calls lack coordinates ({examples})")
        if extra:
            examples = ", ".join(f"{unit}/{gene}" for unit, gene in sorted(extra)[:5])
            details.append(f"{len(extra)} coordinate rows are outside positive-call scope ({examples})")
        raise SchemaError(f"{path}: coordinate/call scope mismatch: {'; '.join(details)}")
    return coordinates


def _read_fasta_lengths(path: str | Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    identifier: str | None = None
    current_length = 0
    try:
        handle = _open_text_auto(path)
    except OSError as exc:
        raise SchemaError(f"cannot open genome FASTA {path}: {exc}") from exc
    with handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    if current_length < 1:
                        raise SchemaError(f"{path}: FASTA record {identifier!r} has no sequence")
                    lengths[identifier] = current_length
                identifier = line[1:].split()[0]
                if not identifier:
                    raise SchemaError(f"{path}:{line_number}: FASTA header has no sequence ID")
                if identifier in lengths:
                    raise SchemaError(f"{path}:{line_number}: duplicate FASTA sequence ID {identifier!r}")
                current_length = 0
            elif identifier is None:
                raise SchemaError(f"{path}:{line_number}: sequence occurs before the first FASTA header")
            else:
                current_length += len("".join(line.split()))
    if identifier is not None:
        if current_length < 1:
            raise SchemaError(f"{path}: FASTA record {identifier!r} has no sequence")
        if identifier in lengths:
            raise SchemaError(f"{path}: duplicate FASTA sequence ID {identifier!r}")
        lengths[identifier] = current_length
    if not lengths:
        raise SchemaError(f"{path}: no FASTA records")
    return lengths


def _read_gene_features(path: str | Path, gene_feature: str) -> list[dict[str, object]]:
    genes: list[dict[str, object]] = []
    seen: set[str] = set()
    try:
        handle = _open_text_auto(path)
    except OSError as exc:
        raise SchemaError(f"cannot open GFF {path}: {exc}") from exc
    with handle:
        for line_number, raw in enumerate(handle, start=1):
            if raw.startswith("##FASTA"):
                break
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise SchemaError(f"{path}:{line_number}: expected 9 GFF columns, found {len(fields)}")
            if fields[2] != gene_feature:
                continue
            start = parse_int(fields[3], "start", f"{path}:{line_number}")
            end = parse_int(fields[4], "end", f"{path}:{line_number}")
            if start < 1 or end < start:
                raise SchemaError(f"{path}:{line_number}: invalid gene interval {start}-{end}")
            gene_id = parse_attributes(fields[8]).get("ID", "")
            if not gene_id:
                continue
            if gene_id in seen:
                raise SchemaError(f"{path}:{line_number}: duplicate {gene_feature} ID {gene_id!r}")
            seen.add(gene_id)
            genes.append(
                {
                    "gene_id": gene_id,
                    "chromosome": fields[0],
                    "start": start,
                    "end": end,
                    "midpoint": (start + end) / 2,
                }
            )
    if not genes:
        raise SchemaError(f"{path}: no {gene_feature!r} features with an ID attribute")
    return genes


def _load_centromeres(path: str | Path | None) -> dict[tuple[str, str], _CentromereInterval]:
    if path is None:
        return {}
    rows, fields = _read_tsv_auto(path)
    required = [
        "assembly_unit_id",
        "chromosome",
        "centromere_start",
        "centromere_end",
        "evidence_source",
    ]
    _require_columns(path, fields, required)
    intervals: dict[tuple[str, str], _CentromereInterval] = {}
    for row in rows:
        line_number = row["__line_number"]
        unit = row["assembly_unit_id"]
        chromosome = row["chromosome"]
        evidence = row["evidence_source"]
        if not unit or not chromosome or not evidence:
            raise SchemaError(
                f"{path}:{line_number}: assembly_unit_id, chromosome, and evidence_source must be nonempty"
            )
        start = parse_int(row["centromere_start"], "centromere_start", f"{path}:{line_number}")
        end = parse_int(row["centromere_end"], "centromere_end", f"{path}:{line_number}")
        start, end = min(start, end), max(start, end)
        if start < 1:
            raise SchemaError(f"{path}:{line_number}: centromere interval begins before coordinate 1")
        key = (unit, chromosome)
        if key in intervals:
            raise SchemaError(f"{path}:{line_number}: duplicate centromere interval for {unit}/{chromosome}")
        intervals[key] = _CentromereInterval(unit, chromosome, start, end, evidence)
    return intervals


def _normalized_end_distance(midpoint: float, chromosome_length: int) -> tuple[float, float, str]:
    """Return midpoint distance (bp), 0..1 half-chromosome distance, and nearest end."""
    if chromosome_length < 2:
        raise SchemaError("chromosome/sequence length must be at least 2 for normalized end distance")
    left = midpoint - 1
    right = chromosome_length - midpoint
    distance = min(left, right)
    nearest = "left" if left < right else "right" if right < left else "equidistant"
    half_span = (chromosome_length - 1) / 2
    normalized = min(1.0, max(0.0, distance / half_span))
    return distance, normalized, nearest


def _unit_interval_bin(value: float, number_of_bins: int) -> int:
    return min(number_of_bins, max(1, math.floor(value * number_of_bins) + 1))


def _distance_to_interval(midpoint: float, start: int, end: int) -> float:
    if midpoint < start:
        return start - midpoint
    if midpoint > end:
        return midpoint - end
    return 0.0


def _mean(values: Sequence[float]) -> float | str:
    return statistics.fmean(values) if values else ""


def _median(values: Sequence[float]) -> float | str:
    return statistics.median(values) if values else ""


def analyze_loss_positions(
    positive_calls_path: str | Path,
    feature_coordinates_path: str | Path,
    assembly_manifest_path: str | Path,
    output_dir: str | Path,
    *,
    analysis_label: str,
    positive_classes: Iterable[str] = ("pseudogenized",),
    number_of_bins: int = 5,
    centromeres_path: str | Path | None = None,
    require_complete_centromeres: bool = False,
    legacy_reproduction: bool = False,
    gene_feature: str = "gene",
    call_unit_column: str = "target_haplotype",
    call_gene_column: str = "reference_gene_id",
    call_classification_column: str = "classification",
    coordinate_unit_column: str = "target_haplotype",
    coordinate_gene_column: str = "reference_gene_id",
    coordinate_chromosome_column: str = "target_chromosome",
    coordinate_start_column: str = "target_start",
    coordinate_end_column: str = "target_end",
    coordinate_classification_column: str | None = "classification",
    manifest_unit_column: str = "assembly_unit_id",
    manifest_species_column: str = "biological_species",
    manifest_haplotype_column: str = "haplotype_or_subgenome",
    manifest_scope_column: str = "assembly_scope",
    manifest_genome_column: str = "genome",
    manifest_gff_column: str = "gff",
    manifest_include_column: str | None = None,
) -> dict[str, Path]:
    """Analyze accepted target-fragment positions across explicit assembly units.

    The primary bins are mutually exclusive.  The optional legacy output is an
    exact, clearly labelled reproduction of the manuscript-era cumulative
    midpoint intervals; because those intervals overlap, it is never used for
    an inferential test.  Positive call IDs, feature-coordinate IDs, assembly
    units, FASTA sequence IDs, and GFF sequence IDs are reconciled before any
    output is written.
    """
    analysis_label = analysis_label.strip()
    if not analysis_label:
        raise SchemaError("analysis_label must be nonempty")
    if number_of_bins < 2:
        raise SchemaError("number_of_bins must be >=2")
    class_values = (
        positive_classes.split(",")
        if isinstance(positive_classes, str)
        else positive_classes
    )
    classes = {value.strip().lower() for value in class_values if value.strip()}
    if not classes:
        raise SchemaError("positive_classes must contain at least one nonempty class")
    if not gene_feature.strip():
        raise SchemaError("gene_feature must be nonempty")

    positive_calls, expected_units = _load_positive_calls(
        positive_calls_path,
        classes,
        unit_column=call_unit_column,
        gene_column=call_gene_column,
        classification_column=call_classification_column,
    )
    coordinates = _load_feature_coordinates(
        feature_coordinates_path,
        positive_calls,
        classes,
        unit_column=coordinate_unit_column,
        gene_column=coordinate_gene_column,
        chromosome_column=coordinate_chromosome_column,
        start_column=coordinate_start_column,
        end_column=coordinate_end_column,
        classification_column=coordinate_classification_column,
    )
    assemblies = _load_spatial_assemblies(
        assembly_manifest_path,
        expected_units,
        unit_column=manifest_unit_column,
        species_column=manifest_species_column,
        haplotype_column=manifest_haplotype_column,
        scope_column=manifest_scope_column,
        genome_column=manifest_genome_column,
        gff_column=manifest_gff_column,
        include_column=manifest_include_column,
    )
    centromeres = _load_centromeres(centromeres_path)
    unknown_centromere_units = {unit for unit, _ in centromeres}.difference(expected_units)
    if unknown_centromere_units:
        raise SchemaError(
            f"{centromeres_path}: centromere assembly units outside analysis scope: "
            f"{', '.join(sorted(unknown_centromere_units, key=natural_key))}"
        )

    positions_rows: list[dict[str, object]] = []
    chromosome_rows: list[dict[str, object]] = []
    equal_bin_rows: list[dict[str, object]] = []
    end_bin_rows: list[dict[str, object]] = []
    unit_rows: list[dict[str, object]] = []
    legacy_rows: list[dict[str, object]] = []
    metadata_units: list[dict[str, object]] = []

    for unit in sorted(assemblies, key=natural_key):
        assembly = assemblies[unit]
        lengths = _read_fasta_lengths(assembly.genome)
        genes = _read_gene_features(assembly.gff, gene_feature)
        genes_by_chromosome: dict[str, list[dict[str, object]]] = defaultdict(list)
        for gene in genes:
            chromosome = str(gene["chromosome"])
            if chromosome not in lengths:
                raise SchemaError(
                    f"{assembly.gff}: GFF sequence {chromosome!r} is absent from genome FASTA "
                    f"for assembly unit {unit}"
                )
            if int(gene["end"]) > lengths[chromosome]:
                raise SchemaError(
                    f"{assembly.gff}: gene {gene['gene_id']!r} ends beyond FASTA sequence "
                    f"{chromosome!r} ({gene['end']} > {lengths[chromosome]})"
                )
            genes_by_chromosome[chromosome].append(gene)

        unit_coordinates = [
            coordinate for key, coordinate in coordinates.items() if key[0] == unit
        ]
        coordinates_by_chromosome: dict[str, list[_FeatureCoordinate]] = defaultdict(list)
        for coordinate in unit_coordinates:
            if coordinate.chromosome not in lengths:
                raise SchemaError(
                    f"{feature_coordinates_path}: {unit}/{coordinate.reference_gene_id}: "
                    f"sequence {coordinate.chromosome!r} is absent from its genome FASTA"
                )
            if coordinate.chromosome not in genes_by_chromosome:
                raise SchemaError(
                    f"{feature_coordinates_path}: {unit}/{coordinate.reference_gene_id}: "
                    f"sequence {coordinate.chromosome!r} has no {gene_feature!r} denominator in its GFF"
                )
            if coordinate.end > lengths[coordinate.chromosome]:
                raise SchemaError(
                    f"{feature_coordinates_path}: {unit}/{coordinate.reference_gene_id}: feature end "
                    f"{coordinate.end} exceeds {coordinate.chromosome} length {lengths[coordinate.chromosome]}"
                )
            coordinates_by_chromosome[coordinate.chromosome].append(coordinate)

        analyzed_chromosomes = set(genes_by_chromosome)
        unit_centromere_keys = {key for key in centromeres if key[0] == unit}
        extra_centromere_chromosomes = {
            chromosome for _, chromosome in unit_centromere_keys if chromosome not in analyzed_chromosomes
        }
        if extra_centromere_chromosomes:
            raise SchemaError(
                f"{centromeres_path}: {unit}: centromere sequence(s) outside the GFF gene-bearing "
                f"analysis scope: {', '.join(sorted(extra_centromere_chromosomes, key=natural_key))}"
            )
        if require_complete_centromeres and centromeres_path is None:
            raise SchemaError("require_complete_centromeres is true but no centromere table was supplied")
        if require_complete_centromeres:
            missing_centromeres = analyzed_chromosomes.difference(
                chromosome for _, chromosome in unit_centromere_keys
            )
            if missing_centromeres:
                raise SchemaError(
                    f"{centromeres_path}: {unit}: missing centromere intervals for "
                    f"{', '.join(sorted(missing_centromeres, key=natural_key))}"
                )

        end_gene_counts = {index: 0 for index in range(1, number_of_bins + 1)}
        end_loss_counts = {index: 0 for index in range(1, number_of_bins + 1)}
        unit_end_distances: list[float] = []
        unit_normalized_end_distances: list[float] = []
        supplied_centromere_chromosomes = 0

        for chromosome in sorted(analyzed_chromosomes, key=natural_key):
            length = lengths[chromosome]
            if length < number_of_bins:
                raise SchemaError(
                    f"{assembly.genome}: sequence {chromosome!r} length {length} is shorter than "
                    f"number_of_bins={number_of_bins}"
                )
            gene_list = genes_by_chromosome[chromosome]
            loss_list = coordinates_by_chromosome.get(chromosome, [])
            centromere = centromeres.get((unit, chromosome))
            if centromere is not None:
                if centromere.end > length:
                    raise SchemaError(
                        f"{centromeres_path}: {unit}/{chromosome}: centromere end {centromere.end} "
                        f"exceeds chromosome length {length}"
                    )
                supplied_centromere_chromosomes += 1

            linear_gene_counts = {index: 0 for index in range(1, number_of_bins + 1)}
            linear_loss_counts = {index: 0 for index in range(1, number_of_bins + 1)}
            for gene in gene_list:
                midpoint = float(gene["midpoint"])
                linear_gene_counts[_equal_width_bin(midpoint, length, number_of_bins)] += 1
                _, normalized, _ = _normalized_end_distance(midpoint, length)
                end_gene_counts[_unit_interval_bin(normalized, number_of_bins)] += 1

            chromosome_end_distances: list[float] = []
            chromosome_normalized_end_distances: list[float] = []
            chromosome_centromere_distances: list[float] = []
            for coordinate in sorted(
                loss_list,
                key=lambda item: (item.start, item.end, item.reference_gene_id),
            ):
                midpoint = coordinate.midpoint
                linear_bin = _equal_width_bin(midpoint, length, number_of_bins)
                linear_loss_counts[linear_bin] += 1
                end_distance, normalized_end_distance, nearest_end = _normalized_end_distance(
                    midpoint, length
                )
                end_distance_bin = _unit_interval_bin(normalized_end_distance, number_of_bins)
                end_loss_counts[end_distance_bin] += 1
                chromosome_end_distances.append(end_distance)
                chromosome_normalized_end_distances.append(normalized_end_distance)
                unit_end_distances.append(end_distance)
                unit_normalized_end_distances.append(normalized_end_distance)
                if centromere is None:
                    centromere_distance: float | str = ""
                    centromere_distance_fraction: float | str = ""
                    centromere_start: int | str = ""
                    centromere_end: int | str = ""
                    centromere_evidence = ""
                    centromere_status = "not_supplied_for_chromosome"
                else:
                    centromere_distance = _distance_to_interval(
                        midpoint, centromere.start, centromere.end
                    )
                    centromere_distance_fraction = centromere_distance / length
                    chromosome_centromere_distances.append(centromere_distance)
                    centromere_start = centromere.start
                    centromere_end = centromere.end
                    centromere_evidence = centromere.evidence_source
                    centromere_status = "independently_supplied_interval"
                positions_rows.append(
                    {
                        "analysis_label": analysis_label,
                        "assembly_unit_id": unit,
                        "biological_species": assembly.biological_species,
                        "haplotype_or_subgenome": assembly.haplotype_or_subgenome,
                        "assembly_scope": assembly.assembly_scope,
                        "reference_gene_id": coordinate.reference_gene_id,
                        "classification": coordinate.classification,
                        "chromosome": chromosome,
                        "feature_start": coordinate.start,
                        "feature_end": coordinate.end,
                        "feature_midpoint": format_number(midpoint),
                        "chromosome_length": length,
                        "equal_width_bin": linear_bin,
                        "midpoint_distance_to_nearest_end_bp": format_number(end_distance),
                        "normalized_end_distance_0_end_1_center": format_number(normalized_end_distance),
                        "end_distance_bin_1_end_to_n_center": end_distance_bin,
                        "nearest_chromosome_end": nearest_end,
                        "centromere_start": centromere_start,
                        "centromere_end": centromere_end,
                        "midpoint_distance_to_centromere_bp": format_number(centromere_distance)
                        if centromere_distance != ""
                        else "",
                        "centromere_distance_fraction_of_chromosome": format_number(
                            centromere_distance_fraction
                        )
                        if centromere_distance_fraction != ""
                        else "",
                        "centromere_status": centromere_status,
                        "centromere_evidence_source": centromere_evidence,
                    }
                )

            for bin_number, bin_start, bin_end in _bin_boundaries(length, number_of_bins):
                denominator = linear_gene_counts[bin_number]
                count = linear_loss_counts[bin_number]
                equal_bin_rows.append(
                    {
                        "analysis_label": analysis_label,
                        "analysis_mode": "primary_mutually_exclusive_equal_width",
                        "assembly_unit_id": unit,
                        "biological_species": assembly.biological_species,
                        "haplotype_or_subgenome": assembly.haplotype_or_subgenome,
                        "chromosome": chromosome,
                        "chromosome_length": length,
                        "bin": bin_number,
                        "bin_start_1based_inclusive": bin_start,
                        "bin_end_1based_inclusive": bin_end,
                        "gff_gene_opportunities": denominator,
                        "positive_loss_fragments": count,
                        "positive_loss_fragments_per_gff_gene": count / denominator
                        if denominator
                        else "",
                    }
                )

            if legacy_reproduction:
                for interval_number, interval_start, interval_end in _legacy_nested_bins(
                    length, number_of_bins
                ):
                    denominator = sum(
                        int(gene["start"]) >= interval_start
                        and int(gene["end"]) <= interval_end
                        for gene in gene_list
                    )
                    count = sum(
                        coordinate.start >= interval_start and coordinate.end <= interval_end
                        for coordinate in loss_list
                    )
                    legacy_rows.append(
                        {
                            "analysis_label": analysis_label,
                            "analysis_mode": "manuscript_era_nested_midpoint_reproduction_only",
                            "assembly_unit_id": unit,
                            "biological_species": assembly.biological_species,
                            "haplotype_or_subgenome": assembly.haplotype_or_subgenome,
                            "chromosome": chromosome,
                            "chromosome_length": length,
                            "nested_interval": interval_number,
                            "interval_start_legacy_coordinate": interval_start,
                            "interval_end_legacy_coordinate": interval_end,
                            "gff_gene_opportunities": denominator,
                            "positive_loss_fragments": count,
                            "positive_loss_fragments_per_gff_gene": count / denominator
                            if denominator
                            else "",
                            "intervals_are_mutually_exclusive": "false",
                            "inferential_test_permitted": "false",
                        }
                    )

            chromosome_rows.append(
                {
                    "analysis_label": analysis_label,
                    "assembly_unit_id": unit,
                    "biological_species": assembly.biological_species,
                    "haplotype_or_subgenome": assembly.haplotype_or_subgenome,
                    "assembly_scope": assembly.assembly_scope,
                    "chromosome": chromosome,
                    "chromosome_length": length,
                    "gff_gene_opportunities": len(gene_list),
                    "positive_loss_fragments": len(loss_list),
                    "positive_loss_fragments_per_gff_gene": len(loss_list) / len(gene_list),
                    "mean_midpoint_distance_to_nearest_end_bp": _mean(chromosome_end_distances),
                    "median_midpoint_distance_to_nearest_end_bp": _median(chromosome_end_distances),
                    "mean_normalized_end_distance_0_end_1_center": _mean(
                        chromosome_normalized_end_distances
                    ),
                    "centromere_status": "independently_supplied_interval"
                    if centromere
                    else "not_supplied_for_chromosome",
                    "centromere_start": centromere.start if centromere else "",
                    "centromere_end": centromere.end if centromere else "",
                    "centromere_evidence_source": centromere.evidence_source if centromere else "",
                    "mean_midpoint_distance_to_centromere_bp": _mean(
                        chromosome_centromere_distances
                    ),
                }
            )

        for bin_number in range(1, number_of_bins + 1):
            denominator = end_gene_counts[bin_number]
            count = end_loss_counts[bin_number]
            end_bin_rows.append(
                {
                    "analysis_label": analysis_label,
                    "analysis_mode": "primary_mutually_exclusive_normalized_end_distance",
                    "assembly_unit_id": unit,
                    "biological_species": assembly.biological_species,
                    "haplotype_or_subgenome": assembly.haplotype_or_subgenome,
                    "assembly_scope": assembly.assembly_scope,
                    "end_distance_bin": bin_number,
                    "normalized_end_distance_start_inclusive": format_number(
                        (bin_number - 1) / number_of_bins
                    ),
                    "normalized_end_distance_end_inclusive_only_for_last_bin": format_number(
                        bin_number / number_of_bins
                    ),
                    "gff_gene_opportunities": denominator,
                    "positive_loss_fragments": count,
                    "positive_loss_fragments_per_gff_gene": count / denominator
                    if denominator
                    else "",
                }
            )

        total_genes = len(genes)
        total_losses = len(unit_coordinates)
        unit_rows.append(
            {
                "analysis_label": analysis_label,
                "assembly_unit_id": unit,
                "biological_species": assembly.biological_species,
                "haplotype_or_subgenome": assembly.haplotype_or_subgenome,
                "assembly_scope": assembly.assembly_scope,
                "genome_sequence_count": len(lengths),
                "gff_gene_bearing_sequence_count": len(genes_by_chromosome),
                "gff_gene_opportunities": total_genes,
                "positive_loss_fragments": total_losses,
                "positive_loss_fragments_per_gff_gene": total_losses / total_genes,
                "mean_midpoint_distance_to_nearest_end_bp": _mean(unit_end_distances),
                "median_midpoint_distance_to_nearest_end_bp": _median(unit_end_distances),
                "mean_normalized_end_distance_0_end_1_center": _mean(
                    unit_normalized_end_distances
                ),
                "centromere_intervals_supplied": supplied_centromere_chromosomes,
                "centromere_intervals_missing": len(genes_by_chromosome)
                - supplied_centromere_chromosomes,
            }
        )
        metadata_units.append(
            {
                "assembly_unit_id": unit,
                "biological_species": assembly.biological_species,
                "haplotype_or_subgenome": assembly.haplotype_or_subgenome,
                "assembly_scope": assembly.assembly_scope,
                "genome_filename": assembly.genome.name,
                "gff_filename": assembly.gff.name,
                "genome_declared_sha256": assembly.genome_declared_sha256,
                "gff_declared_sha256": assembly.gff_declared_sha256,
                "genome_sequence_count": len(lengths),
                "gff_gene_bearing_sequence_count": len(genes_by_chromosome),
                "gff_gene_count": len(genes),
                "positive_loss_fragment_count": len(unit_coordinates),
                "centromere_intervals_supplied": supplied_centromere_chromosomes,
            }
        )

    if len(positions_rows) != len(positive_calls):
        raise SchemaError(
            f"internal reconciliation failure: emitted {len(positions_rows)} loss positions for "
            f"{len(positive_calls)} positive calls"
        )
    if sum(int(row["positive_loss_fragments"]) for row in equal_bin_rows) != len(
        positive_calls
    ):
        raise SchemaError("internal reconciliation failure: equal-width bins do not sum to positive calls")
    if sum(int(row["positive_loss_fragments"]) for row in end_bin_rows) != len(
        positive_calls
    ):
        raise SchemaError("internal reconciliation failure: end-distance bins do not sum to positive calls")

    output = Path(output_dir).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise SchemaError(f"output directory is not empty; refusing to mix runs: {output}")
    output.mkdir(parents=True, exist_ok=True)
    positions_path = write_tsv(
        output / "loss_positions.tsv",
        positions_rows,
        [
            "analysis_label",
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "assembly_scope",
            "reference_gene_id",
            "classification",
            "chromosome",
            "feature_start",
            "feature_end",
            "feature_midpoint",
            "chromosome_length",
            "equal_width_bin",
            "midpoint_distance_to_nearest_end_bp",
            "normalized_end_distance_0_end_1_center",
            "end_distance_bin_1_end_to_n_center",
            "nearest_chromosome_end",
            "centromere_start",
            "centromere_end",
            "midpoint_distance_to_centromere_bp",
            "centromere_distance_fraction_of_chromosome",
            "centromere_status",
            "centromere_evidence_source",
        ],
    )
    chromosomes_path = write_tsv(
        output / "chromosome_summary.tsv",
        chromosome_rows,
        [
            "analysis_label",
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "assembly_scope",
            "chromosome",
            "chromosome_length",
            "gff_gene_opportunities",
            "positive_loss_fragments",
            "positive_loss_fragments_per_gff_gene",
            "mean_midpoint_distance_to_nearest_end_bp",
            "median_midpoint_distance_to_nearest_end_bp",
            "mean_normalized_end_distance_0_end_1_center",
            "centromere_status",
            "centromere_start",
            "centromere_end",
            "centromere_evidence_source",
            "mean_midpoint_distance_to_centromere_bp",
        ],
    )
    equal_bins_path = write_tsv(
        output / "equal_width_bins.tsv",
        equal_bin_rows,
        [
            "analysis_label",
            "analysis_mode",
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "chromosome",
            "chromosome_length",
            "bin",
            "bin_start_1based_inclusive",
            "bin_end_1based_inclusive",
            "gff_gene_opportunities",
            "positive_loss_fragments",
            "positive_loss_fragments_per_gff_gene",
        ],
    )
    end_bins_path = write_tsv(
        output / "end_distance_bins.tsv",
        end_bin_rows,
        [
            "analysis_label",
            "analysis_mode",
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "assembly_scope",
            "end_distance_bin",
            "normalized_end_distance_start_inclusive",
            "normalized_end_distance_end_inclusive_only_for_last_bin",
            "gff_gene_opportunities",
            "positive_loss_fragments",
            "positive_loss_fragments_per_gff_gene",
        ],
    )
    units_path = write_tsv(
        output / "assembly_unit_summary.tsv",
        unit_rows,
        [
            "analysis_label",
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "assembly_scope",
            "genome_sequence_count",
            "gff_gene_bearing_sequence_count",
            "gff_gene_opportunities",
            "positive_loss_fragments",
            "positive_loss_fragments_per_gff_gene",
            "mean_midpoint_distance_to_nearest_end_bp",
            "median_midpoint_distance_to_nearest_end_bp",
            "mean_normalized_end_distance_0_end_1_center",
            "centromere_intervals_supplied",
            "centromere_intervals_missing",
        ],
    )
    outputs: dict[str, Path] = {
        "positions": positions_path,
        "chromosomes": chromosomes_path,
        "equal_width_bins": equal_bins_path,
        "end_distance_bins": end_bins_path,
        "assembly_units": units_path,
    }
    if legacy_reproduction:
        legacy_path = write_tsv(
            output / "legacy_nested_midpoint_intervals.tsv",
            legacy_rows,
            [
                "analysis_label",
                "analysis_mode",
                "assembly_unit_id",
                "biological_species",
                "haplotype_or_subgenome",
                "chromosome",
                "chromosome_length",
                "nested_interval",
                "interval_start_legacy_coordinate",
                "interval_end_legacy_coordinate",
                "gff_gene_opportunities",
                "positive_loss_fragments",
                "positive_loss_fragments_per_gff_gene",
                "intervals_are_mutually_exclusive",
                "inferential_test_permitted",
            ],
        )
        outputs["legacy_reproduction"] = legacy_path

    input_records = {
        "positive_calls": {
            "filename": Path(positive_calls_path).name,
            "sha256": _sha256_small_input(positive_calls_path),
        },
        "feature_coordinates": {
            "filename": Path(feature_coordinates_path).name,
            "sha256": _sha256_small_input(feature_coordinates_path),
        },
        "assembly_manifest": {
            "filename": Path(assembly_manifest_path).name,
            "sha256": _sha256_small_input(assembly_manifest_path),
        },
    }
    if centromeres_path is not None:
        input_records["centromeres"] = {
            "filename": Path(centromeres_path).name,
            "sha256": _sha256_small_input(centromeres_path),
        }
    metadata: dict[str, object] = {
        "schema_version": "1.0",
        "analysis_label": analysis_label,
        "primary_analysis_mode": "mutually_exclusive_equal_width",
        "positive_classes": sorted(classes),
        "number_of_bins": number_of_bins,
        "coordinate_convention": "1-based inclusive; target-fragment midpoint used for primary bins and distances",
        "normalized_end_distance_definition": (
            "min(midpoint-1, chromosome_length-midpoint) / ((chromosome_length-1)/2); "
            "0 is a chromosome end and 1 is the center"
        ),
        "denominator_definition": f"all unique GFF {gene_feature!r} features with ID on gene-bearing sequences",
        "centromere_policy": (
            "distance computed only from independently supplied intervals; blank otherwise"
            if centromeres_path is not None
            else "no independently supplied centromere intervals; all centromere distances blank"
        ),
        "require_complete_centromeres": require_complete_centromeres,
        "legacy_reproduction": {
            "enabled": legacy_reproduction,
            "mode": "manuscript_era_nested_midpoint_reproduction_only",
            "overlapping_intervals": legacy_reproduction,
            "inferential_test_permitted": False,
        },
        "reconciliation": {
            "assembly_unit_count": len(assemblies),
            "positive_call_count": len(positive_calls),
            "emitted_position_count": len(positions_rows),
            "equal_width_positive_count": sum(
                int(row["positive_loss_fragments"]) for row in equal_bin_rows
            ),
            "end_distance_positive_count": sum(
                int(row["positive_loss_fragments"]) for row in end_bin_rows
            ),
        },
        "inputs": input_records,
        "assembly_units": metadata_units,
        "outputs": {
            key: {"filename": path.name, "sha256": _sha256_small_input(path)}
            for key, path in outputs.items()
        },
    }
    metadata_path = _atomic_json(output / "run_summary.json", metadata)
    outputs["metadata"] = metadata_path
    return outputs
