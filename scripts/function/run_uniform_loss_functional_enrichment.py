#!/usr/bin/env python3
"""Run GO and KEGG enrichment on the unified species-level loss calls.

The primary foregrounds use the same ``deleted + strict pseudogenized``
positive-complete definition as the uniform shared/non-shared analysis.  A
separate sensitivity contains only species-gene calls for which every required
assembly unit is strict pseudogenized.  Uncertain and partial calls never enter
a foreground.  Tests are one-sided hypergeometric over-representation with
Benjamini-Hochberg correction independently within each foreground/ontology.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping


GO_RE = re.compile(r"^GO:\d{7}$")
KO_RE = re.compile(r"^(?:ko:)?(K\d{5})$")
PATHWAY_RE = re.compile(r"^(?:ko|map)(\d{5})$")
GO_ROOTS = {"GO:0003674", "GO:0005575", "GO:0008150"}
ONTOLOGIES = ("GO", "KEGG_KO", "KEGG_PATHWAY")


class EnrichmentError(ValueError):
    """Raised when a frozen input cannot support the declared analysis."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-protein", required=True, type=Path)
    parser.add_argument("--species-matrix", required=True, type=Path)
    parser.add_argument("--shared-genes", required=True, type=Path)
    parser.add_argument("--species-summary", required=True, type=Path)
    parser.add_argument("--emapper-annotations", required=True, type=Path)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--ko-descriptions", required=True, type=Path)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
    parser.add_argument("--expected-species", type=int, default=13)
    parser.add_argument("--expected-shared-genes", type=int, default=287)
    parser.add_argument("--expected-species-matrix-rows", type=int, default=462111)
    parser.add_argument("--expected-reference-protein-sha256", default="")
    parser.add_argument("--expected-emapper-sha256", default="")
    parser.add_argument("--expected-go-obo-sha256", default="")
    parser.add_argument("--expected-ko-descriptions-sha256", default="")
    parser.add_argument("--minimum-background-count", type=int, default=5)
    parser.add_argument("--minimum-significant-study-count", type=int, default=2)
    parser.add_argument("--fdr-cutoff", type=float, default=0.05)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.startswith(">"):
                continue
            identifier = line[1:].split(None, 1)[0]
            if not identifier or identifier in seen:
                raise EnrichmentError(
                    f"{path.name}:{line_number}: empty or duplicate FASTA identifier"
                )
            seen.add(identifier)
            identifiers.append(identifier)
    if not identifiers:
        raise EnrichmentError(f"{path.name}: no FASTA records")
    return identifiers


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields or len(fields) != len(set(fields)):
            raise EnrichmentError(f"{path.name}: invalid or duplicate TSV header")
        rows = [dict(row) for row in reader]
    if not rows:
        raise EnrichmentError(f"{path.name}: no data rows")
    return rows, fields


def require_columns(path: Path, fields: Iterable[str], required: Iterable[str]) -> None:
    missing = sorted(set(required).difference(fields))
    if missing:
        raise EnrichmentError(f"{path.name}: missing columns: {', '.join(missing)}")


def parse_nonnegative(value: str, *, context: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise EnrichmentError(f"{context}: expected integer, found {value!r}") from error
    if parsed < 0:
        raise EnrichmentError(f"{context}: expected nonnegative integer")
    return parsed


def read_shared(path: Path, expected: int) -> set[str]:
    rows, fields = read_tsv(path)
    require_columns(path, fields, ("reference_gene_id", "shared_positive_complete"))
    shared: set[str] = set()
    for line_number, row in enumerate(rows, 2):
        gene = row["reference_gene_id"].strip()
        if not gene or gene in shared:
            raise EnrichmentError(f"{path.name}:{line_number}: empty or duplicate gene")
        if row["shared_positive_complete"].strip().lower() != "true":
            raise EnrichmentError(f"{path.name}:{line_number}: non-shared row present")
        shared.add(gene)
    if len(shared) != expected:
        raise EnrichmentError(f"{path.name}: observed {len(shared)} shared genes; expected {expected}")
    return shared


def build_foregrounds(
    path: Path,
    shared: set[str],
    *,
    expected_species: int,
    expected_rows: int,
) -> tuple[dict[str, set[str]], dict[str, dict[str, str]], list[str]]:
    rows, fields = read_tsv(path)
    require_columns(
        path,
        fields,
        (
            "reference_gene_id",
            "biological_species",
            "species_gene_status",
            "species_positive_by_rule",
            "assembly_unit_count",
            "pseudogenized_unit_count",
            "deleted_unit_count",
        ),
    )
    if len(rows) != expected_rows:
        raise EnrichmentError(f"{path.name}: observed {len(rows)} rows; expected {expected_rows}")
    species = sorted({row["biological_species"].strip() for row in rows})
    if len(species) != expected_species or any(not value for value in species):
        raise EnrichmentError(
            f"{path.name}: observed {len(species)} biological species; expected {expected_species}"
        )
    combined = {taxon: set() for taxon in species}
    strict_pseudo = {taxon: set() for taxon in species}
    seen: set[tuple[str, str]] = set()
    shared_rows = 0
    for line_number, row in enumerate(rows, 2):
        gene = row["reference_gene_id"].strip()
        taxon = row["biological_species"].strip()
        pair = (gene, taxon)
        if not gene or pair in seen:
            raise EnrichmentError(f"{path.name}:{line_number}: empty or duplicate species-gene row")
        seen.add(pair)
        status = row["species_gene_status"].strip()
        positive = row["species_positive_by_rule"].strip().lower()
        units = parse_nonnegative(row["assembly_unit_count"], context=f"{path.name}:{line_number}")
        pseudogenized = parse_nonnegative(
            row["pseudogenized_unit_count"], context=f"{path.name}:{line_number}"
        )
        deleted = parse_nonnegative(row["deleted_unit_count"], context=f"{path.name}:{line_number}")
        if units < 1 or pseudogenized + deleted > units:
            raise EnrichmentError(f"{path.name}:{line_number}: inconsistent unit counts")
        is_positive = status == "positive_complete"
        if positive not in {"true", "false"} or (positive == "true") != is_positive:
            raise EnrichmentError(f"{path.name}:{line_number}: positive flag/status mismatch")
        if gene in shared:
            if not is_positive:
                raise EnrichmentError(f"{path.name}:{line_number}: shared gene is not positive_complete")
            shared_rows += 1
            continue
        if is_positive:
            if pseudogenized + deleted != units:
                raise EnrichmentError(f"{path.name}:{line_number}: positive_complete does not close")
            combined[taxon].add(gene)
            if pseudogenized == units:
                strict_pseudo[taxon].add(gene)
    if shared_rows != len(shared) * len(species):
        raise EnrichmentError(f"{path.name}: shared species-gene grid is incomplete")

    gene_sets: dict[str, set[str]] = {"shared_combined": set(shared)}
    metadata: dict[str, dict[str, str]] = {
        "shared_combined": {
            "analysis_scope": "shared_positive_complete",
            "biological_species": "all_13_lineages",
            "evidence_mode": "deleted_plus_strict_pseudogenized",
            "background_scope": "all_reference_genes",
        }
    }
    pooled_combined = set().union(*combined.values())
    pooled_pseudo = set().union(*strict_pseudo.values())
    gene_sets["pooled_nonshared_combined"] = pooled_combined
    metadata["pooled_nonshared_combined"] = {
        "analysis_scope": "pooled_nonshared_positive_complete",
        "biological_species": "all_13_lineages",
        "evidence_mode": "deleted_plus_strict_pseudogenized",
        "background_scope": "reference_genes_excluding_shared",
    }
    gene_sets["pooled_nonshared_strict_pseudogenized"] = pooled_pseudo
    metadata["pooled_nonshared_strict_pseudogenized"] = {
        "analysis_scope": "pooled_nonshared_positive_complete_sensitivity",
        "biological_species": "all_13_lineages",
        "evidence_mode": "all_required_units_strict_pseudogenized",
        "background_scope": "reference_genes_excluding_shared",
    }
    for index, taxon in enumerate(species, 1):
        combined_id = f"lineage_{index:02d}_nonshared_combined"
        pseudo_id = f"lineage_{index:02d}_nonshared_strict_pseudogenized"
        gene_sets[combined_id] = combined[taxon]
        gene_sets[pseudo_id] = strict_pseudo[taxon]
        metadata[combined_id] = {
            "analysis_scope": "lineage_nonshared_positive_complete",
            "biological_species": taxon,
            "evidence_mode": "deleted_plus_strict_pseudogenized",
            "background_scope": "reference_genes_excluding_shared",
        }
        metadata[pseudo_id] = {
            "analysis_scope": "lineage_nonshared_positive_complete_sensitivity",
            "biological_species": taxon,
            "evidence_mode": "all_required_units_strict_pseudogenized",
            "background_scope": "reference_genes_excluding_shared",
        }
    return gene_sets, metadata, species


def read_emapper(path: Path) -> tuple[dict[str, dict[str, set[str]]], dict[str, int], list[str]]:
    header: list[str] | None = None
    metadata: list[str] = []
    associations: dict[str, dict[str, set[str]]] = {}
    with path.open(encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line.startswith("##"):
                metadata.append(line[2:].strip())
                continue
            if line.startswith("#"):
                candidate = line[1:].split("\t")
                if {"query", "GOs", "KEGG_ko", "KEGG_Pathway"}.issubset(candidate):
                    header = candidate
                continue
            if not line:
                continue
            if header is None:
                raise EnrichmentError(f"{path.name}:{line_number}: data before annotation header")
            values = line.split("\t")
            if len(values) != len(header):
                raise EnrichmentError(f"{path.name}:{line_number}: column-count mismatch")
            row = dict(zip(header, values))
            gene = row["query"].strip()
            if not gene or gene in associations:
                raise EnrichmentError(f"{path.name}:{line_number}: empty or duplicate query")
            go_terms: set[str] = set()
            for token in row["GOs"].split(",") if row["GOs"] not in {"", "-"} else ():
                token = token.strip()
                if not GO_RE.fullmatch(token):
                    raise EnrichmentError(f"{path.name}:{line_number}: invalid GO term {token!r}")
                if token not in GO_ROOTS:
                    go_terms.add(token)
            ko_terms: set[str] = set()
            for token in row["KEGG_ko"].split(",") if row["KEGG_ko"] not in {"", "-"} else ():
                match = KO_RE.fullmatch(token.strip())
                if match is None:
                    raise EnrichmentError(f"{path.name}:{line_number}: invalid KEGG KO {token!r}")
                ko_terms.add(match.group(1))
            pathways: set[str] = set()
            for token in row["KEGG_Pathway"].split(",") if row["KEGG_Pathway"] not in {"", "-"} else ():
                match = PATHWAY_RE.fullmatch(token.strip())
                if match is None:
                    raise EnrichmentError(f"{path.name}:{line_number}: invalid KEGG pathway {token!r}")
                pathways.add(f"map{match.group(1)}")
            associations[gene] = {
                "GO": go_terms,
                "KEGG_KO": ko_terms,
                "KEGG_PATHWAY": pathways,
            }
    if header is None or not associations:
        raise EnrichmentError(f"{path.name}: no valid eggNOG-mapper annotation rows")
    qc = {"annotation_record_count": len(associations)}
    for ontology in ONTOLOGIES:
        qc[f"genes_with_{ontology.lower()}"] = sum(
            bool(values[ontology]) for values in associations.values()
        )
        qc[f"distinct_{ontology.lower()}_terms"] = len(
            set().union(*(values[ontology] for values in associations.values()))
        )
    return associations, qc, metadata


def read_obo(path: Path) -> tuple[dict[str, tuple[str, str]], str]:
    terms: dict[str, tuple[str, str]] = {}
    current: dict[str, object] | None = None
    data_version = ""
    for line in path.read_text(encoding="utf-8", errors="strict").splitlines():
        if line.startswith("data-version:"):
            data_version = line.split(":", 1)[1].strip()
        if line == "[Term]":
            if current is not None:
                _store_obo_term(current, terms)
            current = {"alt_id": [], "obsolete": False}
            continue
        if line.startswith("["):
            if current is not None:
                _store_obo_term(current, terms)
            current = None
            continue
        if current is None or not line or line.startswith("!"):
            continue
        key, separator, value = line.partition(": ")
        if not separator:
            continue
        if key == "alt_id":
            current["alt_id"].append(value)  # type: ignore[union-attr]
        elif key == "is_obsolete":
            current["obsolete"] = value == "true"
        elif key in {"id", "name", "namespace"}:
            current[key] = value
    if current is not None:
        _store_obo_term(current, terms)
    if not terms or not data_version:
        raise EnrichmentError(f"{path.name}: invalid GO OBO or missing data-version")
    return terms, data_version


def _store_obo_term(current: Mapping[str, object], terms: dict[str, tuple[str, str]]) -> None:
    if current.get("obsolete"):
        return
    identifier = current.get("id")
    name = current.get("name")
    namespace = current.get("namespace")
    if not all(isinstance(value, str) and value for value in (identifier, name, namespace)):
        return
    assert isinstance(identifier, str) and isinstance(name, str) and isinstance(namespace, str)
    record = (name, namespace)
    terms[identifier] = record
    for alternate in current.get("alt_id", []):  # type: ignore[assignment]
        terms[str(alternate)] = record


def read_ko_descriptions(path: Path) -> dict[str, str]:
    descriptions: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, 1):
            columns = line.rstrip("\r\n").split("\t", 1)
            if len(columns) != 2 or not re.fullmatch(r"K\d{5}", columns[0]):
                raise EnrichmentError(f"{path.name}:{line_number}: expected Kxxxxx<TAB>description")
            if columns[0] in descriptions:
                raise EnrichmentError(f"{path.name}:{line_number}: duplicate KO")
            descriptions[columns[0]] = columns[1]
    if not descriptions:
        raise EnrichmentError(f"{path.name}: no KEGG KO descriptions")
    return descriptions


def bh_adjust(values: list[float]) -> list[float]:
    if not values:
        return []
    ordered = sorted(enumerate(values), key=lambda pair: pair[1])
    adjusted = [1.0] * len(values)
    running = 1.0
    for reverse_index in range(len(values) - 1, -1, -1):
        original_index, value = ordered[reverse_index]
        rank = reverse_index + 1
        running = min(running, value * len(values) / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def write_tsv(path: Path, fields: list[str], rows: Iterable[Mapping[str, object]], *, gz: bool = False) -> None:
    opener = gzip.open if gz else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def run() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise EnrichmentError(f"output directory already exists: {args.output_dir}")
    if args.minimum_background_count < 1 or args.minimum_significant_study_count < 1:
        raise EnrichmentError("minimum counts must be positive")
    if not 0 < args.fdr_cutoff <= 1:
        raise EnrichmentError("fdr cutoff must be in (0, 1]")
    inputs = {
        "reference_protein": args.reference_protein,
        "species_matrix": args.species_matrix,
        "shared_genes": args.shared_genes,
        "species_summary": args.species_summary,
        "emapper_annotations": args.emapper_annotations,
        "go_obo": args.go_obo,
        "ko_descriptions": args.ko_descriptions,
    }
    for role, path in inputs.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise EnrichmentError(f"missing or empty {role}: {path}")
    input_hashes = {role: sha256(path) for role, path in inputs.items()}
    expected_hashes = {
        "reference_protein": args.expected_reference_protein_sha256,
        "emapper_annotations": args.expected_emapper_sha256,
        "go_obo": args.expected_go_obo_sha256,
        "ko_descriptions": args.expected_ko_descriptions_sha256,
    }
    for role, expected in expected_hashes.items():
        if expected and input_hashes[role] != expected.lower():
            raise EnrichmentError(
                f"{role} SHA-256 mismatch: expected {expected.lower()}, observed {input_hashes[role]}"
            )

    reference_ids = read_fasta_ids(args.reference_protein)
    if len(reference_ids) != args.expected_reference_genes:
        raise EnrichmentError(
            f"reference protein has {len(reference_ids)} records; expected {args.expected_reference_genes}"
        )
    reference = set(reference_ids)
    shared = read_shared(args.shared_genes, args.expected_shared_genes)
    gene_sets, foreground_metadata, species = build_foregrounds(
        args.species_matrix,
        shared,
        expected_species=args.expected_species,
        expected_rows=args.expected_species_matrix_rows,
    )
    if not set().union(*gene_sets.values()).issubset(reference):
        raise EnrichmentError("foreground contains a gene outside the reference protein set")
    summary_payload = json.loads(args.species_summary.read_text(encoding="utf-8"))
    if summary_payload.get("status") != "PASS" or summary_payload.get("shared_positive_complete_gene_count") != len(shared):
        raise EnrichmentError("species summary is not a matching PASS aggregation")
    associations, annotation_qc, emapper_metadata = read_emapper(args.emapper_annotations)
    if not set(associations).issubset(reference):
        raise EnrichmentError("eggNOG annotation contains IDs outside the reference protein set")
    go_terms, go_version = read_obo(args.go_obo)
    ko_descriptions = read_ko_descriptions(args.ko_descriptions)
    unknown_go = sorted(
        set().union(*(row["GO"] for row in associations.values())).difference(go_terms)
    )
    if unknown_go:
        unresolved_assignments = 0
        for values in associations.values():
            unresolved_assignments += len(values["GO"].intersection(unknown_go))
            values["GO"].difference_update(unknown_go)
        annotation_qc["unresolved_go_terms_excluded"] = len(unknown_go)
        annotation_qc["unresolved_gene_go_assignments_excluded"] = unresolved_assignments

    try:
        import scipy
        from scipy.stats import hypergeom
    except ImportError as error:
        raise EnrichmentError("scipy is required for enrichment") from error

    output_parent = args.output_dir.parent
    output_parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.tmp.", dir=output_parent))
    all_fields = [
        "foreground_id", "analysis_scope", "biological_species", "evidence_mode",
        "ontology", "term_id", "term_name", "go_namespace", "study_count", "study_size",
        "background_count", "background_size", "p_overrepresentation", "p_fdr_bh",
        "fold_enrichment", "significant_fdr", "study_gene_ids",
    ]
    foreground_rows: list[dict[str, object]] = []
    significant_rows: list[dict[str, object]] = []
    all_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    try:
        for foreground_id in sorted(gene_sets):
            meta = foreground_metadata[foreground_id]
            requested = gene_sets[foreground_id]
            for gene in sorted(requested):
                gene_rows.append({"foreground_id": foreground_id, "reference_gene_id": gene})
            for ontology in ONTOLOGIES:
                excluded = shared if meta["background_scope"] == "reference_genes_excluding_shared" else set()
                background = {
                    gene for gene in reference.difference(excluded)
                    if gene in associations and associations[gene][ontology]
                }
                study = requested & background
                term_to_genes: dict[str, set[str]] = defaultdict(set)
                for gene in background:
                    for term in associations[gene][ontology]:
                        term_to_genes[term].add(gene)
                tested_terms = sorted(
                    term for term, genes in term_to_genes.items()
                    if len(genes) >= args.minimum_background_count
                )
                foreground_rows.append(
                    {
                        "foreground_id": foreground_id,
                        **meta,
                        "ontology": ontology,
                        "requested_gene_count": len(requested),
                        "annotated_study_gene_count": len(study),
                        "annotated_background_gene_count": len(background),
                        "tested_term_count": len(tested_terms),
                        "annotation_coverage": len(study) / len(requested) if requested else "",
                    }
                )
                if not study:
                    continue
                raw: list[dict[str, object]] = []
                p_values: list[float] = []
                for term in tested_terms:
                    term_background = term_to_genes[term]
                    hits = study & term_background
                    p_value = float(
                        hypergeom.sf(
                            len(hits) - 1,
                            len(background),
                            len(term_background),
                            len(study),
                        )
                    )
                    fold = (
                        (len(hits) / len(study)) /
                        (len(term_background) / len(background))
                    ) if hits else 0.0
                    if ontology == "GO":
                        term_name, namespace = go_terms[term]
                    elif ontology == "KEGG_KO":
                        term_name, namespace = ko_descriptions.get(term, ""), ""
                    else:
                        term_name, namespace = term, ""
                    raw.append(
                        {
                            "foreground_id": foreground_id,
                            **meta,
                            "ontology": ontology,
                            "term_id": term,
                            "term_name": term_name,
                            "go_namespace": namespace,
                            "study_count": len(hits),
                            "study_size": len(study),
                            "background_count": len(term_background),
                            "background_size": len(background),
                            "p_overrepresentation": p_value,
                            "fold_enrichment": fold,
                            "study_gene_ids": ";".join(sorted(hits)),
                        }
                    )
                    p_values.append(p_value)
                for row, adjusted in zip(raw, bh_adjust(p_values)):
                    row["p_fdr_bh"] = adjusted
                    significant = (
                        adjusted <= args.fdr_cutoff
                        and int(row["study_count"]) >= args.minimum_significant_study_count
                        and float(row["fold_enrichment"]) > 1.0
                    )
                    row["significant_fdr"] = str(significant).lower()
                    all_rows.append(row)
                    if significant:
                        significant_rows.append(row)
        all_rows.sort(key=lambda row: (str(row["foreground_id"]), str(row["ontology"]), float(row["p_fdr_bh"]), str(row["term_id"])))
        significant_rows.sort(key=lambda row: (str(row["foreground_id"]), str(row["ontology"]), float(row["p_fdr_bh"]), str(row["term_id"])))
        foreground_fields = [
            "foreground_id", "analysis_scope", "biological_species", "evidence_mode",
            "background_scope", "ontology", "requested_gene_count", "annotated_study_gene_count",
            "annotated_background_gene_count", "tested_term_count", "annotation_coverage",
        ]
        write_tsv(temporary / "foreground_summary.tsv", foreground_fields, foreground_rows)
        write_tsv(temporary / "foreground_gene_ids.tsv.gz", ["foreground_id", "reference_gene_id"], gene_rows, gz=True)
        write_tsv(temporary / "enrichment_all_terms.tsv.gz", all_fields, all_rows, gz=True)
        write_tsv(temporary / "enrichment_significant.tsv", all_fields, significant_rows)
        output_files = [
            temporary / "foreground_summary.tsv",
            temporary / "foreground_gene_ids.tsv.gz",
            temporary / "enrichment_all_terms.tsv.gz",
            temporary / "enrichment_significant.tsv",
        ]
        significant_counts: dict[str, int] = defaultdict(int)
        for row in significant_rows:
            significant_counts[f"{row['foreground_id']}|{row['ontology']}"] += 1
        summary = {
            "schema_version": "1.0",
            "status": "PASS_UNIFORM_GO_KEGG_ENRICHMENT",
            "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "definitions": {
                "primary_positive": "species positive-complete from deleted + strict pseudogenized unit calls",
                "strict_pseudogenized_sensitivity": "every required assembly unit is strict pseudogenized",
                "uncertain_and_partial": "excluded from every foreground",
                "go_terms": "direct frozen eggNOG-mapper GO field; roots excluded; no additional ancestor propagation",
                "unresolved_go_terms": "annotation terms absent from the checksum-bound GO OBO are excluded and counted",
                "kegg_ko": "offline frozen eggNOG-mapper KEGG_ko field",
                "kegg_pathway": "offline frozen eggNOG-mapper KEGG_Pathway field; ko/map duplicate prefixes collapsed to mapNNNNN",
                "test": "one-sided hypergeometric over-representation; BH within foreground and ontology",
            },
            "expected_counts": {
                "reference_gene_count": len(reference),
                "biological_species_count": len(species),
                "shared_positive_complete_gene_count": len(shared),
                "foreground_count": len(gene_sets),
            },
            "annotation_qc": annotation_qc,
            "go_obo_data_version": go_version,
            "emapper_metadata": [
                line for line in emapper_metadata
                if line and "/" not in line and "--" not in line
            ],
            "parameters": {
                "minimum_background_count": args.minimum_background_count,
                "minimum_significant_study_count": args.minimum_significant_study_count,
                "fdr_cutoff": args.fdr_cutoff,
            },
            "significant_term_counts": dict(sorted(significant_counts.items())),
            "inputs": [
                {"role": role, "basename": path.name, "sha256": input_hashes[role]}
                for role, path in inputs.items()
            ],
            "outputs": [
                {"basename": path.name, "sha256": sha256(path)} for path in output_files
            ],
            "software": {"python": os.sys.version.split()[0], "scipy": scipy.__version__},
            "checks": {
                "reference_gene_universe_exact": True,
                "species_grid_exact": True,
                "shared_set_exact": True,
                "foregrounds_exclude_uncertain_and_partial": True,
                "strict_pseudogenized_sensitivity_separate": True,
                "annotation_ids_subset_reference": True,
                "offline_frozen_annotations_only": True,
            },
        }
        (temporary / "functional_enrichment_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        temporary.replace(args.output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def main() -> int:
    try:
        run()
    except (EnrichmentError, OSError, UnicodeError, csv.Error, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
