#!/usr/bin/env python3
"""Test GO/KEGG enrichment among 23 terminal loss-event foregrounds.

For each assembly-unit terminal, the response is whether a reference gene was
assigned to that terminal by the frozen maximal-positive-clade event model.
Genes already lost on an ancestor of the terminal are removed from that
terminal's risk set.  The nuisance model adjusts for log2 four-tissue mean TPM
and log2 reference CD-HIT family size, including quadratic terms.  A one-sided
efficient score test is then applied to every eligible GO/KEGG term.  Terms
passing BH q <= 0.05 are refit with the complete logistic model to report an
adjusted odds ratio and Wald confidence interval.
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

import numpy as np
import scipy
from scipy.special import expit
from scipy.stats import hypergeom, norm

import run_uniform_loss_functional_enrichment as base


CLUSTER_MEMBER = re.compile(r">(.+?)\.\.\.")
ONTOLOGIES = ("GO", "KEGG_KO", "KEGG_PATHWAY")
GO_ROOTS = {"GO:0003674", "GO:0005575", "GO:0008150"}


class AdjustedEnrichmentError(ValueError):
    """Raised when exact inputs cannot support the adjusted analysis."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", required=True, type=Path)
    parser.add_argument("--nodes", required=True, type=Path)
    parser.add_argument("--resolved-background", required=True, type=Path)
    parser.add_argument("--unit-metadata", required=True, type=Path)
    parser.add_argument("--four-tissue-tpm", required=True, type=Path)
    parser.add_argument("--cdhit-clusters", required=True, type=Path)
    parser.add_argument("--emapper-annotations", required=True, type=Path)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--ko-descriptions", required=True, type=Path)
    parser.add_argument("--pathway-names", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--analysis-level",
        choices=("assembly_unit", "biological_species"),
        default="assembly_unit",
    )
    parser.add_argument("--terminal-node-type", default="unit_terminal")
    parser.add_argument("--expected-background-genes", type=int, default=33998)
    parser.add_argument("--expected-terminals", type=int, default=23)
    parser.add_argument("--expected-terminal-event-memberships", type=int, default=41302)
    parser.add_argument("--minimum-background-count", type=int, default=5)
    parser.add_argument("--minimum-terminal-loss-count", type=int, default=2)
    parser.add_argument("--fdr-cutoff", type=float, default=0.05)
    parser.add_argument("--expected-events-sha256", default="")
    parser.add_argument("--expected-tpm-sha256", default="")
    parser.add_argument("--expected-clusters-sha256", default="")
    parser.add_argument("--expected-emapper-sha256", default="")
    parser.add_argument("--expected-go-obo-sha256", default="")
    parser.add_argument("--expected-ko-descriptions-sha256", default="")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = [dict(row) for row in reader]
    if not fields or len(fields) != len(set(fields)) or not rows:
        raise AdjustedEnrichmentError(f"{path.name}: invalid or empty TSV")
    return rows, fields


def require_columns(
    path: Path, fields: Iterable[str], required: Iterable[str]
) -> None:
    missing = sorted(set(required).difference(fields))
    if missing:
        raise AdjustedEnrichmentError(
            f"{path.name}: missing columns: {', '.join(missing)}"
        )


def write_tsv(
    path: Path,
    fields: list[str],
    rows: Iterable[Mapping[str, object]],
    *,
    gz: bool = False,
) -> None:
    opener = gzip.open if gz else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_background(path: Path, expected: int) -> set[str]:
    genes = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(genes) != expected or len(set(genes)) != expected:
        raise AdjustedEnrichmentError(
            f"{path.name}: expected {expected} unique genes, found {len(set(genes))}"
        )
    return set(genes)


def read_metadata(
    path: Path, expected_units: int
) -> dict[str, dict[str, str]]:
    rows, fields = read_tsv(path)
    require_columns(
        path,
        fields,
        (
            "assembly_unit_id",
            "biological_species",
            "haplotype_or_subgenome",
            "include",
        ),
    )
    selected = [row for row in rows if row["include"].strip().lower() == "true"]
    result = {row["assembly_unit_id"]: row for row in selected}
    if len(result) != expected_units or len(selected) != expected_units:
        raise AdjustedEnrichmentError(
            f"{path.name}: expected {expected_units} unique included units"
        )
    return result


def read_tree_inputs(
    events_path: Path,
    nodes_path: Path,
    background: set[str],
    metadata: Mapping[str, Mapping[str, str]],
    expected_terminals: int,
    expected_memberships: int,
    terminal_node_type: str,
) -> tuple[
    list[dict[str, str]],
    dict[str, set[str]],
    dict[str, set[str]],
    dict[str, str],
]:
    node_rows, node_fields = read_tsv(nodes_path)
    require_columns(
        nodes_path,
        node_fields,
        ("node_id", "parent_node_id", "node_type", "descendant_units"),
    )
    if len({row["node_id"] for row in node_rows}) != len(node_rows):
        raise AdjustedEnrichmentError(f"{nodes_path.name}: duplicate node ID")
    node_by_id = {row["node_id"]: row for row in node_rows}
    parent = {row["node_id"]: row["parent_node_id"] for row in node_rows}
    terminals = [
        row for row in node_rows if row["node_type"] == terminal_node_type
    ]
    if len(terminals) != expected_terminals:
        raise AdjustedEnrichmentError(
            f"{nodes_path.name}: expected {expected_terminals} terminal nodes"
        )
    terminal_units: dict[str, str] = {}
    for row in terminals:
        units = [token for token in row["descendant_units"].split(";") if token]
        if len(units) != 1:
            raise AdjustedEnrichmentError(
                f"{row['node_id']}: terminal must contain exactly one assembly unit"
            )
        terminal_units[row["node_id"]] = units[0]
    if set(terminal_units.values()) != set(metadata):
        raise AdjustedEnrichmentError(
            "terminal-node unit set differs from included unit metadata"
        )

    event_rows, event_fields = read_tsv(events_path)
    require_columns(
        events_path,
        event_fields,
        ("reference_gene_id", "node_id", "node_type", "descendant_units"),
    )
    events_by_node: dict[str, set[str]] = defaultdict(set)
    seen: set[tuple[str, str]] = set()
    for line_number, row in enumerate(event_rows, 2):
        gene = row["reference_gene_id"]
        node_id = row["node_id"]
        key = (gene, node_id)
        if gene not in background:
            raise AdjustedEnrichmentError(
                f"{events_path.name}:{line_number}: event gene outside resolved background"
            )
        if node_id not in node_by_id or key in seen:
            raise AdjustedEnrichmentError(
                f"{events_path.name}:{line_number}: invalid or duplicate event"
            )
        seen.add(key)
        events_by_node[node_id].add(gene)
    terminal_memberships = sum(
        len(events_by_node[row["node_id"]]) for row in terminals
    )
    if terminal_memberships != expected_memberships:
        raise AdjustedEnrichmentError(
            f"terminal event memberships changed: {terminal_memberships} != "
            f"{expected_memberships}"
        )

    risk_by_terminal: dict[str, set[str]] = {}
    for row in terminals:
        node_id = row["node_id"]
        ancestor_losses: set[str] = set()
        ancestor = parent[node_id]
        visited: set[str] = set()
        while ancestor:
            if ancestor in visited or ancestor not in node_by_id:
                raise AdjustedEnrichmentError(f"{node_id}: invalid ancestor chain")
            visited.add(ancestor)
            ancestor_losses.update(events_by_node.get(ancestor, set()))
            ancestor = parent[ancestor]
        risk = background - ancestor_losses
        foreground = events_by_node.get(node_id, set())
        if not foreground or not foreground.issubset(risk):
            raise AdjustedEnrichmentError(
                f"{node_id}: empty foreground or foreground outside terminal risk set"
            )
        risk_by_terminal[node_id] = risk
    return terminals, events_by_node, risk_by_terminal, terminal_units


def read_tpm(path: Path) -> dict[str, float]:
    rows, fields = read_tsv(path)
    require_columns(path, fields, ("reference_gene_id", "four_tissue_mean_tpm"))
    result: dict[str, float] = {}
    for line_number, row in enumerate(rows, 2):
        gene = row["reference_gene_id"]
        try:
            value = float(row["four_tissue_mean_tpm"])
        except ValueError as error:
            raise AdjustedEnrichmentError(
                f"{path.name}:{line_number}: invalid mean TPM"
            ) from error
        if not gene or gene in result or not math.isfinite(value) or value < 0:
            raise AdjustedEnrichmentError(
                f"{path.name}:{line_number}: invalid or duplicate TPM row"
            )
        result[gene] = value
    return result


def read_family_sizes(path: Path) -> dict[str, int]:
    clusters: list[list[str]] = []
    current: list[str] | None = None
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">Cluster"):
                current = []
                clusters.append(current)
                continue
            match = CLUSTER_MEMBER.search(line)
            if current is None or match is None:
                raise AdjustedEnrichmentError(
                    f"{path.name}:{line_number}: malformed CD-HIT record"
                )
            current.append(match.group(1))
    result: dict[str, int] = {}
    for cluster in clusters:
        if not cluster:
            raise AdjustedEnrichmentError(f"{path.name}: empty CD-HIT cluster")
        size = len(cluster)
        for gene in cluster:
            if gene in result:
                raise AdjustedEnrichmentError(
                    f"{path.name}: duplicate CD-HIT member {gene}"
                )
            result[gene] = size
    if not result:
        raise AdjustedEnrichmentError(f"{path.name}: no CD-HIT members")
    return result


def read_pathway_names(path: Path) -> dict[str, str]:
    rows, fields = read_tsv(path)
    require_columns(path, fields, ("term_id", "term_name"))
    result: dict[str, str] = {}
    for row in rows:
        term = row["term_id"]
        name = row["term_name"]
        if not re.fullmatch(r"map\d{5}", term) or not name or term in result:
            raise AdjustedEnrichmentError(f"{path.name}: invalid pathway-name row")
        result[term] = name
    return result


def build_term_maps(
    associations: dict[str, dict[str, set[str]]],
    go_terms: Mapping[str, tuple[str, str]],
    ko_descriptions: Mapping[str, str],
    pathway_names: Mapping[str, str],
) -> tuple[
    dict[str, dict[str, set[str]]],
    dict[tuple[str, str], tuple[str, str]],
    dict[str, set[str]],
    int,
]:
    unknown_go = set().union(
        *(values["GO"] for values in associations.values())
    ).difference(go_terms)
    unresolved_assignments = 0
    if unknown_go:
        for values in associations.values():
            unresolved_assignments += len(values["GO"].intersection(unknown_go))
            values["GO"].difference_update(unknown_go)
    term_genes: dict[str, dict[str, set[str]]] = {
        ontology: defaultdict(set) for ontology in ONTOLOGIES
    }
    annotated_genes: dict[str, set[str]] = {
        ontology: set() for ontology in ONTOLOGIES
    }
    term_meta: dict[tuple[str, str], tuple[str, str]] = {}
    for gene, values in associations.items():
        for ontology in ONTOLOGIES:
            terms = values[ontology]
            if terms:
                annotated_genes[ontology].add(gene)
            for term in terms:
                term_genes[ontology][term].add(gene)
                if ontology == "GO":
                    term_meta[(ontology, term)] = go_terms[term]
                elif ontology == "KEGG_KO":
                    term_meta[(ontology, term)] = (
                        ko_descriptions.get(term, term),
                        "",
                    )
                else:
                    term_meta[(ontology, term)] = (
                        pathway_names.get(term, term),
                        "",
                    )
    return term_genes, term_meta, annotated_genes, unresolved_assignments


def standardize(values: np.ndarray) -> np.ndarray:
    mean = float(values.mean())
    scale = float(values.std(ddof=0))
    if not math.isfinite(scale) or scale <= 0:
        raise AdjustedEnrichmentError("covariate has zero or invalid variance")
    return (values - mean) / scale


def design_matrix(
    genes: list[str],
    tpm: Mapping[str, float],
    family_sizes: Mapping[str, int],
) -> np.ndarray:
    expression = np.asarray(
        [math.log2(tpm[gene] + 0.1) for gene in genes], dtype=float
    )
    family = np.asarray(
        [math.log2(family_sizes[gene]) for gene in genes], dtype=float
    )
    expression_z = standardize(expression)
    family_z = standardize(family)
    expression_sq = standardize(expression_z**2)
    family_sq = standardize(family_z**2)
    return np.column_stack(
        [
            np.ones(len(genes), dtype=float),
            expression_z,
            expression_sq,
            family_z,
            family_sq,
        ]
    )


def log_likelihood(design: np.ndarray, outcome: np.ndarray, beta: np.ndarray) -> float:
    eta = design @ beta
    return float(np.sum(outcome * eta - np.logaddexp(0.0, eta)))


def fit_logistic(
    design: np.ndarray,
    outcome: np.ndarray,
    *,
    start: np.ndarray | None = None,
    max_iterations: int = 100,
    tolerance: float = 1e-9,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, bool, int]:
    if len(np.unique(outcome)) != 2:
        raise AdjustedEnrichmentError("logistic outcome lacks both classes")
    beta = (
        np.zeros(design.shape[1], dtype=float)
        if start is None
        else np.asarray(start, dtype=float).copy()
    )
    converged = False
    old_ll = log_likelihood(design, outcome, beta)
    for iteration in range(1, max_iterations + 1):
        probability = np.clip(expit(design @ beta), 1e-10, 1 - 1e-10)
        weight = np.clip(probability * (1 - probability), 1e-10, None)
        information = design.T @ (weight[:, None] * design)
        score = design.T @ (outcome - probability)
        try:
            step = np.linalg.solve(information, score)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(information, rcond=1e-10) @ score
        factor = 1.0
        accepted = False
        while factor >= 1 / 1024:
            candidate = beta + factor * step
            new_ll = log_likelihood(design, outcome, candidate)
            if new_ll >= old_ll - 1e-10:
                beta = candidate
                old_ll = new_ll
                accepted = True
                break
            factor /= 2
        if not accepted:
            break
        if float(np.max(np.abs(factor * step))) < tolerance:
            converged = True
            break
    probability = np.clip(expit(design @ beta), 1e-10, 1 - 1e-10)
    weight = np.clip(probability * (1 - probability), 1e-10, None)
    information = design.T @ (weight[:, None] * design)
    covariance = np.linalg.pinv(information, rcond=1e-10)
    return beta, probability, covariance, converged, iteration


def score_term(
    indices: np.ndarray,
    residual: np.ndarray,
    weight: np.ndarray,
    design: np.ndarray,
    nuisance_inverse: np.ndarray,
) -> tuple[float, float, float, float]:
    score = float(residual[indices].sum())
    cross = design[indices].T @ weight[indices]
    information = float(weight[indices].sum() - cross @ nuisance_inverse @ cross)
    if not math.isfinite(information) or information <= 1e-12:
        return score, information, math.nan, math.nan
    z_value = score / math.sqrt(information)
    p_value = float(norm.sf(z_value))
    beta_one_step = score / information
    return score, information, z_value, p_value if math.isfinite(p_value) else math.nan


def safe_exp(value: float) -> float:
    return float(math.exp(max(-30.0, min(30.0, value))))


def format_float(value: float | int | str) -> object:
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.12g}"
    return value


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise AdjustedEnrichmentError(
            f"output directory already exists: {args.output_dir}"
        )
    if not 0 < args.fdr_cutoff <= 1:
        raise AdjustedEnrichmentError("FDR cutoff must be in (0, 1]")
    inputs = {
        "events": args.events,
        "nodes": args.nodes,
        "resolved_background": args.resolved_background,
        "unit_metadata": args.unit_metadata,
        "four_tissue_tpm": args.four_tissue_tpm,
        "cdhit_clusters": args.cdhit_clusters,
        "emapper_annotations": args.emapper_annotations,
        "go_obo": args.go_obo,
        "ko_descriptions": args.ko_descriptions,
        "pathway_names": args.pathway_names,
    }
    for role, path in inputs.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise AdjustedEnrichmentError(f"missing or empty {role}: {path}")
    hashes = {role: sha256(path) for role, path in inputs.items()}
    expected_hashes = {
        "events": args.expected_events_sha256,
        "four_tissue_tpm": args.expected_tpm_sha256,
        "cdhit_clusters": args.expected_clusters_sha256,
        "emapper_annotations": args.expected_emapper_sha256,
        "go_obo": args.expected_go_obo_sha256,
        "ko_descriptions": args.expected_ko_descriptions_sha256,
    }
    for role, expected in expected_hashes.items():
        if expected and hashes[role] != expected.lower():
            raise AdjustedEnrichmentError(
                f"{role} SHA-256 mismatch: {hashes[role]} != {expected.lower()}"
            )

    background = read_background(
        args.resolved_background, args.expected_background_genes
    )
    metadata = read_metadata(args.unit_metadata, args.expected_terminals)
    terminals, events_by_node, risk_by_terminal, terminal_units = read_tree_inputs(
        args.events,
        args.nodes,
        background,
        metadata,
        args.expected_terminals,
        args.expected_terminal_event_memberships,
        args.terminal_node_type,
    )
    tpm = read_tpm(args.four_tissue_tpm)
    family_sizes = read_family_sizes(args.cdhit_clusters)
    complete_covariates = background.intersection(tpm).intersection(family_sizes)
    missing_tpm = sorted(background.difference(tpm))
    missing_family = sorted(background.difference(family_sizes))
    if len(complete_covariates) < len(background) - 100:
        raise AdjustedEnrichmentError("unexpectedly many genes lack model covariates")

    associations, annotation_qc, emapper_metadata = base.read_emapper(
        args.emapper_annotations
    )
    go_terms, go_version = base.read_obo(args.go_obo)
    ko_descriptions = base.read_ko_descriptions(args.ko_descriptions)
    pathway_names = read_pathway_names(args.pathway_names)
    term_genes, term_meta, annotated_genes, unresolved_go_assignments = (
        build_term_maps(
            associations,
            go_terms,
            ko_descriptions,
            pathway_names,
        )
    )

    all_rows: list[dict[str, object]] = []
    model_rows: list[dict[str, object]] = []
    significant_gene_rows: list[dict[str, object]] = []
    terminal_order = {
        row["node_id"]: int(row["minimum_leaf_plot_order"]) for row in terminals
    }
    terminals.sort(key=lambda row: terminal_order[row["node_id"]])

    for terminal in terminals:
        node_id = terminal["node_id"]
        unit = terminal_units[node_id]
        unit_meta = metadata[unit]
        foreground = events_by_node[node_id]
        risk = risk_by_terminal[node_id]
        covariate_risk = risk.intersection(complete_covariates)
        for ontology in ONTOLOGIES:
            eligible_set = covariate_risk.intersection(annotated_genes[ontology])
            genes = sorted(eligible_set)
            gene_index = {gene: index for index, gene in enumerate(genes)}
            outcome = np.asarray(
                [1.0 if gene in foreground else 0.0 for gene in genes], dtype=float
            )
            if int(outcome.sum()) < args.minimum_terminal_loss_count:
                raise AdjustedEnrichmentError(
                    f"{node_id}/{ontology}: too few annotated terminal losses"
                )
            design = design_matrix(genes, tpm, family_sizes)
            beta, probability, covariance, converged, iterations = fit_logistic(
                design, outcome
            )
            weight = np.clip(probability * (1 - probability), 1e-10, None)
            nuisance_information = design.T @ (weight[:, None] * design)
            nuisance_inverse = np.linalg.pinv(nuisance_information, rcond=1e-10)
            residual = outcome - probability
            model_rows.append(
                {
                    "node_id": node_id,
                    "assembly_unit_id": unit,
                    "biological_species": unit_meta["biological_species"],
                    "haplotype_or_subgenome": unit_meta[
                        "haplotype_or_subgenome"
                    ],
                    "ontology": ontology,
                    "resolved_global_background": len(background),
                    "ancestor_loss_genes_excluded": len(background - risk),
                    "terminal_risk_genes": len(risk),
                    "covariate_complete_risk_genes": len(covariate_risk),
                    "annotated_model_background": len(genes),
                    "terminal_event_genes": len(foreground),
                    "annotated_terminal_event_genes": int(outcome.sum()),
                    "event_rate": float(outcome.mean()),
                    "null_model_converged": str(converged).lower(),
                    "null_model_iterations": iterations,
                    "null_log_likelihood": log_likelihood(design, outcome, beta),
                    "intercept": beta[0],
                    "expression_linear_beta": beta[1],
                    "expression_quadratic_beta": beta[2],
                    "family_size_linear_beta": beta[3],
                    "family_size_quadratic_beta": beta[4],
                }
            )

            rows: list[dict[str, object]] = []
            for term_id, members in term_genes[ontology].items():
                indices = np.fromiter(
                    (
                        gene_index[gene]
                        for gene in members
                        if gene in gene_index
                    ),
                    dtype=int,
                )
                background_count = len(indices)
                if background_count < args.minimum_background_count:
                    continue
                study_count = int(outcome[indices].sum())
                if study_count < args.minimum_terminal_loss_count:
                    continue
                score, information, z_value, p_score = score_term(
                    indices,
                    residual,
                    weight,
                    design,
                    nuisance_inverse,
                )
                if not math.isfinite(p_score):
                    continue
                foreground_size = int(outcome.sum())
                background_size = len(genes)
                fold = (study_count / foreground_size) / (
                    background_count / background_size
                )
                p_unadjusted = float(
                    hypergeom.sf(
                        study_count - 1,
                        background_size,
                        background_count,
                        foreground_size,
                    )
                )
                beta_one_step = score / information
                se_one_step = 1 / math.sqrt(information)
                term_name, namespace = term_meta[(ontology, term_id)]
                rows.append(
                    {
                        "node_id": node_id,
                        "assembly_unit_id": unit,
                        "biological_species": unit_meta["biological_species"],
                        "haplotype_or_subgenome": unit_meta[
                            "haplotype_or_subgenome"
                        ],
                        "ontology": ontology,
                        "go_namespace": namespace,
                        "term_id": term_id,
                        "term_name": term_name,
                        "terminal_loss_count": study_count,
                        "terminal_loss_annotated_total": foreground_size,
                        "background_count": background_count,
                        "background_size": background_size,
                        "fold_enrichment_unadjusted": fold,
                        "p_hypergeometric": p_unadjusted,
                        "score": score,
                        "efficient_information": information,
                        "score_z": z_value,
                        "p_score_one_sided": p_score,
                        "score_beta_one_step": beta_one_step,
                        "score_odds_ratio_one_step": safe_exp(beta_one_step),
                        "score_ci95_low_one_step": safe_exp(
                            beta_one_step - 1.959963984540054 * se_one_step
                        ),
                        "score_ci95_high_one_step": safe_exp(
                            beta_one_step + 1.959963984540054 * se_one_step
                        ),
                        "_indices": indices,
                        "_design": design,
                        "_outcome": outcome,
                        "_null_beta": beta,
                        "_genes": genes,
                    }
                )
            adjusted = base.bh_adjust(
                [float(row["p_score_one_sided"]) for row in rows]
            )
            unadjusted = base.bh_adjust(
                [float(row["p_hypergeometric"]) for row in rows]
            )
            for row, q_score, q_unadjusted in zip(rows, adjusted, unadjusted):
                row["q_score_bh"] = q_score
                row["q_hypergeometric_bh"] = q_unadjusted
                significant = (
                    q_score <= args.fdr_cutoff
                    and float(row["score_z"]) > 0
                    and float(row["score_odds_ratio_one_step"]) > 1
                )
                row["significant_adjusted"] = str(significant).lower()
                row["full_fit_converged"] = ""
                row["adjusted_odds_ratio"] = ""
                row["adjusted_ci95_low"] = ""
                row["adjusted_ci95_high"] = ""
                row["full_fit_beta"] = ""
                row["full_fit_se"] = ""
                row["full_fit_wald_p_one_sided"] = ""
                if significant:
                    indices = row["_indices"]
                    design = row["_design"]
                    outcome = row["_outcome"]
                    indicator = np.zeros(len(outcome), dtype=float)
                    indicator[indices] = 1.0
                    full_design = np.column_stack([design, indicator])
                    start = np.append(row["_null_beta"], 0.0)
                    (
                        full_beta,
                        _,
                        full_covariance,
                        full_converged,
                        _,
                    ) = fit_logistic(full_design, outcome, start=start)
                    term_beta = float(full_beta[-1])
                    term_se = math.sqrt(max(0.0, float(full_covariance[-1, -1])))
                    wald_z = (
                        term_beta / term_se if term_se > 0 else math.copysign(math.inf, term_beta)
                    )
                    row["full_fit_converged"] = str(full_converged).lower()
                    row["full_fit_beta"] = term_beta
                    row["full_fit_se"] = term_se
                    row["full_fit_wald_p_one_sided"] = float(norm.sf(wald_z))
                    row["adjusted_odds_ratio"] = safe_exp(term_beta)
                    row["adjusted_ci95_low"] = safe_exp(
                        term_beta - 1.959963984540054 * term_se
                    )
                    row["adjusted_ci95_high"] = safe_exp(
                        term_beta + 1.959963984540054 * term_se
                    )
                    for index in indices:
                        if outcome[index] == 1:
                            significant_gene_rows.append(
                                {
                                    "node_id": node_id,
                                    "assembly_unit_id": unit,
                                    "ontology": ontology,
                                    "term_id": row["term_id"],
                                    "reference_gene_id": genes[index],
                                }
                            )
                for key in ("_indices", "_design", "_outcome", "_null_beta", "_genes"):
                    row.pop(key, None)
                all_rows.append(row)

    all_rows.sort(
        key=lambda row: (
            terminal_order[str(row["node_id"])],
            str(row["ontology"]),
            float(row["q_score_bh"]),
            str(row["term_id"]),
        )
    )
    significant_rows = [
        row for row in all_rows if row["significant_adjusted"] == "true"
    ]
    model_rows.sort(
        key=lambda row: (
            terminal_order[str(row["node_id"])],
            ONTOLOGIES.index(str(row["ontology"])),
        )
    )

    output_fields = [
        "node_id",
        "assembly_unit_id",
        "biological_species",
        "haplotype_or_subgenome",
        "ontology",
        "go_namespace",
        "term_id",
        "term_name",
        "terminal_loss_count",
        "terminal_loss_annotated_total",
        "background_count",
        "background_size",
        "fold_enrichment_unadjusted",
        "p_hypergeometric",
        "q_hypergeometric_bh",
        "score",
        "efficient_information",
        "score_z",
        "p_score_one_sided",
        "q_score_bh",
        "score_beta_one_step",
        "score_odds_ratio_one_step",
        "score_ci95_low_one_step",
        "score_ci95_high_one_step",
        "significant_adjusted",
        "full_fit_converged",
        "full_fit_beta",
        "full_fit_se",
        "full_fit_wald_p_one_sided",
        "adjusted_odds_ratio",
        "adjusted_ci95_low",
        "adjusted_ci95_high",
    ]
    model_fields = list(model_rows[0])
    significant_gene_fields = [
        "node_id",
        "assembly_unit_id",
        "ontology",
        "term_id",
        "reference_gene_id",
    ]

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.tmp.", dir=args.output_dir.parent
        )
    )
    try:
        write_tsv(
            temporary / "terminal_model_summary.tsv",
            model_fields,
            (
                {key: format_float(value) for key, value in row.items()}
                for row in model_rows
            ),
        )
        write_tsv(
            temporary / "terminal_adjusted_enrichment_all_terms.tsv.gz",
            output_fields,
            (
                {key: format_float(value) for key, value in row.items()}
                for row in all_rows
            ),
            gz=True,
        )
        write_tsv(
            temporary / "terminal_adjusted_enrichment_significant.tsv",
            output_fields,
            (
                {key: format_float(value) for key, value in row.items()}
                for row in significant_rows
            ),
        )
        write_tsv(
            temporary / "terminal_adjusted_significant_gene_memberships.tsv.gz",
            significant_gene_fields,
            significant_gene_rows,
            gz=True,
        )
        missing_rows = [
            {"reference_gene_id": gene, "missing_covariate": "four_tissue_mean_tpm"}
            for gene in missing_tpm
        ] + [
            {"reference_gene_id": gene, "missing_covariate": "reference_family_size"}
            for gene in missing_family
        ]
        write_tsv(
            temporary / "missing_covariates.tsv",
            ["reference_gene_id", "missing_covariate"],
            missing_rows,
        )
        outputs = [
            temporary / "terminal_model_summary.tsv",
            temporary / "terminal_adjusted_enrichment_all_terms.tsv.gz",
            temporary / "terminal_adjusted_enrichment_significant.tsv",
            temporary / "terminal_adjusted_significant_gene_memberships.tsv.gz",
            temporary / "missing_covariates.tsv",
        ]
        significant_by_ontology = {
            ontology: sum(
                str(row["ontology"]) == ontology for row in significant_rows
            )
            for ontology in ONTOLOGIES
        }
        manifest = {
            "schema_version": "1.0",
            "status": "PASS_TERMINAL_COVARIATE_ADJUSTED_GO_KEGG",
            "analysis_level": args.analysis_level,
            "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "definitions": {
                "positive_response": (
                    (
                        "article-method complete species loss (all constituent "
                        "assembly units decayed/deleted) assigned to one biological-"
                        "species terminal by the maximal-positive-clade model"
                    )
                    if args.analysis_level == "biological_species"
                    else (
                        "article-method decayed + deleted event assigned to one "
                        "assembly-unit terminal by the frozen maximal-positive-clade "
                        "scaffold"
                    )
                ),
                "terminal_risk_set": (
                    "genes resolved in all 23 assembly units, excluding genes assigned "
                    "to an ancestor loss event of the focal terminal"
                ),
                "expression_covariate": (
                    "standardized log2(four-tissue arithmetic mean TPM + 0.1), "
                    "linear and quadratic terms"
                ),
                "family_size_covariate": (
                    "standardized log2(reference CD-HIT 90% family size), "
                    "linear and quadratic terms"
                ),
                "function_terms": (
                    "direct frozen eggNOG-mapper GO, KEGG KO and KEGG pathway "
                    "assignments; GO roots and unresolved OBO IDs excluded"
                ),
                "test": (
                    "one-sided efficient logistic score test adjusted for expression "
                    "and reference family size; BH within terminal and ontology"
                ),
                "effect_estimate": (
                    "complete logistic refit for score-BH-significant terms; "
                    "score one-step estimate retained for all tested terms"
                ),
            },
            "counts": {
                "terminals": len(terminals),
                "resolved_background_genes": len(background),
                "complete_covariate_genes": len(complete_covariates),
                "terminal_event_memberships": sum(
                    len(events_by_node[row["node_id"]]) for row in terminals
                ),
                "tested_term_rows": len(all_rows),
                "significant_adjusted_rows": len(significant_rows),
                "significant_gene_memberships": len(significant_gene_rows),
                "missing_tpm_genes": len(missing_tpm),
                "missing_family_size_genes": len(missing_family),
            },
            "significant_by_ontology": significant_by_ontology,
            "parameters": {
                "minimum_background_count": args.minimum_background_count,
                "minimum_terminal_loss_count": args.minimum_terminal_loss_count,
                "fdr_cutoff": args.fdr_cutoff,
            },
            "annotation_qc": {
                **annotation_qc,
                "unresolved_go_assignments_excluded": unresolved_go_assignments,
                "go_obo_data_version": go_version,
                "emapper_metadata": [
                    line
                    for line in emapper_metadata
                    if line and "/" not in line and "--" not in line
                ],
            },
            "inputs": [
                {"role": role, "basename": path.name, "sha256": hashes[role]}
                for role, path in inputs.items()
            ],
            "outputs": [
                {
                    "basename": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in outputs
            ],
            "software": {
                "python": os.sys.version.split()[0],
                "numpy": np.__version__,
                "scipy": scipy.__version__,
            },
        }
        (temporary / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(args.output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def main() -> int:
    try:
        run(parse_args())
    except (
        AdjustedEnrichmentError,
        base.EnrichmentError,
        OSError,
        UnicodeError,
        csv.Error,
        json.JSONDecodeError,
        np.linalg.LinAlgError,
    ) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
