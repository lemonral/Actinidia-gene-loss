#!/usr/bin/env python3
"""Run GO/KEGG enrichment for topology-aware article-method loss sets.

Foreground membership and branch placement are prepared by
``prepare_manuscript_method_downstream.py``.  This step performs no new gene
loss classification and no sequence search.  It uses the frozen reference
eggNOG annotations and tests GO, KEGG KO, and KEGG pathway over-representation
with one-sided hypergeometric tests and BH correction within each
foreground/ontology.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import shutil
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping


HERE = Path(__file__).resolve().parent
BASE_SCRIPT = HERE / "run_uniform_loss_functional_enrichment.py"
SPEC = importlib.util.spec_from_file_location("_uniform_enrichment_base", BASE_SCRIPT)
BASE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(BASE)


class TreeEnrichmentError(ValueError):
    """Raised when a frozen foreground or annotation input is inconsistent."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-protein", required=True, type=Path)
    parser.add_argument("--foreground-genes", required=True, type=Path)
    parser.add_argument("--foreground-metadata", required=True, type=Path)
    parser.add_argument("--resolved-background-genes", type=Path)
    parser.add_argument(
        "--resolved-background-scope",
        default="all_13_lineages_resolved_article_method",
    )
    parser.add_argument("--foreground-background-genes", type=Path)
    parser.add_argument(
        "--analysis-profile",
        choices=("tree", "assembly_unit", "scaffold"),
        default="tree",
    )
    parser.add_argument("--emapper-annotations", required=True, type=Path)
    parser.add_argument("--go-obo", required=True, type=Path)
    parser.add_argument("--ko-descriptions", required=True, type=Path)
    parser.add_argument("--expected-reference-genes", type=int, default=35547)
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
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.name.endswith(".gz") else path.open(encoding="utf-8", newline="")


def read_table(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields or len(fields) != len(set(fields)):
            raise TreeEnrichmentError(f"{path.name}: invalid TSV header")
        rows = [dict(row) for row in reader]
    if not rows:
        raise TreeEnrichmentError(f"{path.name}: no data rows")
    return rows, fields


def require(fields: Iterable[str], needed: set[str], label: str) -> None:
    missing = sorted(needed - set(fields))
    if missing:
        raise TreeEnrichmentError(f"{label} missing columns: {', '.join(missing)}")


def read_ids(path: Path) -> set[str]:
    identifiers = {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}
    if not identifiers:
        raise TreeEnrichmentError(f"{path.name}: empty gene-ID set")
    return identifiers


def read_backgrounds(path: Path) -> dict[str, set[str]]:
    rows, fields = read_table(path)
    require(
        fields,
        {"background_scope", "reference_gene_id"},
        path.name,
    )
    backgrounds: dict[str, set[str]] = defaultdict(set)
    seen: set[tuple[str, str]] = set()
    for line_number, row in enumerate(rows, 2):
        key = (
            row["background_scope"].strip(),
            row["reference_gene_id"].strip(),
        )
        if not all(key) or key in seen:
            raise TreeEnrichmentError(
                f"{path.name}:{line_number}: empty or duplicate background row"
            )
        seen.add(key)
        backgrounds[key[0]].add(key[1])
    return dict(backgrounds)


def read_foregrounds(
    genes_path: Path,
    metadata_path: Path,
) -> tuple[dict[str, set[str]], dict[str, dict[str, str]]]:
    gene_rows, gene_fields = read_table(genes_path)
    require(gene_fields, {"foreground_id", "reference_gene_id"}, genes_path.name)
    foregrounds: dict[str, set[str]] = defaultdict(set)
    seen: set[tuple[str, str]] = set()
    for line_number, row in enumerate(gene_rows, 2):
        key = (row["foreground_id"].strip(), row["reference_gene_id"].strip())
        if not all(key) or key in seen:
            raise TreeEnrichmentError(f"{genes_path.name}:{line_number}: empty or duplicate membership")
        seen.add(key)
        foregrounds[key[0]].add(key[1])
    metadata_rows, metadata_fields = read_table(metadata_path)
    require(
        metadata_fields,
        {
            "foreground_id", "analysis_scope", "background_scope",
            "foreground_gene_count", "branch_id", "descendant_lineage_count",
            "descendant_lineages",
        },
        metadata_path.name,
    )
    metadata: dict[str, dict[str, str]] = {}
    for line_number, row in enumerate(metadata_rows, 2):
        foreground_id = row["foreground_id"].strip()
        if not foreground_id or foreground_id in metadata:
            raise TreeEnrichmentError(f"{metadata_path.name}:{line_number}: empty or duplicate foreground")
        try:
            expected = int(row["foreground_gene_count"])
        except ValueError as error:
            raise TreeEnrichmentError(f"{metadata_path.name}:{line_number}: invalid foreground count") from error
        if expected != len(foregrounds.get(foreground_id, set())):
            raise TreeEnrichmentError(f"{metadata_path.name}:{line_number}: membership count mismatch")
        metadata[foreground_id] = row
    if set(metadata) != set(foregrounds):
        raise TreeEnrichmentError("foreground metadata and membership universes differ")
    return dict(foregrounds), metadata


def write_tsv(path: Path, fields: list[str], rows: Iterable[Mapping[str, object]], *, gz: bool = False) -> None:
    opener = gzip.open if gz else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise TreeEnrichmentError(f"output directory already exists: {args.output_dir}")
    if args.minimum_background_count < 1 or args.minimum_significant_study_count < 1:
        raise TreeEnrichmentError("minimum counts must be positive")
    if not 0 < args.fdr_cutoff <= 1:
        raise TreeEnrichmentError("FDR cutoff must be in (0, 1]")
    inputs = {
        "reference_protein": args.reference_protein,
        "foreground_genes": args.foreground_genes,
        "foreground_metadata": args.foreground_metadata,
        "emapper_annotations": args.emapper_annotations,
        "go_obo": args.go_obo,
        "ko_descriptions": args.ko_descriptions,
    }
    if args.resolved_background_genes is not None:
        inputs["resolved_background_genes"] = args.resolved_background_genes
    if args.foreground_background_genes is not None:
        inputs["foreground_background_genes"] = (
            args.foreground_background_genes
        )
    for role, path in inputs.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise TreeEnrichmentError(f"missing or empty {role}: {path}")
    hashes = {role: sha256(path) for role, path in inputs.items()}
    expected_hashes = {
        "reference_protein": args.expected_reference_protein_sha256,
        "emapper_annotations": args.expected_emapper_sha256,
        "go_obo": args.expected_go_obo_sha256,
        "ko_descriptions": args.expected_ko_descriptions_sha256,
    }
    for role, expected in expected_hashes.items():
        if expected and hashes[role] != expected.lower():
            raise TreeEnrichmentError(
                f"{role} SHA-256 mismatch: expected {expected.lower()}, observed {hashes[role]}"
            )

    reference_order = BASE.read_fasta_ids(args.reference_protein)
    if len(reference_order) != args.expected_reference_genes:
        raise TreeEnrichmentError(
            f"reference protein has {len(reference_order)} genes; expected {args.expected_reference_genes}"
        )
    reference = set(reference_order)
    if (
        args.analysis_profile in {"tree", "scaffold"}
        and args.resolved_background_genes is None
    ):
        raise TreeEnrichmentError(
            "tree/scaffold profile requires --resolved-background-genes"
        )
    resolved_background = (
        read_ids(args.resolved_background_genes)
        if args.resolved_background_genes is not None
        else set()
    )
    if not args.resolved_background_scope.strip():
        raise TreeEnrichmentError("resolved background scope is empty")
    if not resolved_background.issubset(reference):
        raise TreeEnrichmentError(
            "resolved tree background contains non-reference genes"
        )
    foregrounds, metadata = read_foregrounds(args.foreground_genes, args.foreground_metadata)
    if not set().union(*foregrounds.values()).issubset(reference):
        raise TreeEnrichmentError("foreground contains non-reference genes")
    allowed_backgrounds = {"all_reference_genes": reference}
    if args.resolved_background_genes is not None:
        allowed_backgrounds[
            args.resolved_background_scope
        ] = resolved_background
    if args.foreground_background_genes is not None:
        unit_backgrounds = read_backgrounds(
            args.foreground_background_genes
        )
        if any(
            not genes.issubset(reference)
            for genes in unit_backgrounds.values()
        ):
            raise TreeEnrichmentError(
                "custom foreground background contains non-reference genes"
            )
        overlap = set(allowed_backgrounds).intersection(unit_backgrounds)
        if overlap:
            raise TreeEnrichmentError(
                f"duplicate background scopes: {sorted(overlap)}"
            )
        allowed_backgrounds.update(unit_backgrounds)
    elif args.analysis_profile == "assembly_unit":
        raise TreeEnrichmentError(
            "assembly-unit profile requires --foreground-background-genes"
        )
    for foreground_id, genes in foregrounds.items():
        scope = metadata[foreground_id]["background_scope"]
        if scope not in allowed_backgrounds:
            raise TreeEnrichmentError(f"{foreground_id}: unsupported background scope {scope!r}")
        if not genes.issubset(allowed_backgrounds[scope]):
            raise TreeEnrichmentError(f"{foreground_id}: foreground is not a subset of its background")

    associations, annotation_qc, emapper_metadata = BASE.read_emapper(args.emapper_annotations)
    if not set(associations).issubset(reference):
        raise TreeEnrichmentError("eggNOG annotations contain non-reference genes")
    go_terms, go_version = BASE.read_obo(args.go_obo)
    ko_descriptions = BASE.read_ko_descriptions(args.ko_descriptions)
    unknown_go = sorted(set().union(*(row["GO"] for row in associations.values())) - set(go_terms))
    if unknown_go:
        assignments = 0
        for values in associations.values():
            assignments += len(values["GO"].intersection(unknown_go))
            values["GO"].difference_update(unknown_go)
        annotation_qc["unresolved_go_terms_excluded"] = len(unknown_go)
        annotation_qc["unresolved_gene_go_assignments_excluded"] = assignments

    try:
        import scipy
        from scipy.stats import hypergeom
    except ImportError as error:
        raise TreeEnrichmentError("scipy is required for enrichment") from error

    all_fields = [
        "foreground_id", "analysis_scope", "background_scope", "branch_id",
        "descendant_lineage_count", "descendant_lineages", "ontology", "term_id",
        "term_name", "go_namespace", "study_count", "study_size", "background_count",
        "background_size", "p_overrepresentation", "p_fdr_bh", "fold_enrichment",
        "significant_fdr", "study_gene_ids",
    ]
    summary_fields = [
        "foreground_id", "analysis_scope", "background_scope", "branch_id",
        "descendant_lineage_count", "descendant_lineages", "ontology",
        "requested_gene_count", "annotated_study_gene_count",
        "annotated_background_gene_count", "tested_term_count", "annotation_coverage",
    ]
    all_rows: list[dict[str, object]] = []
    significant: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    membership_rows = [
        {"foreground_id": foreground_id, "reference_gene_id": gene}
        for foreground_id in sorted(foregrounds)
        for gene in sorted(foregrounds[foreground_id])
    ]
    for foreground_id in sorted(foregrounds):
        meta = metadata[foreground_id]
        meta_out = {field: meta.get(field, "") for field in [
            "analysis_scope", "background_scope", "branch_id",
            "descendant_lineage_count", "descendant_lineages",
        ]}
        requested = foregrounds[foreground_id]
        background_universe = allowed_backgrounds[meta["background_scope"]]
        for ontology in BASE.ONTOLOGIES:
            background = {
                gene for gene in background_universe
                if gene in associations and associations[gene][ontology]
            }
            study = requested & background
            term_to_genes: dict[str, set[str]] = defaultdict(set)
            for gene in background:
                for term in associations[gene][ontology]:
                    term_to_genes[term].add(gene)
            tested = sorted(
                term for term, genes in term_to_genes.items()
                if len(genes) >= args.minimum_background_count
            )
            summary_rows.append(
                {
                    "foreground_id": foreground_id,
                    **meta_out,
                    "ontology": ontology,
                    "requested_gene_count": len(requested),
                    "annotated_study_gene_count": len(study),
                    "annotated_background_gene_count": len(background),
                    "tested_term_count": len(tested),
                    "annotation_coverage": len(study) / len(requested) if requested else "",
                }
            )
            raw_rows: list[dict[str, object]] = []
            p_values: list[float] = []
            for term in tested:
                term_background = term_to_genes[term]
                hits = study & term_background
                p_value = float(
                    hypergeom.sf(len(hits) - 1, len(background), len(term_background), len(study))
                ) if study else 1.0
                fold = (
                    (len(hits) / len(study)) / (len(term_background) / len(background))
                    if study and hits else 0.0
                )
                if ontology == "GO":
                    term_name, namespace = go_terms[term]
                elif ontology == "KEGG_KO":
                    term_name, namespace = ko_descriptions.get(term, ""), ""
                else:
                    term_name, namespace = term, ""
                raw_rows.append(
                    {
                        "foreground_id": foreground_id,
                        **meta_out,
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
            for row, adjusted in zip(raw_rows, BASE.bh_adjust(p_values)):
                row["p_fdr_bh"] = adjusted
                is_significant = (
                    adjusted <= args.fdr_cutoff
                    and int(row["study_count"]) >= args.minimum_significant_study_count
                    and float(row["fold_enrichment"]) > 1
                )
                row["significant_fdr"] = str(is_significant).lower()
                all_rows.append(row)
                if is_significant:
                    significant.append(row)

    all_rows.sort(key=lambda row: (str(row["foreground_id"]), str(row["ontology"]), float(row["p_fdr_bh"]), str(row["term_id"])))
    significant.sort(key=lambda row: (str(row["foreground_id"]), str(row["ontology"]), float(row["p_fdr_bh"]), str(row["term_id"])))
    parent = args.output_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.tmp.", dir=parent))
    try:
        write_tsv(temporary / "foreground_summary.tsv", summary_fields, summary_rows)
        write_tsv(temporary / "foreground_gene_ids.tsv.gz", ["foreground_id", "reference_gene_id"], membership_rows, gz=True)
        write_tsv(temporary / "enrichment_all_terms.tsv.gz", all_fields, all_rows, gz=True)
        write_tsv(temporary / "enrichment_significant.tsv", all_fields, significant)
        outputs = [
            temporary / "foreground_summary.tsv",
            temporary / "foreground_gene_ids.tsv.gz",
            temporary / "enrichment_all_terms.tsv.gz",
            temporary / "enrichment_significant.tsv",
        ]
        significant_counts: dict[str, int] = defaultdict(int)
        for row in significant:
            significant_counts[f"{row['foreground_id']}|{row['ontology']}"] += 1
        is_unit_profile = args.analysis_profile == "assembly_unit"
        is_scaffold_profile = args.analysis_profile == "scaffold"
        if is_unit_profile:
            status = "PASS_UNIT_ARTICLE_METHOD_GO_KEGG"
        elif is_scaffold_profile:
            status = "PASS_UNIT_SCAFFOLD_GO_KEGG"
        else:
            status = "PASS_TREE_AWARE_MANUSCRIPT_GO_KEGG"
        definitions = {
            "loss_classification": "article-method decayed + deleted",
            "test": (
                "one-sided hypergeometric over-representation; "
                "BH within foreground and ontology"
            ),
            "go": (
                "direct frozen eggNOG-mapper GO terms; roots and unresolved "
                "OBO IDs excluded; no ancestor propagation"
            ),
            "kegg": "offline frozen eggNOG-mapper KEGG KO and pathway fields",
        }
        if is_unit_profile:
            definitions.update(
                {
                    "foregrounds": (
                        "23 independent assembly-unit decayed + deleted sets; "
                        "haplotypes and subgenomes are not aggregated"
                    ),
                    "backgrounds": (
                        "matching per-unit retained + decayed + deleted genes; "
                        "not_called_loss excluded"
                    ),
                }
            )
        elif is_scaffold_profile:
            definitions["scaffold_foregrounds"] = (
                "maximal decayed-plus-deleted event genes placed on the "
                "topology-only 23-unit scaffold; assembly units remain "
                "separate and the scaffold is not a 23-species phylogeny"
            )
            definitions["scaffold_background"] = (
                "reference genes with resolved article-method states in all "
                "23 assembly units"
            )
        else:
            definitions["tree_foregrounds"] = (
                "complete biological-lineage losses placed on the exact "
                "matching Actinidia topology; partial and unknown lineage "
                "states are not ancestral events"
            )
        manifest = {
            "schema_version": "1.0",
            "status": status,
            "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "analysis_profile": args.analysis_profile,
            "definitions": definitions,
            "counts": {
                "reference_genes": len(reference),
                "foregrounds": len(foregrounds),
                "foreground_memberships": sum(
                    len(genes) for genes in foregrounds.values()
                ),
                "significant_terms": len(significant),
                **(
                    (
                        {
                            "resolved_scaffold_background_genes": len(
                                resolved_background
                            )
                        }
                        if is_scaffold_profile
                        else {
                            "resolved_tree_background_genes": len(
                                resolved_background
                            )
                        }
                    )
                    if not is_unit_profile else {}
                ),
            },
            "parameters": {
                "minimum_background_count": args.minimum_background_count,
                "minimum_significant_study_count": args.minimum_significant_study_count,
                "fdr_cutoff": args.fdr_cutoff,
            },
            "annotation_qc": annotation_qc,
            "go_obo_data_version": go_version,
            "emapper_metadata": [line for line in emapper_metadata if line and "/" not in line and "--" not in line],
            "significant_term_counts": dict(sorted(significant_counts.items())),
            "inputs": [
                {"role": role, "basename": path.name, "sha256": hashes[role]}
                for role, path in inputs.items()
            ],
            "outputs": [
                {"basename": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
                for path in outputs
            ],
            "software": {"python": os.sys.version.split()[0], "scipy": scipy.__version__},
        }
        if is_unit_profile:
            summary_name = "unit_functional_enrichment_summary.json"
        elif is_scaffold_profile:
            summary_name = "scaffold_functional_enrichment_summary.json"
        else:
            summary_name = "tree_functional_enrichment_summary.json"
        (temporary / summary_name).write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, args.output_dir)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (TreeEnrichmentError, BASE.EnrichmentError, OSError, csv.Error, UnicodeError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    print(f"PASS\t{args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
