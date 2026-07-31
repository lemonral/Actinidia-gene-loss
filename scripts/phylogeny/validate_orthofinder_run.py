#!/usr/bin/env python3
"""Validate an OrthoFinder 3 result against its exact proteome bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path

from phylo_io import DataError, fasta_source_stem, load_sequence_id_map, read_fasta


INTERNAL_NODE = re.compile(r"\)(N\d+)(?=[:),;])")
ANSI_SGR = re.compile(r"\x1b\[[0-9;]*m")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_genes(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def read_statistics(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line:
                break
            fields = line.split("\t")
            if len(fields) != 2 or not fields[0] or fields[0] in values:
                raise DataError(f"{path}: malformed or duplicate summary statistic {line!r}")
            values[fields[0]] = fields[1]
    return values


def read_orthogroup_table(
    path: Path,
    expected_species: list[str],
    gene_to_species: dict[str, str],
) -> tuple[dict[str, dict[str, list[str]]], set[str]]:
    groups: dict[str, dict[str, list[str]]] = {}
    genes: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Orthogroup", *expected_species]:
            raise DataError(f"{path}: species header does not exactly match SpeciesIDs.txt")
        for line_number, row in enumerate(reader, start=2):
            orthogroup = (row.get("Orthogroup") or "").strip()
            if not re.fullmatch(r"OG\d{7}", orthogroup) or orthogroup in groups:
                raise DataError(f"{path}:{line_number}: invalid or duplicate orthogroup {orthogroup!r}")
            per_species: dict[str, list[str]] = {}
            for species in expected_species:
                values = split_genes(row.get(species) or "")
                if len(values) != len(set(values)):
                    raise DataError(f"{path}:{line_number}: duplicate gene within {species}")
                for gene in values:
                    observed_species = gene_to_species.get(gene)
                    if observed_species != species:
                        raise DataError(
                            f"{path}:{line_number}: {gene} belongs to {observed_species!r}, not {species!r}"
                        )
                    if gene in genes:
                        raise DataError(f"{path}:{line_number}: gene occurs in more than one row: {gene}")
                    genes.add(gene)
                per_species[species] = values
            groups[orthogroup] = per_species
    return groups, genes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--input-bundle", required=True, type=Path)
    parser.add_argument("--controller-log", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-species", type=int, default=17)
    parser.add_argument("--expected-single-copy", type=int, default=479)
    args = parser.parse_args()

    temporary: Path | None = None
    try:
        if args.output_dir.exists():
            raise DataError(f"refusing to overwrite existing output directory: {args.output_dir}")
        binding_path = args.input_bundle / "input_binding.json"
        manifest_path = args.input_bundle / "orthofinder_input_manifest.tsv"
        binding = json.loads(binding_path.read_text(encoding="utf-8"))
        if binding.get("status") != "PASS":
            raise DataError("OrthoFinder input bundle is not PASS")
        with manifest_path.open("r", encoding="utf-8", newline="") as handle:
            manifest = list(csv.DictReader(handle, delimiter="\t"))
        if len(manifest) != args.expected_species:
            raise DataError(f"input manifest has {len(manifest)} species; expected {args.expected_species}")
        manifest_by_stem: dict[str, dict[str, str]] = {}
        for row in manifest:
            link = args.input_bundle / row["proteome_link"]
            stem = fasta_source_stem(link.name)
            if stem in manifest_by_stem or not link.is_file():
                raise DataError(f"invalid or duplicate input proteome stem {stem!r}")
            if sha256_file(link) != row["protein_sha256"]:
                raise DataError(f"input proteome checksum mismatch for {stem}")
            manifest_by_stem[stem] = row

        log_text = ANSI_SGR.sub(
            "", args.controller_log.read_text(encoding="utf-8", errors="strict")
        )
        if "OrthoFinder finished" not in log_text or "Results directory:" not in log_text:
            raise DataError("controller log does not contain a clean OrthoFinder completion")

        working = args.results / "WorkingDirectory"
        sequence_ids = working / "SequenceIDs.txt"
        species_ids = working / "SpeciesIDs.txt"
        sequence_to_source = load_sequence_id_map(sequence_ids, species_ids=species_ids)
        gene_to_species = {
            gene: fasta_source_stem(source) for gene, source in sequence_to_source.items()
        }
        if set(gene_to_species.values()) != set(manifest_by_stem):
            raise DataError("SpeciesIDs.txt does not close to the exact input proteome set")
        observed_counts = Counter(gene_to_species.values())
        for stem, row in manifest_by_stem.items():
            if observed_counts[stem] != int(row["protein_records"]):
                raise DataError(
                    f"{stem}: SequenceIDs count {observed_counts[stem]} != {row['protein_records']}"
                )

        statistics_path = args.results / "Comparative_Genomics_Statistics" / "Statistics_Overall.tsv"
        stats = read_statistics(statistics_path)
        species_count = int(stats["Number of species"])
        total_genes = int(stats["Number of genes"])
        assigned_count = int(stats["Number of genes in orthogroups"])
        unassigned_count = int(stats["Number of unassigned genes"])
        orthogroup_count = int(stats["Number of orthogroups"])
        single_copy_count = int(stats["Number of single-copy orthogroups"])
        if species_count != args.expected_species or total_genes != len(gene_to_species):
            raise DataError("statistics species/gene counts do not close to SpeciesIDs/SequenceIDs")
        if assigned_count + unassigned_count != total_genes:
            raise DataError("assigned plus unassigned genes do not equal total genes")
        if single_copy_count != args.expected_single_copy:
            raise DataError(f"statistics report {single_copy_count} SCOs; expected {args.expected_single_copy}")

        expected_species = list(manifest_by_stem)
        species_file_order: list[str] = []
        with species_ids.open("r", encoding="utf-8") as handle:
            for raw in handle:
                if ": " in raw:
                    species_file_order.append(fasta_source_stem(raw.rstrip().split(": ", 1)[1]))
        expected_species = species_file_order
        orthogroups_path = args.results / "Orthogroups" / "Orthogroups.tsv"
        groups, assigned = read_orthogroup_table(
            orthogroups_path, expected_species, gene_to_species
        )
        if len(groups) != orthogroup_count or len(assigned) != assigned_count:
            raise DataError("orthogroup table counts do not match Statistics_Overall.tsv")
        unassigned_path = args.results / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv"
        unassigned_groups, unassigned = read_orthogroup_table(
            unassigned_path, expected_species, gene_to_species
        )
        if len(unassigned) != unassigned_count or assigned.intersection(unassigned):
            raise DataError("unassigned-gene table count/disjointness failed")
        if assigned.union(unassigned) != set(gene_to_species):
            raise DataError("assigned and unassigned tables do not close to all input genes")

        sco_path = args.results / "Orthogroups" / "Orthogroups_SingleCopyOrthologues.txt"
        sco = [line.strip() for line in sco_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        if len(sco) != args.expected_single_copy or len(set(sco)) != len(sco):
            raise DataError("single-copy orthogroup list count/uniqueness failed")
        sco_directory = args.results / "Single_Copy_Orthologue_Sequences"
        for orthogroup in sco:
            if orthogroup not in groups:
                raise DataError(f"single-copy orthogroup is absent from Orthogroups.tsv: {orthogroup}")
            expected_genes = {values[0] for values in groups[orthogroup].values() if len(values) == 1}
            if len(expected_genes) != args.expected_species or any(
                len(values) != 1 for values in groups[orthogroup].values()
            ):
                raise DataError(f"{orthogroup}: does not contain exactly one gene per species")
            records = read_fasta(sco_directory / f"{orthogroup}.fa")
            if set(records) != expected_genes:
                raise DataError(f"{orthogroup}: SCO FASTA does not match Orthogroups.tsv")

        rooted_tree_path = args.results / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt"
        internal_nodes = set(INTERNAL_NODE.findall(rooted_tree_path.read_text(encoding="utf-8")))
        hog_directory = args.results / "Phylogenetic_Hierarchical_Orthogroups"
        hog_paths = sorted(hog_directory.glob("N*.tsv"))
        if {path.stem for path in hog_paths} != internal_nodes - {"N0"}:
            raise DataError("HOG node files do not close to rooted internal node labels")
        for path in hog_paths:
            with path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != [
                    "HOG", "OG", "Gene Tree Parent Clade", *expected_species
                ]:
                    raise DataError(f"{path}: HOG header closure failed")
                seen_hogs: set[str] = set()
                rows = 0
                for line_number, row in enumerate(reader, start=2):
                    rows += 1
                    hog = (row.get("HOG") or "").strip()
                    orthogroup = (row.get("OG") or "").strip()
                    if not hog.startswith(f"{path.stem}.HOG") or hog in seen_hogs:
                        raise DataError(f"{path}:{line_number}: invalid or duplicate HOG ID")
                    if orthogroup not in groups:
                        raise DataError(f"{path}:{line_number}: unknown orthogroup {orthogroup}")
                    seen_hogs.add(hog)
                    for species in expected_species:
                        for gene in split_genes(row.get(species) or ""):
                            if gene_to_species.get(gene) != species:
                                raise DataError(f"{path}:{line_number}: HOG gene/species mismatch")
                if rows == 0:
                    raise DataError(f"{path}: empty HOG table")

        key_paths = [
            binding_path,
            manifest_path,
            args.controller_log,
            sequence_ids,
            species_ids,
            statistics_path,
            orthogroups_path,
            unassigned_path,
            sco_path,
            rooted_tree_path,
            *hog_paths,
        ]
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        temporary = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=args.output_dir.parent))
        checksum_manifest = temporary / "key_output_checksums.tsv"
        with checksum_manifest.open("w", encoding="utf-8", newline="") as handle:
            handle.write("basename\tbytes\tsha256\n")
            for path in key_paths:
                handle.write(f"{path.name}\t{path.stat().st_size}\t{sha256_file(path)}\n")
        result = {
            "schema_version": 1,
            "status": "PASS",
            "workflow": "orthofinder3_exact_run_validation",
            "species": species_count,
            "genes": total_genes,
            "assigned_genes": assigned_count,
            "unassigned_genes": unassigned_count,
            "orthogroups": orthogroup_count,
            "all_species_orthogroups": int(stats["Number of orthogroups with all species present"]),
            "single_copy_orthogroups": single_copy_count,
            "hog_node_tables": len(hog_paths),
            "input_binding_sha256": sha256_file(binding_path),
            "key_output_checksums_sha256": sha256_file(checksum_manifest),
        }
        (temporary / "validation.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(temporary, args.output_dir)
        temporary = None
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except (DataError, OSError, csv.Error, json.JSONDecodeError, KeyError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        if temporary is not None:
            shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
