#!/usr/bin/env python3
"""Manifest-driven, read-only quality checks for phylogeny and CAFE outputs.

The program deliberately does *not* launch MAFFT, RAxML-NG or CAFE.
It validates the hand-off between those tools before a new run is accepted.
All source files are opened read-only; a report is printed to stdout unless an
explicit ``--report`` path is supplied.

Examples
--------
python scripts/phylogeny/phylo_qc.py metadata \
  --manifest config/phylogeny/taxa.tsv

python scripts/phylogeny/phylo_qc.py tree \
  --tree results/phylogeny/<run_id>/species_tree/raxml.support \
  --manifest config/phylogeny/terminals.selected.tsv

"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

from phylo_io import (
    DataError,
    bool_value,
    load_samples,
    load_sequence_id_map,
    load_terminal_manifest,
    read_fasta,
    read_tsv,
    report_tsv,
    resolve_record_samples,
    source_label_to_sample,
    validate_equal_lengths,
)


LEAF_RE = re.compile(r"(?<=[(,])\s*('(?:[^']|'')*'|[^\s():,;\[\]]+)")


def result(rows: list[dict[str, object]], check: str, severity: str, status: str, detail: str, **extra: object) -> None:
    row: dict[str, object] = {
        "check": check,
        "severity": severity,
        "status": status,
        "detail": detail,
    }
    row.update(extra)
    rows.append(row)


def finish(rows: list[dict[str, object]], report: str | None) -> int:
    print(report_tsv(rows, report), end="")
    return 2 if any(row["severity"] == "ERROR" for row in rows) else 0


def load_aliases(path: str | Path | None) -> list[dict[str, str]]:
    """Load an optional context-specific alias table.

    Both ``terminal_id`` and the previous ``sample_id`` target-column name are
    accepted.  Alias files are optional for newly generated canonical labels.
    """

    if path is None:
        return []
    aliases = read_tsv(path, ("analysis_context", "legacy_label", "alias_type"))
    for row in aliases:
        terminal_id = row.get("terminal_id", row.get("sample_id", ""))
        if not terminal_id:
            raise DataError(
                f"{path}:{row['__line__']}: alias row needs terminal_id or sample_id"
            )
        row["sample_id"] = terminal_id
        row["terminal_id"] = terminal_id
    return aliases


def alias_mapping(
    samples: list[dict[str, str]], aliases: list[dict[str, str]], context: str
) -> tuple[dict[str, str], list[str]]:
    """Map labels valid in one context to stable sample IDs.

    Canonical tree labels and stable terminal IDs are always accepted.  Context aliases may be added
    only once; conflicting aliases are an unsafe configuration error.
    """

    mapping: dict[str, str] = {}
    errors: list[str] = []
    for row in samples:
        for label in (row["canonical_tree_label"], row["sample_id"]):
            prior = mapping.get(label)
            if prior is not None and prior != row["sample_id"]:
                errors.append(
                    f"canonical label {label!r} maps to both {prior!r} and {row['sample_id']!r}"
                )
            else:
                mapping[label] = row["sample_id"]
    valid_sample_ids = {row["sample_id"] for row in samples}
    for row in aliases:
        if row["analysis_context"] not in {context, "all"}:
            continue
        label = row["legacy_label"]
        sample_id = row["sample_id"]
        if sample_id not in valid_sample_ids:
            errors.append(
                f"aliases.tsv line {row['__line__']}: unknown sample_id {sample_id!r}"
            )
            continue
        prior = mapping.get(label)
        if prior is not None and prior != sample_id:
            errors.append(
                f"alias {label!r} maps to both {prior!r} and {sample_id!r} in context {context}"
            )
        else:
            mapping[label] = sample_id
    return mapping, errors


def parse_newick_leaves(path: str | Path) -> list[str]:
    """Extract terminal labels from a conventional Newick/NHX file.

    This is intentionally a restricted parser: it reliably handles the tree
    files in this project, including branch lengths and NHX annotations, but it
    is not intended to interpret topology or arbitrary comments.  Labels with
    whitespace must be single-quoted Newick labels.
    """

    text = Path(path).read_text(encoding="utf-8").strip()
    if not text.endswith(";"):
        raise DataError(f"{path}: Newick tree must end with ';'")
    # Drop square-bracket annotations before applying the terminal-label regex.
    cleaned = re.sub(r"\[[^\]]*\]", "", text)
    leaves: list[str] = []
    for match in LEAF_RE.finditer(cleaned):
        label = match.group(1).strip()
        if label.startswith("'") and label.endswith("'"):
            label = label[1:-1].replace("''", "'")
        if label:
            leaves.append(label)
    if not leaves:
        raise DataError(f"{path}: no terminal labels parsed")
    return leaves


def command_metadata(args: argparse.Namespace) -> int:
    rows: list[dict[str, object]] = []
    try:
        manifest = load_terminal_manifest(args.manifest)
        samples = manifest.selected
        aliases = load_aliases(args.aliases)
    except DataError as error:
        result(rows, "metadata", "ERROR", "FAIL", str(error))
        return finish(rows, args.report)

    result(
        rows,
        "manifest_schema",
        "INFO",
        "PASS",
        f"loaded {manifest.schema} with {len(manifest.all_rows)} registered rows",
        manifest=str(args.manifest),
    )
    if not samples:
        result(
            rows,
            "terminal_selection",
            "WARN",
            "ACTION_REQUIRED",
            "no row is currently selected with include_species_tree=true",
        )
    else:
        result(
            rows,
            "terminal_selection",
            "INFO",
            "PASS",
            f"{len(samples)} terminals are explicitly selected",
        )
    if manifest.candidates:
        result(
            rows,
            "candidate_taxa",
            "WARN",
            "ACTION_REQUIRED",
            f"{len(manifest.candidates)} candidate rows remain outside the selected tree: "
            + ", ".join(row["sample_id"] for row in manifest.candidates),
        )
    if manifest.excluded:
        result(
            rows,
            "excluded_taxa",
            "INFO",
            "PASS",
            f"{len(manifest.excluded)} rows are explicitly excluded",
        )

    expected_terminals = args.expected_terminals
    terminal_mismatch = expected_terminals is not None and len(samples) != expected_terminals
    result(
        rows,
        "terminal_count",
        "ERROR" if terminal_mismatch else "INFO",
        "FAIL" if terminal_mismatch else "PASS",
        (
            f"found {len(samples)} selected terminals; explicit expectation is {expected_terminals}"
            if expected_terminals is not None
            else f"found {len(samples)} selected terminals; count is derived from the manifest"
        ),
    )
    biological_species = {row["biological_species_id"] for row in samples}
    expected_species = args.expected_species
    species_mismatch = expected_species is not None and len(biological_species) != expected_species
    result(
        rows,
        "biological_species_count",
        "ERROR" if species_mismatch else "INFO",
        "FAIL" if species_mismatch else "PASS",
        (
            f"found {len(biological_species)} selected biological-species groups; "
            f"explicit expectation is {expected_species}"
            if expected_species is not None
            else f"found {len(biological_species)} selected biological-species groups"
        ),
    )
    species_counts = Counter(row["biological_species_id"] for row in samples)
    duplicated_species = {species: count for species, count in species_counts.items() if count > 1}
    result(
        rows,
        "multi_terminal_species",
        "INFO",
        "PASS",
        "biological-species groups with multiple technical terminals: "
        + (
            ", ".join(f"{species} ({duplicated_species[species]})" for species in sorted(duplicated_species))
            or "none"
        ),
    )
    try:
        roots = [row["sample_id"] for row in samples if bool_value(row["is_root_outgroup"])]
        root_error = None
    except DataError as error:
        roots, root_error = [], str(error)
    root_failure = bool(root_error) or len(roots) > 1 or (args.require_root and len(roots) != 1)
    result(
        rows,
        "rooting_outgroup",
        "ERROR" if root_failure else ("WARN" if not roots else "INFO"),
        "FAIL" if root_failure else ("ACTION_REQUIRED" if not roots else "PASS"),
        root_error
        or (
            f"configured root outgroup: {', '.join(roots)}"
            if roots
            else "no selected terminal is marked is_root_outgroup=true"
        ),
    )

    pending_identity = []
    for row in manifest.all_rows:
        identity_status = row.get("identity_status", "").lower()
        needs_confirmation = not identity_status.startswith("confirmed") and any(
            token in identity_status
            for token in ("confirm", "pending", "revalid", "select", "audit", "conflict", "unresolved")
        )
        if needs_confirmation:
            pending_identity.append(row["sample_id"])
    if pending_identity:
        result(
            rows,
            "taxon_identity",
            "WARN",
            "ACTION_REQUIRED",
            "author confirmation is still required for: " + ", ".join(pending_identity),
        )

    valid_sample_ids = {row["sample_id"] for row in samples}
    seen_aliases: dict[tuple[str, str], tuple[str, str]] = {}
    for alias in aliases:
        if alias["sample_id"] not in valid_sample_ids:
            result(
                rows,
                "alias_reference",
                "ERROR",
                "FAIL",
                f"line {alias['__line__']}: unknown sample_id {alias['sample_id']!r}",
            )
        key = (alias["analysis_context"], alias["legacy_label"])
        previous = seen_aliases.get(key)
        if previous is not None and previous[0] != alias["sample_id"]:
            result(
                rows,
                "alias_collision",
                "ERROR",
                "FAIL",
                f"{key!r} maps to both {previous[0]!r} (line {previous[1]}) and {alias['sample_id']!r}",
            )
        seen_aliases[key] = (alias["sample_id"], alias["__line__"])
    result(
        rows,
        "alias_table",
        "INFO",
        "PASS",
        f"validated {len(aliases)} context-specific aliases"
        if aliases
        else "no alias table supplied; only canonical labels and stable terminal IDs are accepted",
    )
    return finish(rows, args.report)


def command_tree(args: argparse.Namespace) -> int:
    rows: list[dict[str, object]] = []
    try:
        samples = load_samples(args.manifest)
        aliases = load_aliases(args.aliases)
        mapping, mapping_errors = alias_mapping(samples, aliases, args.context)
        leaves = parse_newick_leaves(args.tree)
    except DataError as error:
        result(rows, "tree", "ERROR", "FAIL", str(error))
        return finish(rows, args.report)
    for error in mapping_errors:
        result(rows, "alias_mapping", "ERROR", "FAIL", error)

    resolved: list[str] = []
    aliases_used: list[str] = []
    canonical_labels = {
        label for row in samples for label in (row["canonical_tree_label"], row["sample_id"])
    }
    unknown: list[str] = []
    for leaf in leaves:
        sample_id = mapping.get(leaf)
        if sample_id is None:
            unknown.append(leaf)
        else:
            resolved.append(sample_id)
            if leaf not in canonical_labels:
                aliases_used.append(leaf)
    expected_terminals = args.expected_terminals if args.expected_terminals is not None else len(samples)
    result(
        rows,
        "tree_tip_count",
        "ERROR" if len(leaves) != expected_terminals else "INFO",
        "FAIL" if len(leaves) != expected_terminals else "PASS",
        f"parsed {len(leaves)} terminal labels; expected {expected_terminals} "
        + ("from the manifest" if args.expected_terminals is None else "from the explicit override"),
        tree=str(args.tree),
        context=args.context,
    )
    if unknown:
        result(
            rows,
            "tree_unknown_labels",
            "ERROR",
            "FAIL",
            "unmapped labels: " + ", ".join(sorted(set(unknown))),
        )
    duplicates = sorted(sample_id for sample_id, count in Counter(resolved).items() if count > 1)
    if duplicates:
        result(
            rows,
            "tree_duplicate_samples",
            "ERROR",
            "FAIL",
            "more than one terminal resolved to: " + ", ".join(duplicates),
        )
    expected = {row["sample_id"] for row in samples}
    missing = sorted(expected.difference(resolved))
    if missing:
        result(
            rows,
            "tree_missing_samples",
            "ERROR",
            "FAIL",
            "configured samples absent from tree: " + ", ".join(missing),
        )
    if aliases_used:
        result(
            rows,
            "tree_aliases_used",
            "WARN",
            "ALIASES_USED",
            "non-canonical labels resolved through the alias table: "
            + ", ".join(sorted(set(aliases_used))),
        )
    if not unknown and not duplicates and not missing:
        result(rows, "tree_sample_set", "INFO", "PASS", "tree resolves one terminal for every configured sample")
    return finish(rows, args.report)


def alignment_files(directory: str | Path, pattern: str) -> list[Path]:
    files = sorted(Path(directory).glob(pattern))
    if not files:
        raise DataError(f"{directory}: no files matched {pattern!r}")
    return files


def command_alignment(args: argparse.Namespace) -> int:
    rows: list[dict[str, object]] = []
    try:
        samples = load_samples(args.manifest)
        files = alignment_files(args.alignment_dir, args.glob)
        all_ids: set[str] = set()
        for path in files:
            all_ids.update(read_fasta(path))
        sequence_to_source = load_sequence_id_map(args.sequence_ids, all_ids, args.species_ids)
        source_to_sample = source_label_to_sample(samples)
    except DataError as error:
        result(rows, "alignment_setup", "ERROR", "FAIL", str(error))
        return finish(rows, args.report)

    result(
        rows,
        "alignment_file_count",
        "ERROR" if args.expected_groups is not None and len(files) != args.expected_groups else "INFO",
        "FAIL" if args.expected_groups is not None and len(files) != args.expected_groups else "PASS",
        f"found {len(files)} files matching {args.glob}",
        directory=str(args.alignment_dir),
    )
    expected_samples = {row["sample_id"] for row in samples}
    failures = 0
    for path in files:
        try:
            records = read_fasta(path)
            length = validate_equal_lengths(records, path)
            record_to_sample = resolve_record_samples(records, sequence_to_source, source_to_sample)
            mapped_samples = list(record_to_sample.values())
            duplicate_samples = sorted(
                sample_id for sample_id, count in Counter(mapped_samples).items() if count > 1
            )
            missing_samples = sorted(expected_samples.difference(mapped_samples))
            checks: list[str] = []
            expected_records = args.expected_records if args.expected_records is not None else len(samples)
            if len(records) != expected_records:
                checks.append(f"records={len(records)}, expected={expected_records}")
            if duplicate_samples:
                checks.append("duplicate_sample=" + ",".join(duplicate_samples))
            if missing_samples:
                checks.append("missing_sample=" + ",".join(missing_samples))
            if args.mode == "codon":
                if length % 3:
                    checks.append(f"alignment_length={length} is not divisible by 3")
                partial_gap_ids = []
                for record_id, (_, sequence) in records.items():
                    for start in range(0, len(sequence), 3):
                        codon = sequence[start : start + 3]
                        if "-" in codon and codon != "---":
                            partial_gap_ids.append(record_id)
                            break
                if partial_gap_ids:
                    checks.append("partial_gap_codon=" + ",".join(partial_gap_ids[:5]))
            if checks:
                failures += 1
                result(rows, "alignment_group", "ERROR", "FAIL", "; ".join(checks), group=path.name)
            else:
                result(
                    rows,
                    "alignment_group",
                    "INFO",
                    "PASS",
                    f"{len(records)} records; {length} aligned columns; one sequence per sample",
                    group=path.name,
                )
        except DataError as error:
            failures += 1
            result(rows, "alignment_group", "ERROR", "FAIL", str(error), group=path.name)
    result(
        rows,
        "alignment_summary",
        "ERROR" if failures else "INFO",
        "FAIL" if failures else "PASS",
        f"{len(files) - failures}/{len(files)} alignment groups passed",
    )
    return finish(rows, args.report)


def strip_trailing_empty(values: list[str]) -> list[str]:
    while values and not values[-1].strip():
        values.pop()
    return values


def read_family_matrix(path: str | Path) -> tuple[list[str], list[tuple[str, list[str]]], int]:
    """Read CAFE-style family matrix and return tips, family rows, start index."""

    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = strip_trailing_empty(next(reader))
        except StopIteration as error:
            raise DataError(f"{path}: empty family matrix") from error
        if not header:
            raise DataError(f"{path}: empty header")
        first = header[0].lstrip("#")
        start = 2 if first.lower() in {"desc", "description"} else 1
        if len(header) <= start:
            raise DataError(f"{path}: no species columns found")
        tips = header[start:]
        rows: list[tuple[str, list[str]]] = []
        for line_number, raw in enumerate(reader, start=2):
            values = strip_trailing_empty(raw)
            if not values:
                continue
            if len(values) != len(header):
                raise DataError(
                    f"{path}:{line_number}: {len(values)} columns; header has {len(header)}"
                )
            family = values[1] if start == 2 else values[0]
            if not family:
                raise DataError(f"{path}:{line_number}: empty FamilyID")
            rows.append((family, values[start:]))
    return tips, rows, start


def command_cafe(args: argparse.Namespace) -> int:
    rows: list[dict[str, object]] = []
    try:
        samples = load_samples(args.manifest)
        aliases = load_aliases(args.aliases)
        mapping, mapping_errors = alias_mapping(samples, aliases, args.context)
        tips, families, _ = read_family_matrix(args.input_matrix)
    except DataError as error:
        result(rows, "cafe_setup", "ERROR", "FAIL", str(error))
        return finish(rows, args.report)
    for error in mapping_errors:
        result(rows, "alias_mapping", "ERROR", "FAIL", error)

    expected = {row["sample_id"] for row in samples}
    resolved = [mapping.get(label) for label in tips]
    unknown = [label for label, sample_id in zip(tips, resolved) if sample_id is None]
    duplicate = sorted(
        sample_id for sample_id, count in Counter(item for item in resolved if item is not None).items() if count > 1
    )
    missing = sorted(expected.difference(item for item in resolved if item is not None))
    if unknown:
        result(rows, "cafe_tip_labels", "ERROR", "FAIL", "unmapped columns: " + ", ".join(unknown))
    if duplicate:
        result(rows, "cafe_tip_labels", "ERROR", "FAIL", "duplicate samples: " + ", ".join(duplicate))
    if missing:
        result(rows, "cafe_tip_labels", "ERROR", "FAIL", "missing samples: " + ", ".join(missing))
    if not unknown and not duplicate and not missing:
        result(
            rows,
            "cafe_tip_labels",
            "INFO",
            "PASS",
            f"{len(tips)} terminal columns resolve to the configured {len(samples)} samples",
        )

    seen_families: set[str] = set()
    invalid_counts = 0
    differential_failures = 0
    for family, values in families:
        if family in seen_families:
            invalid_counts += 1
            result(rows, "cafe_family_id", "ERROR", "FAIL", f"duplicate FamilyID {family}")
        seen_families.add(family)
        try:
            counts = [int(value) for value in values]
        except ValueError:
            invalid_counts += 1
            result(rows, "cafe_counts", "ERROR", "FAIL", f"{family}: non-integer count")
            continue
        if any(value < 0 for value in counts):
            invalid_counts += 1
            result(rows, "cafe_counts", "ERROR", "FAIL", f"{family}: negative count")
        if args.threshold is not None and max(counts) - min(counts) >= args.threshold:
            differential_failures += 1
    result(
        rows,
        "cafe_family_matrix",
        "ERROR" if invalid_counts else "INFO",
        "FAIL" if invalid_counts else "PASS",
        f"{len(families)} family rows; {invalid_counts} invalid rows",
    )
    if args.threshold is not None:
        result(
            rows,
            "cafe_threshold_filter",
            "ERROR" if differential_failures else "INFO",
            "FAIL" if differential_failures else "PASS",
            f"{differential_failures} families have max(count)-min(count) >= {args.threshold}",
        )

    if args.results_matrix:
        try:
            result_tips, result_families, _ = read_family_matrix(args.results_matrix)
            result_terminal_labels = [item.split("<", 1)[0] for item in result_tips if "<" in item and not item.startswith("<")]
            result_terminal_samples = [mapping.get(label) for label in result_terminal_labels]
            result_unknown = [
                label for label, sample_id in zip(result_terminal_labels, result_terminal_samples) if sample_id is None
            ]
            result_set = {family for family, _ in result_families}
            missing_from_input = result_set.difference(seen_families)
            result(
                rows,
                "cafe_results_family_set",
                "ERROR" if missing_from_input else "INFO",
                "FAIL" if missing_from_input else "PASS",
                f"{len(result_families)} result families; {len(missing_from_input)} absent from declared input",
            )
            if result_unknown:
                result(
                    rows,
                    "cafe_results_tip_labels",
                    "ERROR",
                    "FAIL",
                    "unmapped result tips: " + ", ".join(sorted(set(result_unknown))),
                )
            if len(result_families) != len(families):
                result(
                    rows,
                    "cafe_input_output_closure",
                    "WARN",
                    "ACTION_REQUIRED",
                    f"input has {len(families)} families but results contain {len(result_families)}; preserve an exclusion table and command log",
                )
        except DataError as error:
            result(rows, "cafe_results_matrix", "ERROR", "FAIL", str(error))
    return finish(rows, args.report)


def parser() -> argparse.ArgumentParser:
    main_parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = main_parser.add_subparsers(dest="command", required=True)

    metadata = subparsers.add_parser("metadata", help="validate taxon selection, grouping, rooting, and aliases")
    metadata.add_argument(
        "--manifest", "--samples", dest="manifest", required=True,
        help="config/phylogeny/taxa.tsv or an explicit terminal manifest",
    )
    metadata.add_argument("--aliases", default=None, help="optional context-specific label aliases")
    metadata.add_argument(
        "--expected-terminals", type=int, default=None,
        help="optional explicit assertion; otherwise the count is derived from the manifest",
    )
    metadata.add_argument(
        "--expected-species", type=int, default=None,
        help="optional explicit assertion; otherwise biological-species groups are counted from the manifest",
    )
    metadata.add_argument(
        "--require-root", action="store_true",
        help="fail if exactly one selected terminal is not marked as the root outgroup",
    )
    metadata.add_argument("--report", default=None, help="TSV path; default is stdout")
    metadata.set_defaults(func=command_metadata)

    tree = subparsers.add_parser("tree", help="validate Newick terminal labels against context aliases")
    tree.add_argument("--tree", required=True)
    tree.add_argument(
        "--manifest", "--samples", dest="manifest", required=True,
        help="taxon registry or explicit terminal manifest",
    )
    tree.add_argument("--aliases", default=None, help="optional context-specific label aliases")
    tree.add_argument("--context", default="all", help="analysis_context in the optional alias table")
    tree.add_argument(
        "--expected-terminals", type=int, default=None,
        help="optional explicit assertion; default is the selected manifest count",
    )
    tree.add_argument("--report", default=None, help="TSV path; default is stdout")
    tree.set_defaults(func=command_tree)

    alignment = subparsers.add_parser("alignment", help="validate SCO protein or codon alignments by SequenceIDs mapping")
    alignment.add_argument("--alignment-dir", required=True)
    alignment.add_argument("--glob", default="*.fa")
    alignment.add_argument("--sequence-ids", required=True, help="OrthoFinder WorkingDirectory/SequenceIDs.txt")
    alignment.add_argument(
        "--species-ids", default=None,
        help="OrthoFinder WorkingDirectory/SpeciesIDs.txt; defaults to the SequenceIDs.txt sibling",
    )
    alignment.add_argument(
        "--manifest", "--samples", dest="manifest", required=True,
        help="taxon registry or explicit terminal manifest",
    )
    alignment.add_argument("--mode", choices=("protein", "codon"), default="protein")
    alignment.add_argument("--expected-groups", type=int, default=None)
    alignment.add_argument(
        "--expected-records", type=int, default=None,
        help="optional explicit assertion; default is the selected manifest count",
    )
    alignment.add_argument("--report", default=None, help="TSV path; default is stdout")
    alignment.set_defaults(func=command_alignment)

    cafe = subparsers.add_parser("cafe", help="validate CAFE input/results matrices and alias closure")
    cafe.add_argument("--input-matrix", required=True)
    cafe.add_argument(
        "--manifest", "--samples", dest="manifest", required=True,
        help="taxon registry or explicit terminal manifest",
    )
    cafe.add_argument("--aliases", default=None, help="optional context-specific label aliases")
    cafe.add_argument("--context", default="all", help="analysis_context in the optional alias table")
    cafe.add_argument("--threshold", type=int, default=None, help="expect max-min < this threshold for every family")
    cafe.add_argument("--results-matrix", default=None, help="optional CAFE Base_count.tab")
    cafe.add_argument("--report", default=None, help="TSV path; default is stdout")
    cafe.set_defaults(func=command_cafe)

    return main_parser


def main() -> int:
    args = parser().parse_args()
    try:
        return args.func(args)
    except BrokenPipeError:
        return 0


if __name__ == "__main__":
    sys.exit(main())
