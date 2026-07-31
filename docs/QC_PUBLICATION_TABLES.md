# Assembly and annotation QC publication tables

`scripts/qc/build_qc_publication_tables.py` is the publication gate between
the large, private QC workspace and the small reviewer-facing tables. It does
not run BUSCO or recalculate FASTA statistics. It reconciles and validates the
completed outputs before publishing a path-free bundle.

## Inputs

The command requires four TSVs:

1. An explicit metadata/decision table keyed by `assembly_unit_id`.
2. The exact output schema of `scripts/qc/basic_stats.py`, keyed by `sample`.
3. The exact output schema of `scripts/qc/collect_busco.py` for genome BUSCO.
4. The exact output schema of `scripts/qc/collect_busco.py` for protein BUSCO.

The metadata table is also the one-to-one identity map from `qc_sample` in the
three producer tables to the reader-facing `assembly_unit_id`. Its exact
columns are:

```text
assembly_unit_id
qc_sample
biological_species
individual_id
haplotype_or_subgenome
accession
decision_status
publisher_assembly_scope
local_qc_scope
publisher_assembly_provenance
publisher_annotation_provenance
publisher_protein_provenance
decision_reason
```

`decision_status` is `current`, `candidate`, or `excluded`. Candidate and
excluded rows require an explicit reason. They are retained in every output;
the QC table is an audit of all assessed assemblies, not a filtered analysis
cohort.

`publisher_assembly_scope` describes what the publisher says was released.
`local_qc_scope` describes the exact genome/annotation/protein scope used for
the locally recomputed statistics. These may differ, for example when a
publication describes a primary assembly but only chromosome-anchored records
have matched public annotation files. The three publisher-provenance fields
should contain stable accessions, DOIs, or public URLs, never local paths.

The accession in the basic-statistics input must exactly match the accession in
the metadata row. This catches a common error in which the metrics from one
assembly are paired with the decision record for another.

## Validation gate

The command stops without creating an output directory unless all of the
following checks pass:

- both `assembly_unit_id` and `qc_sample` are nonempty and unique;
- basic statistics, genome BUSCO, and protein BUSCO contain exactly the
  declared `qc_sample` set;
- the BUSCO version, lineage dataset, dataset creation date, and lineage count
  are uniform across genome and protein runs;
- the genome and protein BUSCO modes are uniform within their own table and
  represent different run roles;
- BUSCO `C = S + D` and `S + D + F + M = n` for every row;
- every BUSCO percentage agrees with its count and denominator at the printed
  precision;
- genome length, sequence count, N50/L50, N and GC percentages, GFF feature
  counts, and protein-set metrics are internally consistent;
- retained output fields contain no recognizable private runtime path.

`C_count` may be blank in a raw BUSCO collector row because some BUSCO v5
short summaries omit the total-complete count line. The builder derives it as
`S_count + D_count`, validates the reported complete percentage against that
value, and writes the derived count explicitly in the normalized outputs.

## Command

```bash
python scripts/qc/build_qc_publication_tables.py \
  --metadata config/qc_assembly_decisions.tsv \
  --basic-stats "$DATA_ROOT/qc/basic_stats.tsv" \
  --genome-busco "$DATA_ROOT/qc/genome_busco.tsv" \
  --protein-busco "$DATA_ROOT/qc/protein_busco.tsv" \
  --output-dir results/qc/assembly_annotation_qc
```

The output directory must be absent or empty. All files are first written and
fsynced in a sibling staging directory, then the complete directory is
published with one rename. A nonempty destination is refused rather than
overwritten.

## Outputs

The bundle contains:

- `qc_metadata_public.tsv`: normalized identity, decision, scope, and
  publisher provenance;
- `qc_basic_stats_public.tsv`: all path-free locally recomputed assembly,
  annotation, and protein statistics;
- `qc_genome_busco_public.tsv`: normalized genome BUSCO results;
- `qc_protein_busco_public.tsv`: normalized protein BUSCO results;
- `assembly_annotation_qc_supplementary.tsv`: the combined reviewer-ready
  table, including genome size, sequence count, N50/L50, N content, annotated
  genes, protein count, and genome/protein BUSCO C/S/D/F/M percentages, counts,
  and `n`;
- `qc_publication_validation.json`: the path-free pass report, input basenames,
  SHA-256 checksums, and decision-status counts.

The four normalized TSVs are ordered by the metadata table. The normalized
metadata, basic-statistics, and BUSCO tables can be passed directly to
`scripts/figures/make_qc_figure.py`. Runtime columns named `genome_path`,
`gff_path`, `protein_path`, `input_path`, and `short_summary_path` are never
published. Shell commands are not accepted as input metadata and are not
written to any output.

## Current revised Actinidia bundle

The reviewer-facing revised bundle contains 23 independent analysis units:
four *A. arguta* haplotypes (A-D), two Hongyang v4 *A. chinensis* haplomes
(HY4A/HY4P), six newly selected *A. deliciosa* ADM units (A-F), two newly
selected *A. eriantha* haplotypes (HAP1/HAP2), two *A. × zhejiangensis*
parental-lineage units (A/B), one selected *A. rufa* assembly, and the retained
single-unit legacy species. HAP1/HAP2 are paired haplotypes from one individual
and are not biological replicates.

`scripts/qc/prepare_primary_actinidia_qc_inputs.py` imports legacy values only
after exact genome real-path identity and protein SHA-256 identity checks, then
combines them with QC results for new or changed *Actinidia* inputs. It never
runs BUSCO. The current table therefore uses a single BUSCO version/lineage
while avoiding unnecessary reruns of unchanged assemblies.

Assembly-wide QC and chromosome-dependent analysis scope are separate fields.
For the selected *A. rufa* release, all 38 released nuclear sequences are
included in assembly-wide statistics and BUSCO. The nine unplaced contigs are
retained in the audit, whereas JCVI and chromosome-position analyses use the
declared 29-pseudomolecule subset. Publisher gene counts and analyzed primary
coding-gene counts are also reported separately, so omitted invalid models or
unplaced records cannot be mistaken for unreported data.
