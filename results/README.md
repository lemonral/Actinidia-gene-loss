# Result-bundle format

Study outputs are archived separately from this code repository because the
underlying genomes, alignments, and intermediate files are large and may be
subject to third-party redistribution restrictions. This directory documents
the expected structure of the curated result release and contains only small
tables needed for reproducibility.

## Bundle contents

A curated figure bundle contains:

```text
figure_name.png
figure_name.pdf
figure_name.plot_data.tsv
figure_name.caption.txt
figure_name.validation.json
figure_name.manifest.json
```

The plot-data table contains the exact values supplied to the renderer. The
validation file records schema and biological closure checks. The manifest
contains relative filenames, byte counts, and SHA-256 checksums. Runtime paths
and credentials are never written to a curated bundle.

A curated table bundle normally contains:

```text
summary.tsv
validation.json
manifest.json
```

Additional row-level tables are included only when they are small enough for
redistribution and are required to reproduce a reported summary. The
`tables/ploidy_comparison` directory contains the 23-genome input and exact
statistics reported in Table S12; each haplotype, subgenome, or single-genome
assembly is one observation. Small tables tracked directly in this repository
may use a README plus a deterministic generating script instead of duplicating
the complete archived-bundle metadata.

## Analysis-set conventions

The primary loss matrix evaluates 35,547 *C. scandens* reference genes in 23
*Actinidia* genomes. Its positive state is `decayed + deleted`, and its resolved
denominator is `retained + decayed + deleted`. Strict pseudogenized evidence is
a nested mechanistic subset and is not added to that numerator.

For the primary functional analysis, the 23 genomes are grouped into 13
biological lineages. Complete loss in a multi-genome lineage requires every
constituent haplotype or subgenome to be `decayed` or `deleted`. For each focal
lineage, genes assigned to tree-node losses on its root-to-lineage path are
removed. Complete species-specific lost genes form the foreground; the other
annotation- and covariate-complete genes in that risk set form the background.
Logistic score tests adjust for four-tissue mean expression and *C. scandens*
gene copy number, defined as CD-HIT 90% cluster size.

The expression/copy-number analysis pools shared and non-shared `decayed`
calls. Its denominator is `retained + decayed + deleted` within each genome and
bin or copy-number class; `not_called_loss` is excluded. This CD-HIT cluster
size is not a target-species CNV estimate.

## Analysis groups

The separately archived result release is organized into these broad groups:

- assembly, annotation, BUSCO, chromosome-homology, and JCVI quality control;
- the 23-genome primary loss matrix and shared/non-shared summaries;
- chromosome-level and five-zone analyses of localized decayed loci;
- four-tissue expression and gene-copy-number associations;
- species-specific complete-loss functional analysis across 13 biological
  lineages;
- exact ploidy-group comparisons using all 23 genomes as observations;
- NLR repertoire, loss-class, and species-specific/tree-node event summaries;
- OrthoFinder, species-tree, divergence-time, and CAFE5 outputs; and
- publication-ready figures and supplementary tables.

Definitions of the analysis sets and denominators are given in
[`../docs/LOSS_CLASSIFICATION.md`](../docs/LOSS_CLASSIFICATION.md). The
provenance map in [`../docs/RESULTS_PROVENANCE.md`](../docs/RESULTS_PROVENANCE.md)
links each analysis group to its validation contract.

## Validation requirements

A result bundle is considered complete only when:

1. all required inputs match their declared sample identifiers and checksums;
2. row counts and denominators satisfy the analysis-specific closure rules;
3. output filenames are relative and contain no private runtime paths;
4. every declared output matches the checksum manifest; and
5. figures and tables have been visually inspected after rendering.

Raw BLAST, Miniprot, SynOrths, JCVI, BUSCO, OrthoFinder, NLR-Annotator, and
phylogenetic working directories remain outside the curated release.
