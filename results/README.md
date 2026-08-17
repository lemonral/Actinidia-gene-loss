# Result-bundle format

Study outputs are archived separately from this code repository because the
underlying genomes, alignments, and intermediate files are large and may be
subject to third-party redistribution restrictions. This directory documents
the expected structure of the curated result release. It does not contain the
study data itself.

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
validation file records the schema and biological closure checks. The manifest
contains relative filenames, byte counts, and SHA-256 checksums. Runtime paths
and credentials are never written to a curated bundle.

A curated table bundle normally contains:

```text
summary.tsv
validation.json
manifest.json
```

Additional row-level tables are included only when they are small enough for
redistribution and are required to reproduce a reported summary.

## Analysis groups

The separately archived result release is organized into these broad groups:

- assembly, annotation, BUSCO, chromosome-homology, and JCVI quality control;
- the 23-unit primary loss matrix and shared/non-shared summaries;
- chromosome-level and five-zone analyses of localized decayed loci;
- four-tissue expression and reference-family-size associations;
- terminal complete-loss functional analysis across 13 biological lineages;
- NLR repertoire, loss-class, and branch-event summaries;
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
