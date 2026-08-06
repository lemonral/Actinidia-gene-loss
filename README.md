# Actinidia gene-loss analysis

Reproducible workflows for comparative gene-loss analysis, genome collinearity,
phylogeny, and gene-family evolution in *Actinidia*. The analyses use
*Clematoclethra scandens* as the reference gene set and retain chromosome-level
assemblies, haplotypes, and subgenomes as explicit genome units.

This repository contains the Python code, machine-readable configuration,
workflow documentation, and synthetic tests used in the study. Genome
assemblies, annotations, sequencing reads, and large intermediate files are not
redistributed here.

## Analysis overview

The workflow includes:

- checksum- and manifest-based input validation;
- primary-isoform extraction and assembly/annotation quality control;
- chromosome correspondence assessment with nucleotide similarity and JCVI
  collinearity;
- SynOrths- and tBLASTX-based gene-loss classification for 23 genome units;
- independent Miniprot evidence for frameshifts and in-frame stop codons;
- genome-unit and species-lineage summaries of shared, terminal, internal, and
  recurrent loss events;
- chromosome-level and within-chromosome analysis of decayed loci;
- associations of gene loss with four-tissue mean expression and reference
  protein-family size;
- covariate-adjusted GO and KEGG analysis of terminal loss events;
- NLR repertoire, loss, structural-class, and branch-event analyses;
- OrthoFinder, IQ-TREE, ASTRAL-Pro, MCMCTree, and CAFE5 workflows for species
  phylogeny and gene-family evolution; and
- reproducible generation of publication figures and tables.

## Gene-loss definitions

The primary classification applies the same thresholds to all genome units:

| State | Evidence |
| --- | --- |
| `retained` | Exact SynOrths-supported orthologous anchor |
| `decayed` | A missing-gene candidate with at least one genome-wide tBLASTX hit at identity >= 50%, bit score >= 50, and E-value < 1e-5 |
| `deleted` | A missing-gene candidate with no qualifying tBLASTX hit |
| `not_called_loss` | Outside the resolved candidate scope |

The primary positive-loss state is `decayed + deleted`; its resolved
denominator is `retained + decayed + deleted`. Miniprot evidence is joined
afterward to distinguish frameshift-only, in-frame-stop-only, combined, and
mechanism-unresolved decayed calls. This evidence refines the possible
mechanism but does not change the primary loss class.

Genome units are preserved separately for unit-level comparisons. For
branch-based functional and NLR analyses, units are grouped into biological
lineages. A multi-unit lineage is considered completely lost only when every
constituent haplotype or subgenome is `decayed` or `deleted`; mixed states are
reported as partial or homeolog-specific loss.

## Repository structure

```text
config/     machine-readable cohorts, parameters, manifests, and mappings
docs/       methods, decision rules, data contracts, and workflow descriptions
scripts/    stage-specific command-line programs
src/        reusable Python package (geneloss_repro)
tests/      synthetic unit and regression tests
results/    description of the curated result-release structure
```

The main script groups are:

```text
scripts/qc/          assembly, annotation, BUSCO, SynOrths, JCVI, and chromosome QC
scripts/gene_loss/   loss matrices, genome-unit comparisons, and lineage summaries
scripts/spatial/     chromosome and within-chromosome position analyses
scripts/downstream/  expression, copy-family, and loss-event preparation
scripts/function/    GO and KEGG enrichment analyses
scripts/nlr/         NLR-Annotator preparation, validation, and loss summaries
scripts/phylogeny/   phylogeny, divergence-time, and CAFE5 preparation/validation
scripts/figures/     publication figure renderers
```

## Installation

Python 3.10 or later is required.

```bash
git clone git@github.com:lemonral/Actinidia-gene-loss.git
cd Actinidia-gene-loss

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[all]'
```

Confirm the installation with:

```bash
geneloss --help
pytest -q
```

The complete workflow also uses stage-specific external programs, including
BLAST+, Miniprot, minimap2, JCVI, BUSCO, SynOrths, OrthoFinder, MAFFT,
IQ-TREE 2, ASTRAL-Pro, PAML/MCMCTree, CAFE5, and NLR-Annotator. Exact
tool roles and validation requirements are described in
[`docs/PHYLOGENY_TOOLCHAIN.md`](docs/PHYLOGENY_TOOLCHAIN.md) and
[`docs/PIPELINE.md`](docs/PIPELINE.md).

## Data configuration

Large data are kept outside the Git repository. Create a local configuration
file and point `DATA_ROOT` to an external data directory:

```bash
cp config/project.env.example config/project.env
```

Example layout:

```text
Actinidia-gene-loss/
actinidia_gene_loss_data/
├── downloads/
├── raw/
├── standardized/
├── work/
├── results_large/
└── logs/
```

Public assets, expected checksums, genome-unit identities, biological
groupings, and analysis parameters are declared in `config/`. Local paths and
credentials must not be added to these tracked files. See
[`docs/DATA_LAYOUT.md`](docs/DATA_LAYOUT.md) for the complete storage and
provenance policy.

## Running the workflow

The reusable command-line package provides core operations:

```bash
geneloss extract-annotation --help
geneloss normalize-synorth --help
geneloss call-candidates --help
geneloss classify-tblastx --help
geneloss build-loss-master --help
geneloss spatial-summary --help
```

Study-scale analyses are manifest-driven and are run in validated stages rather
than through a single monolithic command. Each stage checks its input schema,
sample identities, record counts, and checksums before publishing an output.
Start with:

1. [`docs/REPOSITORY_GUIDE.md`](docs/REPOSITORY_GUIDE.md) for project
   orientation;
2. [`docs/PIPELINE.md`](docs/PIPELINE.md) for the end-to-end workflow;
3. [`docs/MANUSCRIPT_METHOD_LOSS_CLASSIFICATION.md`](docs/MANUSCRIPT_METHOD_LOSS_CLASSIFICATION.md)
   for the primary loss-classification rules;
4. [`docs/DECAYED_CHROMOSOME_POSITION_ANALYSIS.md`](docs/DECAYED_CHROMOSOME_POSITION_ANALYSIS.md)
   for spatial analysis;
5. [`docs/COVARIATE_ADJUSTED_FUNCTIONAL_ENRICHMENT.md`](docs/COVARIATE_ADJUSTED_FUNCTIONAL_ENRICHMENT.md)
   for the terminal-loss GO and KEGG models;
6. [`docs/TREE_AWARE_LOSS_FUNCTION.md`](docs/TREE_AWARE_LOSS_FUNCTION.md) for
   branch-event and complementary functional analyses.

A concise record of analysis and figure changes is maintained in
[`docs/ANALYSIS_CHANGELOG.md`](docs/ANALYSIS_CHANGELOG.md).

Every stage-specific script supports `--help`, for example:

```bash
python scripts/qc/extract_primary_annotation.py --help
python scripts/gene_loss/prepare_translated_search_candidates.py --help
python scripts/gene_loss/compare_within_species_genome_units.py --help
python scripts/spatial/analyze_decayed_chromosome_distribution.py --help
python scripts/function/run_terminal_covariate_adjusted_enrichment.py --help
python scripts/nlr/run_nlr_annotator_batch.py --help
```

## Reproducibility

The workflow is designed to fail closed when sample identifiers, schemas,
checksums, sequence sets, or expected row counts do not match. Configuration
files preserve biological groupings and analysis parameters, while output
bundles record the corresponding input and tool checksums. Tests use synthetic
data and do not contain study genomes or manuscript files.

Run the test suite with:

```bash
pytest -q
```

## Results and data availability

Large analysis outputs are stored separately because of file size and
third-party redistribution restrictions. The expected structure and contents
of curated result bundles are documented in
[`results/README.md`](results/README.md). A typical figure bundle contains the
PNG and PDF figure, plotted data, caption, validation summary, and checksum
manifest.

## Citation

If you use this workflow, please cite the associated *Actinidia* comparative
genomics study. Full citation details will be added after publication.
