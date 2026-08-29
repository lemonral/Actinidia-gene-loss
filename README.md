# Synteny-informed gene-loss analysis across Actinidia

Reproducible workflows for comparative gene-loss analysis, genome collinearity,
phylogeny, and gene-family evolution across *Actinidia*. The study uses the
sister-lineage species *Clematoclethra scandens* as a syntenic reference and
analyses 23 chromosome-level assemblies, haplotypes, or subgenomes representing
12 *Actinidia* species. These 23 genomes are retained separately in per-genome
comparisons and grouped into 13 biological lineages for branch-based analyses.

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
- SynOrths- and tBLASTX-based `retained`, `decayed`, and `deleted`
  classification for 23 genomes;
- independent Miniprot evidence used to identify a strict pseudogenized subset
  without changing the primary `decayed` or `deleted` class;
- genome- and lineage-level summaries of shared, species-specific, tree-node,
  and recurrent loss events;
- chromosome-level and within-chromosome analysis of decayed loci;
- associations of gene loss with four-tissue mean expression and *C. scandens*
  gene copy number, defined as CD-HIT 90% cluster size;
- covariate-adjusted GO and KEGG analysis of species-specific loss events across
  13 biological lineages;
- NLR repertoire, loss, structural-class, and branch-event analyses;
- OrthoFinder, IQ-TREE, ASTRAL-Pro, MCMCTree, and CAFE5 workflows for species
  phylogeny and gene-family evolution; and
- reproducible generation of analysis figures and tables.

## Gene-loss definitions

The primary classification applies the same thresholds to all 23 genomes:

| State | Evidence |
| --- | --- |
| `retained` | Exact SynOrths-supported orthologous anchor |
| `decayed` | A missing-gene candidate with at least one genome-wide tBLASTX hit at identity >= 50%, bit score >= 50, and E-value < 1e-5 |
| `deleted` | A missing-gene candidate with no qualifying tBLASTX hit |
| `not_called_loss` | Outside the resolved candidate scope |

The primary positive-loss state is `decayed + deleted`; its resolved
denominator is `retained + decayed + deleted`. Miniprot evidence is joined
afterward to distinguish frameshift-only, in-frame-stop-only, combined, and
mechanism-unresolved decayed calls. Decayed genes passing the coding-disruption
and alignment-quality criteria form the strict pseudogenized subset. This
evidence refines the possible mechanism but does not change the primary loss
class.

The 23 genomes are preserved separately for per-genome comparisons. For
branch-based functional and NLR analyses, haplotypes and subgenomes from the
same species are combined, while the two parental haplomes of the F1 hybrid
*A. zhejiangensis* are retained separately, giving 13 biological lineages. A
multi-genome lineage is considered completely lost only when every constituent
haplotype or subgenome is `decayed` or `deleted`; mixed states are reported as
partial or homeolog-specific loss.

For functional analysis, the foreground is the set of complete losses assigned
specifically to the focal lineage. Genes assigned to tree-node losses along the
focal root-to-tip path are removed first. The remaining risk-set genes form the
comparison background after annotation and covariate filtering. Logistic score
tests adjust for four-tissue mean expression and *C. scandens* gene copy number
(CD-HIT 90% cluster size).

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
scripts/gene_loss/   loss matrices, per-genome comparisons, and lineage summaries
scripts/spatial/     chromosome and within-chromosome position analyses
scripts/downstream/  expression, copy-number, and loss-event preparation
scripts/function/    GO and KEGG enrichment analyses
scripts/statistics/  exact ploidy-group and related statistical comparisons
scripts/nlr/         NLR-Annotator preparation, validation, and loss summaries
scripts/phylogeny/   phylogeny, divergence-time, and CAFE5 preparation/validation
scripts/figures/     publication figure renderers
```

## Installation

Python 3.10 or later is required.

```bash
git clone https://github.com/lemonral/Actinidia-gene-loss.git
cd Actinidia-gene-loss

python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[all]'
```

Confirm the installation with:

```bash
geneloss --help
python -m pytest -q
```

The complete workflow also uses stage-specific external programs, including
BLAST+, Miniprot, minimap2, JCVI, BUSCO, SynOrths, OrthoFinder, MAFFT,
IQ-TREE 2, ASTRAL-Pro, PAML/MCMCTree, CAFE5, and NLR-Annotator. Exact
tool roles and validation requirements are described in
[`docs/PHYLOGENY_TOOLCHAIN.md`](docs/PHYLOGENY_TOOLCHAIN.md) and
[`docs/PIPELINE.md`](docs/PIPELINE.md).

Publication figures require Arial regular, italic, bold, and bold-italic faces.
Arial is not redistributed here. If it is installed outside a standard system
font directory, set `ARIAL_FONT_DIR` before rendering:

```bash
export ARIAL_FONT_DIR=/path/to/licensed/arial-font-files
```

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

Public assets, expected checksums, genome identities, biological
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
3. [`docs/LOSS_CLASSIFICATION.md`](docs/LOSS_CLASSIFICATION.md)
   for the primary loss-classification rules;
4. [`docs/DECAYED_CHROMOSOME_POSITION_ANALYSIS.md`](docs/DECAYED_CHROMOSOME_POSITION_ANALYSIS.md)
   for spatial analysis;
5. [`docs/COVARIATE_ADJUSTED_FUNCTIONAL_ENRICHMENT.md`](docs/COVARIATE_ADJUSTED_FUNCTIONAL_ENRICHMENT.md)
   for the species-specific-loss GO and KEGG models;
6. [`docs/TREE_AWARE_LOSS_FUNCTION.md`](docs/TREE_AWARE_LOSS_FUNCTION.md) for
   branch-event and complementary functional analyses.

A concise record of analysis and figure changes is maintained in
[`CHANGELOG.md`](CHANGELOG.md).

Every stage-specific script supports `--help`, for example:

```bash
python scripts/qc/extract_primary_annotation.py --help
python scripts/gene_loss/prepare_translated_search_candidates.py --help
python scripts/gene_loss/compare_within_species_genome_units.py --help
python scripts/spatial/analyze_decayed_chromosome_distribution.py --help
python scripts/function/run_terminal_covariate_adjusted_enrichment.py --help
python scripts/statistics/compare_ploidy_loss_rates.py --help
python scripts/nlr/run_nlr_annotator_batch.py --help
```

The manuscript functional analysis uses the biological-lineage profile:

```text
--analysis-level biological_species
--terminal-node-type species_terminal
--expected-terminals 13
--expected-terminal-event-memberships 19192
```

These arguments are explicit safeguards: a run made with the default
23-genome exploratory profile is not the analysis reported in the manuscript.

## Reproducibility

The workflow is designed to fail closed when sample identifiers, schemas,
checksums, sequence sets, or expected row counts do not match. Configuration
files preserve biological groupings and analysis parameters, while output
bundles record the corresponding input and tool checksums. Tests use synthetic
data and do not contain study genomes or manuscript files.

Run the test suite with:

```bash
python -m pytest -q
```

## Results and data availability

Large analysis outputs are archived separately because of file size and
third-party redistribution restrictions. The expected structure of a curated
result bundle is documented in
[`results/README.md`](results/README.md). A typical figure bundle contains the
PNG and PDF figure, plotted data, caption, validation summary, and checksum
manifest.

## Citation

If you use this workflow, please cite the associated *Actinidia* comparative
genomics study. Full citation details will be added after publication.

## License

The source code and workflow documentation in this repository are available
under the [MIT License](LICENSE). Third-party genome assemblies, annotations,
software, and other externally sourced materials remain subject to their
original licenses and terms of use.
