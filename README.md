# Reproducible Actinidia gene-loss revision

This repository is a clean, manifest-driven rebuild of the comparative
gene-loss analyses used in the associated Actinidia manuscript. It contains
code, configuration, documentation, tests, and curated small outputs only.
Genome assemblies, annotations, raw alignments, and other large files live in
an external data store and are never committed to Git.

## Current revision goals

1. Evaluate the 2026 haplotype-resolved *Actinidia eriantha* HAP1 and HAP2
   assemblies as a paired biological individual.
2. Replace the previous unversioned *A. deliciosa* A-F inputs with the 2025
   haplotype-resolved hexaploid assembly, if all six matched haplotype bundles
   can be validated.
3. Evaluate the three sequence-distinct downloaded *A. rufa* candidates (ARU,
   Fuchu, and ActinidiaBase v1) and exclude the species from downstream
   gene-loss analyses if chromosome-matched/callable coverage is inadequate.
   The unnamed ActinidiaBase plant is an assembly sensitivity, not a biological
   replicate.
4. Rebuild all shared/non-shared loss summaries from the accepted cohort. No
   historical shared-loss count or denominator is hard-coded.
5. Build a separate biological-species phylogeny and gene-family workflow with
   OrthoFinder 3, IQ-TREE 2, ASTRAL-Pro, MCMCTree, and CAFE 5. The compact
   selected-lineage cohort closes no fossil bracket, so the author-approved
   revised dated tree uses four checksum-bound TimeTree 5 secondary soft
   intervals and is described as secondary-calibrated, never fossil-calibrated.
   The approved minimal-outgroup contract reuses one existing *Rhododendron*,
   *Coffea*, and *Vitis*; treats *Leea* as a conditional root diagnostic and
   *Catharanthus* as a one-for-one *Coffea* sensitivity; prohibits *Saurauia*;
   and downloads no further outgroup taxa.
   Assembly-unit trees remain diagnostic, phylogeny-only outgroups never enter
   the gene-loss denominator, and species-level PGLS runs only after the dated
   biological-species tree and callable shared/non-shared aggregation pass.
   Ordinary Gaussian PGLS is retained as an exploratory phylogenetic trait
   sensitivity and is publication-blocked until a denominator-aware
   phylogenetic count model is predeclared and validated. The validated CAFE 5
   result is the Base Poisson model; its Gamma3 initialization failure is
   retained transparently and is not reported as Gamma support.

## Repository and data layout

On the analysis server the repository and large-data store are siblings:

```text
actinidia_gene_loss_rebuild/   # this Git repository
actinidia_gene_loss_data/      # downloads, standardized inputs, work, large results
```

The repository contains a relative `data_store` symlink for server use. A
clone on another machine may instead set `DATA_ROOT` in `config/project.env`.
See `docs/DATA_LAYOUT.md`.

Legacy *Rhododendron*, *Coffea*, and *Vitis* files are retained as soft-linked
phylogeny assets, with checksums and provenance state recorded separately.
See `docs/OUTGROUPS_AND_CALIBRATIONS.md` and
`config/phylogeny/minimal_outgroup_design.tsv`. Their presence does not add
them to the gene-loss cohort or activate a fossil calibration.

## Analysis units are not biological replicates

Haplotypes and subgenomes are processed separately for assembly QC, synteny,
orthology, and loss calling. Species-level statistics then aggregate those
assembly units according to explicit rules in `config/cohorts.tsv`.

For diploid *A. eriantha*, HAP1 and HAP2 come from the same individual. A
species-supported loss therefore requires evidence from both haplotypes;
one-sided absence is reported as haplotype-specific or uncertain rather than
as an independent lineage loss.

## Two explicit loss-evidence layers

The article-comparable main trend uses one rule for all old and new units. An
exact SynOrths target anchor is `retained`; a missing-gene candidate with a
genome-wide tBLASTX hit at identity >=50%, bit score >=50, and e-value `<1e-5`
is `decayed`; and a candidate without such a hit is `deleted`. There is no
alignment-length minimum. Main trend numerators use `decayed + deleted`, while
resolved denominators additionally include retained. Rows outside the
historical candidate scope are not called. Completed Miniprot evidence only
subdivides possible causes and never rewrites these article classes. See
`docs/MANUSCRIPT_METHOD_LOSS_CLASSIFICATION.md`.

The stricter callable-aware grid remains a sensitivity/mechanism layer. It
requires a callable local interval for deletion and an explicit high-quality
frameshift or in-frame stop for `pseudogenized`; other local sequence and
non-callable cases are uncertain. The chromosome-position analysis uses this
strict layer because it requires an observed target locus: strict
`pseudogenized` loci are compared with observed retained loci, while deleted
expected-locus midpoints remain a labelled sensitivity. See
`docs/UNIFORM_LOSS_EVIDENCE.md`.

GO and KEGG main results use the article-method classes and the frozen matching
topology. Complete lineage losses are separated into terminal-branch,
internal-branch, and repeated independent events; mixed homeolog/subgenome
states remain partial. Frozen eggNOG-mapper annotations, GO definitions, and
KEGG KO descriptions are checksum-bound. The former strict enrichment remains
a separately labelled sensitivity.

Publication figures use concise display labels without altering internal unit
IDs or aggregation. The display-only qualifier `unphased` and the
multiplication sign in *A. zhejiangensis* are omitted consistently. Upright
A/B suffixes remain where needed to distinguish its two curated parental
lineages. Latin binomials remain italic and other informative assembly-unit
suffixes remain upright.

## Quick start

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[all]'
cp config/project.env.example config/project.env
pytest -q
```

Download and analysis commands are enabled only after URLs and scopes in
`config/assemblies.tsv` reconcile with `config/downloads.tsv` and the generated
download report. Publisher checksums, declared byte sizes, and calculated local
SHA-256 values are retained in the download manifest/report rather than copied
into multiple configuration files.

The complete stage-by-stage workflow, decision rules, outputs, and uncertainty
handling are described in `docs/PIPELINE.md`. Frozen numerical and categorical
choices are machine-readable in `config/analysis_parameters.toml` and explained
in `docs/DECISION_RULES.md`. A directory-by-directory orientation and complete
data-flow overview are provided in `docs/REPOSITORY_GUIDE.md`.
Publication-label and final-figure composition rules are documented in
`docs/PUBLICATION_FIGURE_POLICY.md`.
The path-free mapping from revised claims to canonical result bundles and
publication gates is maintained in `docs/REVISION_EVIDENCE_INDEX.md`.
The production primary-isoform/CDS/protein extraction contract, including the
independent gffread sequence gate, is documented in
`docs/PRIMARY_ANNOTATION_STANDARDIZATION.md`.
The subsequent exact-ID and publisher-protein sequence gate is documented in
`docs/PUBLISHED_PROTEIN_COMPATIBILITY.md`.
Publisher bundles whose protein accessions differ from their GFF3 transcript
IDs first use the audited, sequence-preserving primary-subset adapter described
in `docs/PUBLISHER_PRIMARY_PROTEIN_REMAP.md`; IDs are never inferred by
truncation or similarity.
The broader public assembly inventory and its deduplication rules are documented
in `docs/CANDIDATE_ASSEMBLY_CATALOG.md`; inventory membership is not analysis
inclusion.
The chromosome-naming policy and corrected Hongyang v4.0 HY4A/HY4P reference
identities are documented in `docs/CHROMOSOME_HOMOLOGY_RENUMBERING.md`.
Production names are the complete global one-to-one maximum-nucleotide-
similarity assignment to HY4A: every unit must contain each `Chr01`--`Chr29`
label exactly once. Absolute coverage, reciprocal-best/separation measures,
HY4P agreement, and independent JCVI anchors are retained as confidence
diagnostics and do not block that naming rule. Publisher sequence direction is
preserved. The relabelling worker writes new matched FASTA/GFF files, never
edits publisher inputs, and requires exact CDS/protein sequence closure. The
older strict four-matrix `PASS_AUTO` gate remains available as a labelled
diagnostic only; it is not the production naming authority.

## Publication boundary

This repository and its local revision commit remain private until the
manuscript authors approve a public release. Reviewer-response working text may
be retained only in that private history; it must not be pushed or included in
a public release. Manuscript source files, supplementary workbooks, server
paths, credentials, raw genomes, and large tool outputs must never be committed.
Nothing is pushed automatically.
