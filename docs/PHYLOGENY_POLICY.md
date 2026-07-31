# Phylogeny and gene-family policy

This policy separates biological inference from technical assembly diagnostics.
The machine-readable contracts are `config/phylogeny/analysis_designs.tsv`,
`config/phylogeny/representation_policy.tsv`, and
`config/phylogeny/taxa.tsv`.

## Manuscript method and revised production method

The submitted manuscript compared 17 tips with OrthoFinder 2.5.5, DIAMOND,
DendroBLAST gene trees, a concatenated MAFFT/back-translated alignment,
RAxML-NG, `ape::chronos` with TimeTree points, and CAFE 4.2. The selected tips
included manuscript-era genome complements such as *A. arguta* C and
*A. deliciosa* D, while the two parental haplomes of the F1 hybrid
*A. zhejiangensis* were both displayed. That exact method remains available
only as `legacy_manuscript_reproduction`; it is not described as the revised
biological-species analysis.

The revised production design uses OrthoFinder 3, per-locus alignments and
gene trees, IQ-TREE 2 for a partitioned concatenation analysis, ASTRAL-Pro for
a multi-copy coalescent analysis, a fossil-ready biological-species topology,
MCMCTree dating, and CAFE 5. Exact versions, binaries, commands, checksums, thresholds, random
seeds, and support-replicate counts must be frozen in each production run.
The OrthoFinder internal species tree is retained as a diagnostic and is not
silently substituted for either accepted production topology.

The outgroup set is now bounded by
`config/phylogeny/minimal_outgroup_design.tsv`: one QC-selected existing
*Rhododendron*, existing *Coffea arabica*, and existing *Vitis vinifera* form
the required candidates; one *Leea coccinea* haplotype is a conditional
root-stability addition; and *Catharanthus roseus* is a declared one-for-one
*Coffea* rooting swap rather than a default extra tip. No additional outgroup
taxon is acquired. *Saurauia* is prohibited from the nuclear matrix.

## Four products, four tip contracts

### Primary selected-lineage tree

This is the only topology eligible for primary dating, CAFE 5, and PGLS after
its dated-tree and lineage-loss gates pass. It contains one declared
representative haploid complement per nonhybrid biological species, both
author-designated parental lineages of *A. x zhejiangensis* as separately
labelled A and B tips, and one validated proteome per outgroup species. Other
haplotypes, subgenomes, assemblies, accessions, or individuals do not increase
the primary sample size.

### Assembly-unit diagnostic tree

This tree retains HAP1/HAP2, *A. arguta* A-D, *A. deliciosa* A-F, and other
accepted units with visible suffixes. It detects contamination, switches,
unexpected placement, and orthologue-occupancy problems. It is not dated,
does not enter CAFE, and never increases PGLS sample size.

### Representative-swap sensitivities

Each sensitivity replaces one species representative, rebuilds OrthoFinder
and every affected locus, and reruns both topology methods. Results are never
patched into the primary alignment. A sensitivity is dated, passed to CAFE, or
used for PGLS only with its matching topology, sequence set, family counts, and
species-level loss row as applicable.

### Manuscript reproduction

The legacy mixed-tip analysis is retained for traceability. Its RAxML-NG,
`chronos`, TimeTree-point, and CAFE 4.2 outputs must be labelled reproduction,
not evidence that the revised fossil or species-representation gates passed.

## Representative selection

Selection is frozen before looking at the desired topology. The ranking is:

1. verified biological identity, individual, accession, and matched annotation;
2. complete declared assembly/annotation scope;
3. genome and protein BUSCO plus annotation integrity;
4. contamination and sequence-name closure;
5. orthogroup occupancy and recoverable-locus count;
6. the predeclared continuity preference in
   `config/phylogeny/representation_policy.tsv`.

Topology, branch length, expansion count, or agreement with a preferred
biological story is never a representative-selection criterion. If the
continuity preference fails, the QC-ranked replacement and reason are recorded.

## One protein per gene

Every primary-tree terminal must provide one protein and one matching CDS per
selected gene. Newly downloaded and revised legacy inputs use
`scripts/qc/extract_primary_annotation.py`: all annotated isoforms are first
validated against the genome, then the longest valid spliced CDS is chosen;
ties are resolved by genomic span and transcript ID. Production additionally
requires exact independent gffread CDS/protein agreement. This is isoform
normalization of the published annotation, not de novo gene prediction.

The recovered manuscript script selected the longest genomic mRNA span per
gene. A legacy pair is reused only when the GFF contains one mRNA per gene or
that exact legacy selection plus protein--CDS translation has been verified.
`config/phylogeny/legacy_primary_isoform_policy.tsv` records every decision.
The current exceptions are explicit: *A. arguta* C retains its exact legacy
pair because revised extraction produced three gffread protein disagreements,
and *Vitis vinifera* retains its exact legacy pair because the legacy GFF3 has
duplicate `Source` attributes. Neither file is silently repaired.

For *A. latifolia*, the primary pair is re-extracted directly from the frozen
chromosome genome and GFF3. Of 41,317 annotated genes, 40,943 have a valid
coding transcript and pass gffread; 374 genes with no valid coding transcript
are omitted from the protein/CDS set and remain enumerated in the gene and
transcript audit tables. The previous header-only `.t` remap is diagnostic and
is not a production phylogeny input.

Specific rules are:

- *A. arguta*: complement C is preferred only if it passes; A-D and independent
  assemblies are diagnostic and representative-swap sensitivities.
- *A. deliciosa*: choose Qinmei or ADM as one parent bundle, then prefer D only
  if it passes. All A-F units remain diagnostic. Qinmei and ADM are never mixed
  within one primary representation.
- *A. eriantha*: HAP1/HAP2 are one individual. Prefer HAP1 only if both public
  haplotypes pass the declared gates; HAP2 and full annotated non-hybrid
  accessions are sensitivities.
- *A. macrosperma*: one unresolved available complement does not imply four
  recovered haplotypes.
- *A. rufa*: phylogeny selection is independent of gene-loss inclusion. Rank
  ARU, Fuchu, and the sequence-distinct ActinidiaBase v1 download by the same
  QC/occupancy rules. ActinidiaBase is not an exact sequence mirror of either
  named assembly, but its plant accession is unresolved, so it is an assembly
  sensitivity rather than a biological replicate. Retain both an omit-species
  sensitivity and end-to-end representative-assembly swaps.
- *A. x zhejiangensis*: retain A and B as two author-designated parental-lineage
  tips in the primary tree, CAFE table, per-lineage loss table, and PGLS. Their
  suffixes are upright. Do not add a shared-individual correlation constraint
  and do not require an A/B-pruned sensitivity.

## Production topology workflow

1. Validate one primary-isoform protein and matching CDS set for every selected
   representative and outgroup.
2. Run OrthoFinder 3 on the exact frozen proteome set. Its MSA/tree method and
   orthogroup/HOG level are recorded; DendroBLAST is reproduction-only.
3. For concatenation, retain strict one-copy orthogroups across every selected
   biological-species representative. Align proteins with MAFFT, back-translate
   through exact protein-CDS mappings, trim codons with one frozen occupancy
   rule, and reject ID or frame failures.
4. Infer each locus tree with IQ-TREE 2 under a declared ModelFinder and support
   policy. Contract or exclude weak branches only under a predeclared rule.
5. Infer the partitioned concatenation tree with IQ-TREE 2 and retain locus and
   site concordance diagnostics.
6. Infer the coalescent tree with ASTRAL-Pro from validated multi-copy
   gene-family trees and an exact copy-to-biological-species map. ASTRAL-Pro
   handles paralogy; it does not convert phased units into species replicates or
   solve reticulate polyploid history.
7. Compare topology, support, quartet score, occupancy, and representative-swap
   results. Freeze the accepted biological-species topology before dating.
   Genuine concatenation/coalescent disagreement is reported.

The current official tool references are the
[OrthoFinder documentation](https://orthofinder.github.io/OrthoFinder/),
[IQ-TREE documentation](https://www.iqtree.org/doc/),
[ASTRAL/ASTRAL-Pro repository](https://github.com/smirarab/ASTRAL), and
[CAFE 5 repository](https://github.com/hahnlab/CAFE5).

## Fossil-ready dating and CAFE 5

MCMCTree is the production dating engine. The topology is fixed first. A fossil
row may become active only when its required descendant and sister-lineage
bracketing taxa are present in the exact dated alignment and topology and have
passed asset/QC gates. There are currently no active fossil rows.

By explicit author decision on 2026-07-20, the revised dated tree therefore
uses four archived TimeTree 5 summary-confidence intervals as secondary soft
bounds. Each pair maps to a unique MRCA in the frozen 17-tip topology and its
raw API response and checksum are frozen. This result is eligible for the
declared downstream CAFE5 and PGLS analyses, but it must be called TimeTree
secondary-calibrated and never fossil-calibrated.

Neither acquired new outgroup currently activates a fossil. *Catharanthus* is
outside Rubiaceae and cannot bracket an internal Rubiaceae fossil. *Leea* is
outside Vitoideae, while the primary description places *Indovitis* within
Vitoideae; *Leea* plus the lone *Vitis* tip therefore cannot define crown
Vitoideae. These corrected gates are frozen in
`config/phylogeny/calibrations.tsv` and
`config/phylogeny/fossil_bracketing_taxa.tsv`.

Dating requires a prior-only run, at least two independent posterior chains,
convergence and effective-sample-size checks, leave-one-fossil-out analyses,
and every declared placement/bound sensitivity. `ape::chronos` with TimeTree
points remains a manuscript-reproduction comparison only. The recovered five
fixed points and exact artifact hashes are preserved in the two
`legacy_chronos_*` ledgers, but reproduction is blocked because the current
script input has 9 tips and the archived output has 17, and the original query
metadata/software version are missing. A secondary-TimeTree analysis is
allowed only through complete rows in
`config/phylogeny/secondary_timetree_constraints.tsv`; undocumented point ages
are prohibited. The current four rows are active under the author-approved
scope. See `config/phylogeny/dating_designs.tsv`.

CAFE 5 receives the matching dated biological-species tree and exactly one
representative-complement family count per species. Raw sums across A-F or
HAP1/HAP2 are prohibited. The primary result therefore estimates gene-family
dynamics of comparable representative complements, not organism-wide
polyploid dosage. Alternative representatives and any independently validated
phase-collapsed count policy are separate sensitivities. Taxon names and the
tree/count matrix must reconcile exactly. The predeclared Base Poisson and
Gamma3 attempts are retained even when one model fails initialization.

### Current validated execution

The frozen 17-tip IQ-TREE and ASTRAL-Pro topologies are identical (rooted
RF = 0). With no active fossil bracket, MCMCTree completed under the four
author-approved TimeTree 5 secondary intervals. Two 50,001-row posterior
chains passed the declared relaxed gate (minimum combined ESS 143.347;
maximum split-Rhat 1.034279), and the pooled mean tree is ultrametric with root
age 116.6885 Ma. This is TimeTree secondary-calibrated, not fossil-calibrated.

CAFE5 Base Poisson is the sole validated production family model. Of 28,128
prepared primary families, CAFE5's logged root-absence filter retained 15,066;
the fitted lambda is 0.0085698614157905 and 1,539 families have nominal
`p < 0.05`. Gamma3 failed during initialization and is frozen as unavailable.
It must not be retried through additional post hoc family removal or described
as supporting a Gamma model.

## PGLS status

PGLS completed only after the exact dated primary selected-lineage tree and
the callable-aware shared/non-shared lineage aggregation had both passed. It uses one
row per selected nonhybrid species plus separate *A. x zhejiangensis* A and B
rows and a matching ultrametric version of the accepted primary time tree.
HAP1/HAP2, *A. arguta* A--D, and *A. deliciosa* A--F do not otherwise become
independent observations.

For species `i`, let `C_i` be its definitive binary species calls (only
`positive_complete` or fully callable `not_positive`) and `S` the shared
positive-complete set across the frozen biological-species cohort. The primary
response numerator is the lineage-specific/non-shared positive count `L_i`;
the denominator is `|C_i \ S|`. Shared positives are excluded from both terms,
and partial, uncertain, not-called, or non-callable genes never enter the
denominator. The primary predictor is `log2_ploidy`.

The executed bundle includes the requested *A. rufa* exclusion and
leave-one-lineage-out fits. A full *A. rufa* assembly-swap remains a distinct
end-to-end sensitivity because it requires rebuilding the matched loss row,
topology, and time tree. Production release requires these declared
sensitivities: exclude *A. rufa*; replace
its accepted assembly with each other QC-passing sequence-distinct candidate
and rebuild the matching loss row, topology, and time tree; and omit every species
in turn while using the induced marginal covariance of the accepted rooted
time tree. Ordinary Gaussian PGLS is exploratory: unequal callable denominators
are not an observation-specific variance model, so publication remains blocked
until a denominator-aware phylogenetic count model is separately predeclared
and validated. See `docs/SPECIES_PGLS.md` for the executable contract.

## Resource ceiling

The total scientific worker pool across simultaneous OrthoFinder, MAFFT,
IQ-TREE, ASTRAL-Pro, MCMCTree, and CAFE processes must never exceed 15. Fifteen
is a project-wide ceiling, not permission to give each concurrent process 15
threads. Locus-level jobs must be scheduled so the sum of job count multiplied
by threads per job remains at most 15.
