# Frozen decision rules

`config/analysis_parameters.toml` records the methodological choices that must
not drift between candidates. The numerical JCVI gate reproduces the submitted
study's 50% *Clematoclethra scandens* reference-side coverage criterion, but the
rebuild also reports the query-side result and all denominators. Candidate
assemblies are never selected from BUSCO or JCVI alone: identity, matched file
scope, annotation compatibility, unplaced sequence content, BUSCO, JCVI,
SynOrths, and callable loss evidence are considered together in the assembly
decision ledger.

## Assembly candidates

- The 2026 *Actinidia eriantha* HAP1 and HAP2 releases are two haplotypes from
  one biological individual. Both are evaluated, and neither is treated as an
  independent species replicate.
- The Qinmei 2025 and ADM 2026 *A. deliciosa* bundles are independent candidate
  replacements. Genome, GFF, CDS, and protein sequence-name rules must all pass
  before an A-F split is accepted. The better matched candidate is selected;
  the other remains in the QC supplement.
- The ARU, complete Fuchu, and ActinidiaBase *A. rufa* bundles are all evaluated.
  The submitted chromosome-only Fuchu subset remains a labelled reproduction
  input. A complete candidate may replace it only after the same JCVI and
  callable-scope gates used for every other taxon.
- Every candidate, including a later-excluded assembly, remains in the assembly
  and annotation QC table.

## Shared and non-shared loss

Assembly units are aggregated into biological species before the reviewer's
shared-loss question is answered. `positive_complete` means that every selected
unit for that species supports a positive loss. A shared loss is
`positive_complete` in every included biological species. Partial haplotype or
subgenome evidence is retained as `positive_partial`; it is not a shared loss
and cannot by itself be called a lineage-restricted species loss.

`not_called_loss` means only that a gene was absent from an earlier positive
list. It is not retained evidence and is aggregated as `uncertain`.
In the revised primary matrix, `deleted` and strictly supported
`pseudogenized` are the positive unit states. `Deleted` requires a callable
bilateral interval and no qualifying local alignment. `Pseudogenized`
requires an explicit Miniprot frameshift or in-frame stop plus the strict
coverage, identity, and score gate. A local translated alignment without that
event remains `uncertain`; both uncertain states are excluded from resolved
rate denominators.
The historical `decayed=pseudogenized` interpretation is preserved only in a
separately labelled manuscript-era reproduction matrix.

The shared set is reported separately. Expression, copy-number, spatial, NLR,
and comparative functional analyses use the non-shared set unless their output
is explicitly labelled as a submitted-analysis reproduction.

## Chromosome naming

For every newly downloaded or changed 29-chromosome *Actinidia* unit,
chromosome names are harmonized to Hongyang v4 HY4A by a global one-to-one
maximum-nucleotide-similarity assignment. The result must contain each
canonical label `Chr01`--`Chr29` exactly once. Absolute alignment support and
independent HY4P/JCVI evidence are retained as QC diagnostics, but they do not
block naming under the author's simplified naming rule. Publisher chromosome
direction is preserved: chromosome sequences are not reverse-complemented and
GFF coordinates or strands are not flipped merely to match HY4A orientation.

The relabelled FASTA and GFF must remain a matched bundle. Sequence names in
FASTA records, GFF sequence-region directives, and all GFF feature rows are
changed by the same bijection. CDS and protein sequences must close exactly
against the relabelled genome/GFF before the bundle is accepted for JCVI,
gene-loss, or position analyses.

## Spatial interpretation

The primary position analysis uses only observed target-genome loci: a strict
`pseudogenized` disrupted-alignment midpoint is compared with the observed
target midpoint of an exact-SynOrths `retained` gene. It uses mutually
exclusive equal-width bins and this observed-locus gene-opportunity
denominator. A positive deletion has no observed target-gene feature; its
callable expected-locus midpoint, delimited by bilateral
same-target-chromosome SynOrths anchors, is therefore used only in an
explicitly labelled sensitivity. It is never described as an observed remnant
coordinate.

Distance to the nearest chromosome end is a valid coordinate-derived measure
and is normalized from zero at an end to one at the chromosome centre. It must
not be called telomere distance unless an independent telomere interval is
supplied. Likewise, failed or exploratory centromere predictions are not
evidence; centromere distance is emitted only for independently supported
intervals. The manuscript-era nested intervals are available only as a clearly
labelled sensitivity reproduction.

## Terminology

`assembly unit` is used for haplotypes, subgenomes, and unphased assemblies.
`biological species` is used for species-level inference. `terminal` is reserved
for a leaf in a phylogenetic tree and is not used in gene-loss, QC, JCVI,
expression, copy-number, spatial, or NLR results.
