# Gene-loss classification and downstream analysis sets

## Reference set

The *Clematoclethra scandens* annotation contains 35,558 gene models. The
sequence-based loss analysis uses the 35,547 models that have matching CDS and
protein records in the frozen SynOrths reference bundle. The 11 BED-only models
without a matching record in that bundle are documented separately and are not
used for loss inference.

Each reference gene is evaluated independently in 23 *Actinidia* genomes.
This produces a complete grid of 817,581 genome-by-gene comparisons.

## Primary classification

The same rules are applied to every genome:

1. A reference gene with an exact SynOrths-supported ortholog is `retained`.
2. A missing-gene candidate with at least one genome-wide tBLASTX hit is
   `decayed` when the hit has identity of at least 50%, bit score of at least
   50, and E-value below `1e-5`. No alignment-length minimum is imposed.
3. A missing-gene candidate without a qualifying hit is `deleted`.
4. A comparison outside the resolved candidate scope is `not_called_loss`.

The positive-loss numerator is `decayed + deleted`. The resolved denominator
is `retained + decayed + deleted`, so `not_called_loss` is excluded from both.
The completed matrix contains 633,957 retained, 171,866 decayed, 7,961 deleted,
and 3,797 not-called comparisons.

`Decayed` is an operational sequence-similarity class. It indicates that a
homologous residual sequence was detected and does not by itself establish a
specific disruptive mutation or functional pseudogene.

## Coding-disruption evidence

Miniprot alignments are joined after primary classification. Qualifying
frameshift and internal in-frame-stop tags divide supported decayed calls into
frameshift-only, stop-only, and combined classes. Partial or terminally
truncated alignments are recorded as candidates when they meet their declared
coverage rules. These labels refine the evidence associated with a decayed
locus but do not change its primary `decayed` classification.

Within the primary decayed class, 19,888 calls form the strict pseudogenized
subset: 11,559 are frameshift-only, 3,258 are in-frame-stop-only, and 5,071
contain both signals. Candidate and unresolved decayed calls remain decayed and
are reported separately in the evidence catalogue.

The evidence is not treated as an exhaustive catalogue of mutation mechanisms.
Start-codon loss, splice-site disruption, exon deletion, gene fusion or fission,
transposable-element insertion, regulatory loss, and epigenetic silencing are
not inferred when the available inputs do not test them directly.

## Genome and lineage summaries

All 23 genomes remain separate in per-genome summaries and matched
within-species comparisons. A topology-only scaffold places units from the
same biological species as parallel unresolved tips. It is used to describe
species-specific, tree-node, recurrent, partial, and unresolved patterns and is
not an independently inferred 23-species phylogeny.

For branch-based functional and NLR analyses, the units are grouped into 13
biological lineages. Complete loss in a multi-unit lineage requires every
constituent unit to be `decayed` or `deleted`. Mixed retained and positive
states are classified as partial or homeolog-specific loss. Genes assigned to
a tree-node event on a focal root-to-lineage path are removed from that
lineage's risk set.

## Downstream analysis sets

- Loss counts, shared/non-shared summaries, and branch placement use the
  primary `decayed + deleted` positive state.
- Chromosome-position analyses use only `decayed` comparisons with an observed
  residual-sequence coordinate in the corresponding target assembly. Deleted
  and unlocalized comparisons are excluded. Shared and non-shared decayed
  comparisons are included together.
- Expression and gene-copy-number analyses also use `decayed` as the
  numerator, combine shared and non-shared comparisons, and exclude not-called
  comparisons from the resolved denominator. Expression is the arithmetic
  mean TPM across four *C. scandens* tissues. Gene copy number is the
  *C. scandens* CD-HIT 90% cluster size and is not target-genome copy-number
  variation.
- The primary functional analysis tests complete species-specific losses
  across 13 biological lineages. For each focal lineage, tree-node losses on
  the root-to-lineage path are excluded from the risk set. Complete
  species-specific lost genes form the foreground; the other annotation- and
  covariate-complete risk-set genes form the background. Logistic score tests
  account for the linear and quadratic effects of four-tissue mean expression
  and *C. scandens* gene copy number. Per-genome hypergeometric and
  topology-scaffold summaries are descriptive companion analyses, not
  independent species-level replicates.
- NLR analyses use the same primary loss matrix. Shared reference-NLR losses,
  non-shared unit-level calls, structural classes, complete target repertoires,
  and branch events are reported as distinct quantities with their matching
  denominators.

Some scripts and result-bundle directories retain `article_method` or
`manuscript_method` in their internal names for compatibility with validated
manifests. In public documentation, these stable identifiers refer to the
primary threshold-based classification defined above.
