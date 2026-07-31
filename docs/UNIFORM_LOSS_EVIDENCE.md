# Uniform old/new loss evidence

The production V2 analysis applies one evidence chain to all 23 assembly
units.  It does not copy the manuscript-era `decayed` label into the primary
matrix.

## Evidence states

| State | Exact rule | Positive loss? |
|---|---|---:|
| `retained` | The reference gene has an exact accepted SynOrths target anchor. | No |
| `deleted` | Two unambiguous same-target-chromosome SynOrths flanks define a callable interval, but no local Miniprot alignment passes query coverage 0.50, exact identity 0.50, and alignment score 50. | Yes |
| `pseudogenized` | A qualifying local Miniprot alignment contains `fs>0` or `st>0` and also passes the strict disruption gate: query coverage at least 0.80, exact identity at least 0.70, and alignment score at least 100. | Yes |
| `uncertain, callable=true` | A reliable interval exists, but local sequence has no supported disruptive event or an event falls below the strict disruption-quality gate. | No |
| `uncertain, callable=false` | No reliable bilateral interval/evidence chain exists. | No |

`fs` and `st` are Miniprot's explicit frameshift and in-frame-stop tags.
Noncanonical splice sites alone are not positive.  Because no matched raw-read
validation is available for every assembly, `pseudogenized` means
assembly-sequence-supported disruption; it is not claimed to be a
population-fixed allele.

## Chromosome-position analysis

The primary spatial numerator is strict `pseudogenized` at its observed target
alignment locus. Its denominator is the same unit's observed loci: strict
`pseudogenized` plus exact-SynOrths `retained`. Both uncertain states and
`deleted` are excluded from this primary model. A positive deletion has no
observed target-gene feature, so its bilateral-anchor expected-interval
midpoint appears only in explicitly labelled deleted-only and combined
sensitivity models.

Position is normalized from chromosome end (0) to chromosome center (1).
Binomial logistic models include assembly-unit fixed effects, reference-gene
clustered sandwich errors, a legacy/new slope interaction, and separate
legacy/new sensitivities.  A callable-candidate positive-versus-uncertain
diagnostic reports whether censoring changes with position.  Centromeres are
not analyzed without independent centromere intervals.

## Historical reproduction

The original manuscript-era `decayed=pseudogenized` tables remain separate
historical evidence.  They are never merged into the V2 primary numerator.

## Frozen V2 outputs

The complete matrix contains 23 x 35,547 = 817,581 unit-gene rows: 57,027
`deleted`, 20,046 strict `pseudogenized`, 633,957 `retained`, and 106,551
`uncertain`.  The 13-lineage aggregation contains 287 genes that are
positive-complete in every lineage.

The primary chromosome-position model uses 625,857 observed target loci
(20,046 strict `pseudogenized` and 605,811 retained). With normalized end
distance increasing from chromosome end (0) to center (1), the pooled strict
pseudogene slope is 0.339 log-odds (odds ratio 1.404; reference-gene clustered
Wald p = 2.25e-6). Source interaction and source-stratified estimates remain in
the production report. The deleted-only and combined expected-locus models,
plus the position-dependent uncertainty diagnostic, are sensitivity analyses
and not the primary target-locus result.

The matching ordinary PGLS was also refit. Its log2-ploidy coefficient is
-0.0832 (two-sided p = 0.628); excluding *A. rufa* gives -0.0406 (p = 0.797).
This remains exploratory because an ordinary Gaussian PGLS does not model the
unequal precision induced by different resolved denominators.

## Expression and copy number

Both Figure 3 panels use the unified non-shared matrix, not the historical
`decayed=pseudogenized` result. Removing the 287 shared genes leaves 35,260
reference genes. Across 23 units, 704,429 `retained`, `deleted`, or strict
`pseudogenized` rows form resolved opportunities; 106,551 `uncertain` rows are
excluded. In every rate, the numerator is `deleted + pseudogenized` and the
denominator is `retained + deleted + pseudogenized`.

The expression panel uses the declared *C. scandens* leaf S23I0033
featureCounts raw-count column. The 35,260 non-shared reference genes are
ranked after `log(raw count + 0.1)` and split deterministically into 14 nearly
equal bins. A separate loss rate is calculated for every assembly unit and
bin; points are coloured by ploidy. This is reference-gene expression, not
expression measured independently in every target species.

The copy-number panel uses CD-HIT 0.90 clusters of the complete reference
protein set. A gene's copy number is its original reference-cluster size,
assigned before shared genes are removed. Classes 1--7 enter the fitted panel
because each contains more than 100 reference genes; larger classes are shown
only in QC. This is a reference gene-family similarity/redundancy measure, not
a target-assembly CNV estimate. Rows without a mapped reference cluster and
all uncertain calls are excluded rather than converted to retained.

## Functional enrichment

The article-comparable primary functional analysis is unit-resolved and is
documented in `MANUSCRIPT_METHOD_LOSS_CLASSIFICATION.md`. It uses 23
independent `decayed + deleted` foregrounds without species or subgenome
aggregation.

This stricter sensitivity uses the unified evidence
semantics. Its non-shared foreground is `deleted + pseudogenized`
after excluding the 287 shared-positive genes; strict pseudogenized-only is a
separate sensitivity. Uncertain and partial calls are excluded. One-sided
hypergeometric tests use annotated reference genes as the opportunity
background, followed by Benjamini-Hochberg correction within foreground and
ontology.

The pooled non-shared primary foreground has 29 significant GO terms, 27 KEGG
KOs, and 30 KEGG pathways at adjusted `q < 0.05`. The strict-pseudogenized-only
sensitivity has 18, 21, and 32, respectively. The shared foreground has no
significant term. Annotation inputs and ontology/description files are
checksum-bound; older GO identifiers absent from the frozen ontology are
counted as unresolved rather than guessed. These results are homology-based
functional hypotheses and do not by themselves establish mechanism.

The category summary separates GO biological process, molecular function,
and cellular component, plus KEGG orthology and KEGG pathway. It reports the
pooled and all 13 lineage foregrounds in separate combined-loss and strict-
pseudogenized-only panels. For the pooled combined set the five category counts
are 12, 17, 0, 27, and 30; for the strict-pseudogenized-only sensitivity they
are 4, 14, 0, 21, and 32. These are counts of significant terms, not effect
sizes, so annotation coverage remains in the plot-data table.

Publication plot labels are shortened independently of analysis semantics:
*A. zhejiangensis* A/B, *A. rufa*, and *A. macrosperma*.
Internal unit IDs, lineage grouping, and every numerator/denominator remain
unchanged.
