# Manuscript-method loss classification and refinement

## Primary article-comparable layer

The article-comparable trend is calculated independently for all 23 assembly
units. It reproduces `lost_type_83.py` exactly at the biological-rule level:

1. an exact SynOrths anchor is `retained`;
2. among missing-gene candidates, any genome-wide tBLASTX hit with percent
   identity at least 50, bit score at least 50, and e-value below `1e-5` is
   `decayed`;
3. there is no alignment-length minimum;
4. a candidate without such a hit is `deleted`;
5. a row outside the historical candidate scope is `not_called_loss`, not a
   forced deletion.

This gives a complete 23 x 35,547 grid with 633,957 retained, 171,866 decayed,
7,961 deleted, and 3,797 not-called rows. The manuscript-positive trend is
`decayed + deleted` over `retained + decayed + deleted` for each assembly unit.

The primary evidence bundle additionally keeps all 23 units as separate
terminals on a topology-only scaffold. It contains 179,827 positive unit–gene
rows, 3,616 reference genes positive in all 23 units, and 56,602 maximal
loss-event rows. This scaffold expands units within a biological species as
parallel terminals and is not a newly inferred 23-species phylogeny.

The word `decayed` means sequence detected by the article threshold. It is not
synonymous with experimentally validated pseudogenization.

## Orthogonal cause refinement

Already-completed Miniprot evidence is joined to the article layer without
rerunning a search and without changing `retained/decayed/deleted` totals.
Explicit frameshift and in-frame-stop tags are separated into frameshift,
stop, and frameshift-plus-stop evidence. Tags below the strict quality gate are
kept as candidates. For qualifying local alignments, missing at least 20% of
the reference protein from an alignment end is reported only as an N-terminal,
C-terminal, or both-terminal truncation candidate.

The resulting cause labels are deliberately evidence labels, not an exhaustive
catalogue of mutation mechanisms. Current inputs do not reliably establish
start-codon loss, terminal-stop loss, splice-site disruption, exon loss, gene
fission/fusion, TE insertion, regulatory loss, or epigenetic silencing. Those
mechanisms are not guessed from low coverage alone.

## Use in downstream trend analyses

Expression, reference-family copy class, and assembly-unit loss summaries may
use the article-comparable `decayed + deleted` trend. Cause-refined categories
must be displayed separately as supporting subdivisions. They must not be used
to silently change the article totals or to imply that every `decayed` locus is
a pseudogene.

The strict callable-aware matrix remains a sensitivity analysis for testing
whether the article-level trends persist under stronger evidence requirements.

The completed article-method downstream bundle retains all 23 assembly units.
It identifies 3,616 reference genes positive in all 23 units and excludes this
shared set from the main non-shared expression/copy panels. A biological
lineage is called complete loss only when every assembly unit assigned to that
lineage is positive. Mixed positive/retained states in *A. arguta*, *A.
chinensis*, *A. deliciosa*, or *A. eriantha* are partial/homeolog-specific,
not complete species loss.

The primary functional analysis does not aggregate assembly units. It contains
23 independent foregrounds, each comprising every `decayed + deleted` gene in
one assembly unit. It does not require that any other unit or lineage be
retained. The matching background is `retained + decayed + deleted` in that
same unit, with `not_called_loss` excluded. This gives 179,827 unit–gene
foreground memberships and 6,420 BH-significant GO/KEGG rows across the 23
units.

The publication-oriented functional detail figure names representative GO
biological-process, molecular-function, cellular-component, and KEGG pathway
terms. Point area is the actual lost-gene count and color is
`-log10(BH q-value)`. Eligible terms require BH `q <= 0.05`, at least two
foreground genes, and fold enrichment greater than one; ranking is frozen by
significant-unit recurrence, median significance, median fold enrichment,
total contributing genes, and stable term identifier. The earlier
significant-term count heatmap is retained as QC rather than the main
biological display.

The matching scaffold-aware functional layer tests all 39 internal-node,
within-species unit-group, and unit-terminal event sets. Its 56,602
event–gene memberships are tested against the 33,998 reference genes with
resolved article-method states in all 23 units. This produces 905
BH-significant GO/KEGG rows. Nodes without a significant term remain explicit
zero rows. This layer describes functions associated with event placement; it
is not a PGLS and does not turn the topology-only scaffold into a newly
inferred 23-species phylogeny.

For supplementary evolutionary interpretation, complete lineage losses are
also mapped to the matching frozen 13-lineage topology under a minimum
irreversible-loss description. A monophyletic complete-loss set is one
terminal or internal branch event; disjoint complete-loss clades are repeated
independent events. Genes with a not-called lineage are excluded from exact
branch placement. This tree-aware layer is not a PGLS and does not replace the
primary 23-unit functional analysis.

The pure species-specific sensitivity applies a much narrower foreground:
one exact terminal-branch loss, resolved retained states in all other
lineages, and no recurrent, internal, partial/homeolog-specific, or unknown
pattern. It is retained only as a supplementary sensitivity because it does
not represent the broader unit-level loss definition requested for the main
analysis. This produces 1,167 non-overlapping species-specific genes. GO
biological process, molecular function, cellular component, KEGG orthology,
and KEGG pathway are corrected separately within each species foreground.

The NLR analysis reuses the completed NLR-Annotator repertoire
calls. Reference-NLR loss uses `decayed + deleted` over
`retained + decayed + deleted` at the 23-unit level, without species
aggregation. Reference NLRs positive in all 23 units are excluded from the
non-shared set. Positive calls are then subdivided without changing that
numerator into no qualifying translated hit, frameshift-supported,
in-frame-stop-supported, combined frameshift-and-stop, truncation/partial-
alignment candidate, and residual-sequence mechanism-unresolved groups. The
strict callable-aware NLR result remains a separately labelled sensitivity.
As an orthogonal classification, every positive non-shared reference-NLR call
also inherits the structural class of its exact *C. scandens* NLR-Annotator
record: CC-NBARC, CC-NBARC-LRR, NBARC, NBARC-LRR, TIR, TIR-LRR, TIR-NBARC, or
TIR-NBARC-LRR. This classifies which NLR architecture was lost; it does not
replace the independent loss-evidence mechanism groups and does not infer a
complete target-gene architecture from a fragmented residual sequence.
The 138 unique reference NLRs lost in all 23 units are shown as a separate
shared-loss composition. Complete per-unit repertoires are counted directly
from each finished NLR-Annotator call set; these additionally retain the rare
TIR-CC-NBARC-LRR class observed in two target calls.
Class-specific loss percentages use reference-NLR opportunities rather than
the extant target repertoire as the denominator. For each unit and class, the
total rate is `(shared loss + non-shared loss) / (shared reference
opportunities + resolved non-shared reference opportunities)`. The
non-shared-only rate is `non-shared loss / resolved non-shared reference
opportunities`. Not-called comparisons are excluded from both denominators.
This avoids comparing absent reference orthologs directly with lineage-specific
target repertoire expansions.

## Mechanism-stratified target-coordinate analysis

Spatial analysis uses residual sequence coordinates in each target
*Actinidia* assembly. Chromosome names are converted to the corresponding HY4A
`Chr01`--`Chr29` label before pooling or comparison. Ordinary `decayed` calls
are included when an already-computed qualifying residual alignment provides a
coordinate. Calls without a residual coordinate remain unlocalized and are not
assigned an expected-interval midpoint.

The six evidence groups are analyzed separately. A best genome-wide residual
placement outside the expected interval is labelled a same-chromosome or
interchromosomal displacement candidate, not a confirmed inversion,
translocation, or orthologue, because paralogous placement remains possible.
Chromosome heterogeneity is evaluated against chromosome-length opportunities
within each assembly unit. The audited *C. scandens* GFF contains no repeat or
transposon features, so a TE-association analysis is unavailable rather than
being inferred from gene annotations.
