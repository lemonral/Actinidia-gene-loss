# Tree-aware functional analysis of article-method loss

This analysis adds phylogenetic interpretation after, not inside, the
article-comparable loss classification. It never changes a unit-level
`retained`, `decayed`, or `deleted` call.

## Lineage states

For each reference gene, the 23 assembly units are grouped into 13 biological
lineages only for tree interpretation:

- complete lineage loss: every assigned unit is `decayed` or `deleted`;
- partial/homeolog-specific loss: at least one assigned unit is positive and
  at least one is retained;
- retained: at least one retained unit and no positive unit;
- unknown: no retained unit and at least one not-called unit.

Partial loss is not treated as absence of the entire lineage. This distinction
is essential for the four-unit *A. arguta*, two-unit *A. chinensis*, six-unit
*A. deliciosa*, and two-unit *A. eriantha* groups.

## Branch events

Only genes with states resolved in all 13 lineages enter exact branch
placement. Complete-loss leaves are covered by the minimum set of maximal
all-loss clades on the exact TimeTree-secondary-calibrated topology. One
single-leaf clade is a terminal-branch event; one multi-lineage clade is an
internal-branch event; two or more disjoint clades are repeated independent
events. Loss in all 13 lineages maps to the Actinidia stem relative to the
retained *Clematoclethra scandens* reference.

Of 35,547 reference genes, 34,637 have complete lineage-state information and
910 retain missing data. Exact patterns contain 3,280 single terminal-branch,
3,846 single internal-branch, and 5,620 repeated-independent-loss genes. A
further 8,472 genes show partial lineage loss without complete lineage loss.

## GO and KEGG

GO, KEGG KO, and KEGG pathway foregrounds are built separately for each exact
branch and for the broader pattern categories. Tests use the frozen reference
eggNOG-mapper annotation, one-sided hypergeometric over-representation, and
Benjamini-Hochberg correction independently within each foreground and
ontology. Branch foregrounds use the 34,637 fully resolved genes as their
background. Unit-level all-shared and any-nonshared summaries use the complete
reference universe.

The complete analysis has 31 foregrounds and 1,826 rows at BH-adjusted
`q <= 0.05`. Nested categories are intentionally not mutually exclusive; for
example, the 3,616 genes positive in all 23 units are also part of the
Actinidia-stem branch foreground. Enrichment is homology-based functional
interpretation, not experimental validation of loss mechanism.

PGLS is not used in this analysis. PGLS addresses association between lineage
traits and loss rates, whereas this workflow asks where complete losses are
most parsimoniously placed and what annotated functions occur in those branch
or pattern sets.
