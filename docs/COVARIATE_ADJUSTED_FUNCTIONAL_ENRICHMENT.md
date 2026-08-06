# Covariate-adjusted functional enrichment

## Analysis set and terminal events

The analysis uses 33,998 *Clematoclethra scandens* reference genes with resolved
states across all 23 genome units. The units are grouped into 13 biological
lineages for branch analysis. A multi-unit lineage is completely lost only when
every constituent haplotype or subgenome is `decayed` or `deleted`. Genes
assigned to an ancestral event on a lineage's root-to-tip path are removed from
that terminal risk set. The terminal foreground therefore contains complete
losses first assigned to the focal terminal lineage.

Reference functions are taken from the checksum-bound eggNOG-mapper annotation
generated with the Viridiplantae taxonomic scope. GO roots and unresolved GO
identifiers are excluded. GO biological process, molecular function, and
cellular component are retained as separate namespaces for interpretation.
KEGG orthology and KEGG pathway membership form two additional annotation
systems.

## Covariate-adjusted tests

For each lineage and annotation system, the background contains risk-set genes
with complete covariates and at least one annotation in that system. The null
logistic model contains z-standardized log2(four-tissue arithmetic mean TPM +
0.1), its squared term, z-standardized log2(reference CD-HIT 90% gene-family
size), and its squared term. Twenty-four genes without a family-size estimate
are excluded, leaving 33,974 covariate-complete genes before lineage- and
annotation-specific filtering.

Each term is tested separately with a binary membership indicator. The
one-sided efficient score statistic compares the observed loss residual among
term members with that expected under the null model. Terms with fewer than
five background genes or fewer than two terminally lost genes are not tested.
Benjamini-Hochberg correction is applied within each lineage and annotation
system. The three GO namespaces form one GO testing family. KEGG orthology and
KEGG pathway form separate testing families. A positive score-test association
with `q <= 0.05` is refitted with the full logistic model to estimate the
adjusted odds ratio and 95% confidence interval. Non-convergent full refits are
retained in the complete table but excluded from figures.

The validated species-complete analysis contains 13 terminal lineages, 19,192
terminal event memberships, 32,591 tested lineage-term rows, and 3,646
score-test significant rows. The significant rows comprise 3,143 GO, 390 KEGG
orthology, and 113 KEGG pathway associations.

## Representative terms in figures

Figure selection does not change the complete statistical results. Eligible
terms must have a significant positive association and a converged full-model
refit in at least one lineage. Terms are ranked by the number of significant
lineages, then by the smallest adjusted q value, and then by the mean adjusted
log2 odds ratio across significant lineages. Term identifiers provide a stable
final tie-break.

GO redundancy is reduced using the registered ontology graph. If two terms have
an ancestor-descendant relationship and the Jaccard similarity between their
sets of significant lineages is at least 0.65, only the higher-ranked term is
retained. The current figures show up to 12 biological-process terms, 12
molecular-function terms, 10 cellular-component terms, 12 KEGG orthology terms,
and 10 KEGG pathway terms. Point size represents score-test significance and
colour represents the full-model adjusted log2 odds ratio. All tested terms are
retained outside the display subset.
