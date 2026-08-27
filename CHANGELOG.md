# Analysis changelog

## 2026-08-27

- Synchronized reader-facing terminology with the revised manuscript:
  `species-specific` and `tree-node` replace terminal/internal labels in final
  figures, while internal schema names remain unchanged for compatibility.
- Renamed the CD-HIT 90% cluster-size predictor to *C. scandens* gene copy
  number. Figure 5 values and fitted statistics are unchanged.
- Documented the primary functional foreground and background explicitly. For
  each of 13 biological lineages, complete species-specific losses form the
  foreground after tree-node losses on the root-to-tip path are removed; the
  remaining annotation- and covariate-complete risk-set genes form the
  comparison background.
- Updated final Figure 2, Figure 4, Figure 5, and Figure 6 renderers, including
  publication Arial validation and the revised Figure 4 lower-panel layout.
- Added a standalone, deterministic Table S12 implementation in which all 23
  genomes are independent observations (11 polyploid and 12 diploid). The
  primary exact Mann-Whitney result is `U = 22`, `P = 0.0056224567`; the exact
  mean-difference permutation result is `P = 0.0041380749`.
- Documented the 13-lineage NLR branch mapping as 188 events: 45
  species-specific and 143 tree-node events, including 138 events on the
  *Actinidia* stem.

## 2026-08-06

- Updated publication annotations in the expression and gene-copy-number
  figure. Statistical symbols are typeset consistently, the redundant
  `Global OLS` prefix is removed, and the panel-b annotation no longer overlaps
  the regression line.
- Updated the between- and within-chromosome decayed-locus figures to report
  very small likelihood-ratio-test probabilities as `P < 0.001` and to typeset
  the Benjamini-Hochberg q value consistently. Figure geometry and data are
  unchanged.
- Documented how the 23 genomes enter the negative-binomial chromosome
  models and how unit-level numerators and gene-opportunity denominators are
  pooled only for descriptive summaries.
- Documented the species-specific-loss logistic score tests, the expression and
  gene-copy-number covariates, and the deterministic selection of
  representative GO and KEGG terms.

## 2026-07-26 to 2026-07-28

- Replaced the expression input with the arithmetic mean TPM across four
  tissues and retained shared plus non-shared `decayed` calls in the primary
  expression and gene-copy-number summaries.
- Added decayed-only between-chromosome and five-zone within-chromosome models
  using target annotated genes as opportunity offsets.
- Added the 13-lineage species-specific-loss analysis with
  covariate-adjusted GO and KEGG score tests.
- Added JCVI support summaries between the 24-chromosome *C. scandens*
  reference and each accepted 29-chromosome *Actinidia* genome unit.
- Standardized publication labels by omitting `unphased` and the multiplication
  sign from *A. zhejiangensis* while retaining analytical identifiers.
