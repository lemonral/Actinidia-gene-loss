# Analysis changelog

## 2026-08-06

- Updated publication annotations in the expression and reference-family-size
  figure. Statistical symbols are typeset consistently, the redundant
  `Global OLS` prefix is removed, and the panel-b annotation no longer overlaps
  the regression line.
- Updated the between- and within-chromosome decayed-locus figures to report
  very small likelihood-ratio-test probabilities as `P < 0.001` and to typeset
  the Benjamini-Hochberg q value consistently. Figure geometry and data are
  unchanged.
- Documented how the 23 genome units enter the negative-binomial chromosome
  models and how unit-level numerators and gene-opportunity denominators are
  pooled only for descriptive summaries.
- Documented the species-complete terminal-loss logistic score tests, the
  expression and family-size covariates, and the deterministic selection of
  representative GO and KEGG terms.

## 2026-07-26 to 2026-07-28

- Replaced the expression input with the arithmetic mean TPM across four
  tissues and retained shared plus non-shared `decayed` calls in the primary
  expression and reference-family-size summaries.
- Added decayed-only between-chromosome and five-zone within-chromosome models
  using target annotated genes as opportunity offsets.
- Added the 13-lineage species-complete terminal-loss analysis with
  covariate-adjusted GO and KEGG score tests.
- Added JCVI support summaries between the 24-chromosome *C. scandens*
  reference and each accepted 29-chromosome *Actinidia* genome unit.
- Standardized publication labels by omitting `unphased` and the multiplication
  sign from *A. zhejiangensis* while retaining analytical identifiers.
