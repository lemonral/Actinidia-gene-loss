# Archived species-level PGLS sensitivity

This sensitivity was completed as an exploratory analysis and was not used for
primary inference. The completed run retains the machine-readable status
`BLOCKED_DENOMINATOR_AWARE_MODEL_REQUIRED` to prevent an ordinary Gaussian PGLS
from being presented as a denominator-aware loss model. The accepted
17-tip protein/codon topology was dated with MCMCTree under four
checksum-bound TimeTree 5 secondary constraints because the compact cohort
has no active fossil bracket. The validated output is ultrametric and must not
be presented as fossil-calibrated dating. Its exact 13-lineage *Actinidia*
subtree, callable-aware species aggregation, and ploidy ledger passed the
input closure gate. Ordinary PGLS was enabled after production prerequisites
pass and completed with the required *A. rufa* exclusion and
topology-pruned leave-one-lineage-out diagnostics.

Execution is permitted only after all exact upstream products pass. The
mandatory upstream products are:

1. the dated, rooted, ultrametric primary biological-species tree;
2. the callable-aware shared/non-shared loss aggregation for the same species;
3. the validated biological-species ploidy ledger; and
4. a checksum-bound PGLS input-builder PASS report reconciling all three.

`enabled` does not mean that a provisional tree or an assembly-unit matrix may
be analysed. The ordinary PGLS tests whether `log2_ploidy` is associated with a
lineage-specific/non-shared loss trait after accounting for shared phylogenetic
history. It is retained as a phylogenetic trait sensitivity, not as a substitute
for a count model that respects unequal denominators.

## Statistical unit and aggregation

There is one row per selected nonhybrid species plus separate
*A. x zhejiangensis* A and B parental-lineage rows. HAP1/HAP2,
*A. arguta* A--D, and *A. deliciosa* A--F do not otherwise become independent
rows. No shared-individual constraint is added for the Zhejiang A/B rows and no
A/B-pruned sensitivity is required. Technical-unit columns such as
`assembly_unit_id`, `haplotype_id`, `subgenome_id`, `sample_id`, and
`terminal_id` remain rejected from the general species rows.

For biological species `i`, define:

- `S` as the reference genes classified `positive_complete` in every included
  biological species in the frozen cohort;
- `C_i` as the reference genes with a definitive binary species call in species
  `i`; and
- `L_i` as the genes in `C_i` with a positive lineage-specific/non-shared loss
  call for species `i`.

The denominator and raw response proportion are:

```text
D_i = |C_i \ S|
loss_proportion_i = |L_i| / D_i
```

`positive_complete` requires positive, callable evidence in every selected
haplotype or subgenome and the production aggregation rule is always
`all_units_positive`. `not_positive` requires callable retained evidence in
every selected unit. `positive_partial`, `uncertain`, non-callable, and
`not_called_loss` rows are excluded from `C_i`; they are not silently counted
as retained. Thus shared positives are excluded from both numerator and
denominator.

This estimand is **complete biological-species loss**, not the mean fraction of
lost haplotype/subgenome copies. It therefore differs from an assembly-unit or
copy-opportunity analysis and must be described that way in the manuscript.
If copy loss is the intended biological question, build and validate a separate
copy-opportunity model rather than changing this aggregation rule after seeing
the result.

The transformed response uses the Haldane--Anscombe correction:

```text
log((|L_i| + 0.5) / (D_i - |L_i| + 0.5))
```

The input builder must bind the exact schema-2.0 PASS aggregation manifest,
species matrix, shared-gene set, ploidy ledger, and ploidy PASS report used to
calculate every row. A historical denominator, an `any_unit_positive` rule, or
a denominator calculated before shared-loss removal is invalid.

## Required inputs and tree semantics

The headered TSV requires exactly one row per biological species:

```text
biological_species  analysis_level      loss_scope                  lineage_specific_nonshared_positive_loss_count  callable_denominator  log2_ploidy
Actinidia arguta    biological_species  lineage_specific_nonshared  731                                              32100                 2
```

Real fields are tab-separated. `analysis_level` must be
`biological_species`, `loss_scope` must be `lineage_specific_nonshared`, and
`callable_denominator` must be `D_i`. Counts are canonical decimal integers
satisfying `0 <= count <= callable_denominator`; scientific notation and
floating-point integer coercion are rejected. Values above the exact signed
64-bit limit are rejected, far above any plausible reference-gene universe.
At least six species are required
so each leave-one-species-out fit retains at least five.
`log2_ploidy` must be a canonical non-negative decimal in `[0,10]`; the input
builder recalculates it from an integer ploidy no greater than 1024.

The Newick file must be the accepted primary selected-lineage MCMCTree result
pruned to exactly the PGLS rows, including *A. x zhejiangensis* A and B. It must contain one rooted, strictly
bifurcating, ultrametric time tree. The top-level node is the accepted
biological-species MRCA root. Every non-root branch needs an explicit finite,
non-negative length. A checksum-bound time-tree PASS report identifies branch
lengths as million years, binds the dating manifest, and confirms the root and
exact tip set. Star trees and trees with a flat Pagel-lambda likelihood are
rejected because lambda is not identifiable.

## Model and command

The model has an intercept and the frozen primary predictor `log2_ploidy`.
Create that predictor from the passed biological-species ploidy ledger; do not
infer ploidy from assembly-unit labels.

First build the exact count table and its two PASS reports from the passed
species-loss bundle. The ploidy ledger must contain one ordered row per species
with `biological_species`, integer `ploidy`, `ploidy_source`, and
`source_reference` columns.

```bash
python scripts/downstream/build_species_pgls_input.py \
  --species-loss-dir results/species_loss/species_aggregation \
  --ploidy-ledger config/species_ploidy.tsv \
  --output-dir results/species_loss/pgls_input_bundle
```

This builder derives `D_i` only from non-shared `positive_complete` and
`not_positive` rows and fails if the species grid, shared intersection,
aggregation rule, ploidy set, or checksums do not reconcile.

```bash
python scripts/downstream/species_pgls.py \
  --data results/species_loss/pgls_input_bundle/pgls_input.tsv \
  --time-tree results/phylogeny/species_time_tree.nwk \
  --input-pass-report results/species_loss/pgls_input_bundle/pgls_input_pass.json \
  --species-loss-manifest results/species_loss/species_aggregation/species_loss_summary.json \
  --ploidy-ledger-pass-report results/species_loss/pgls_input_bundle/ploidy_ledger_pass.json \
  --time-tree-pass-report results/phylogeny/species_time_tree_pass.json \
  --predictor-column log2_ploidy \
  --sensitivity without_rufa='Actinidia rufa' \
  --output-dir results/statistics/pgls_ploidy
```

Brownian covariance is the shared root-to-MRCA branch length. Pagel's lambda
multiplies only off-diagonal covariance and is estimated by ML on `[0, 1]`.
The reported likelihood uses `RSS/n`. Coefficient SEs, 95% CIs, and two-sided
p-values use `RSS/(n-2)` and a t distribution with `n-2` degrees of freedom.

Those intervals and p-values are conditional on the fitted lambda. Ordinary
PGLS gives every lineage the same residual process regardless of its callable
denominator: the denominator changes the logit value but does not enter an
observation-specific variance model. This limitation is reported explicitly;
a denominator-aware phylogenetic count model may be added as a complementary
analysis without blocking the requested PGLS result.

## Mandatory sensitivities

1. **Exclude *A. rufa*.** The named `without_rufa` fit uses the primary data and
   the induced covariance for the remaining species.
2. **Swap the *A. rufa* assembly.** Replace the accepted representative with the
   other QC-passing sequence-distinct candidates; rerun loss calling and
   aggregation; rebuild each matching phylogeny and time tree; and run separate
   PGLS bundles. The unnamed ActinidiaBase plant remains an assembly sensitivity,
   not an independent species observation. Reusing the original loss row or tree
   is prohibited.
3. **Topology-pruned leave-one-species-out.** Omit each biological species and
   use the induced marginal covariance of the accepted rooted time tree. This
   exact covariance submatrix preserves the accepted root and branch depths; it
   deliberately does not reroot every sensitivity at a new MRCA.

Interpret exploratory output from effect size, interval width, lambda stability,
and sensitivity consistency, never from one p-value.

## Output and failure behavior

The output directory is published atomically and must not already exist.

- `analysis_data.tsv`: standardized input values and transformed response.
- `model_summary.tsv`: lambda, ML likelihood, variances, AIC, and sample size.
- `model_coefficients.tsv`: estimates and conditional exploratory inference.
- `fitted_residuals.tsv`: primary fitted values and residuals.
- `leave_one_species_out.tsv`: every omitted-species refit and prediction.
- `named_exclusion_sensitivities.tsv`: named exclusion-model results.
- `publication_gate.tsv`: upstream binding, tree, denominator, and model-fit gate.
- `analysis_manifest.json`: formulas, exact input bindings, runtime, and limits.
- `checksums.sha256.tsv`: checksums for every other bundle file.

The workflow fails if reports do not pass and bind exact bytes, tips differ,
the tree is non-bifurcating/non-ultrametric, lambda is unidentifiable, a matrix
or design is singular or numerically ill-conditioned, counts are invalid, or a
fit has too few species. Inputs
are analysed from immutable byte snapshots and rechecked before atomic
no-replace publication. A changed input or racing output therefore fails
without replacing anything.

The assembly-unit tree, manuscript `chronos` reproduction tree, and a CAFE tree
with another tip set are invalid PGLS inputs. The accepted
*A. x zhejiangensis* A/B primary-lineage tips are valid by the frozen author
decision. The real regression suite requires NumPy, SciPy, and Biopython. An
interpreter missing them produces a failing dependency-contract test; an
all-skipped PGLS log is never accepted as PASS.
