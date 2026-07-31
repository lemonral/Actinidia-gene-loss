# Decayed-only chromosome-position analysis

## Purpose

This analysis follows the article-method positional interpretation. Only
`decayed` calls are used as positive positional observations. A `decayed` call
has a qualifying residual sequence in the target assembly, so its observed
target coordinate can be tested. `Deleted` calls are excluded because absence
from the expected interval does not provide an observed target-genome locus.

All 23 assemblies, haplotypes, and subgenomes remain independent. No
species-level mean or consensus is used.

## Numerators and denominators

- Primary numerator: article-method `decayed` calls with a target-assembly
  coordinate.
- Primary denominator: annotated target genes in the matching assembly,
  chromosome, and, for within-chromosome tests, chromosome zone.
- Strict sensitivity: frameshift-supported, in-frame-stop-supported, and
  combined frameshift-plus-stop calls within the article-method `decayed`
  set. These subsets never replace or enlarge the primary numerator.
- Excluded: `deleted`, `retained`, `not_called_loss`, and spatially unlocalized
  `decayed` calls.

The exact target GFF and genome FASTA for every unit are registered in
`config/decayed_position_gene_denominators.tsv`. Chromosome identities use the
accepted HY4A-standardized `Chr01`--`Chr29` names.

## Five within-chromosome zones

Each locus is represented by its midpoint. Distance to the nearest chromosome
end is divided by half the chromosome length, giving a scale from 0 at either
end to 1 at the chromosome center. The five equal bins are therefore
orientation-independent:

1. `Z1_terminal`: outermost 10% from either chromosome end;
2. `Z2_subterminal`: 10--20% from either end;
3. `Z3_intermediate_outer`: 20--30% from either end;
4. `Z4_intermediate_inner`: 30--40% from either end;
5. `Z5_central`: the central 20% of the chromosome.

Target annotated genes are assigned to the same zones by their own midpoints.
Consequently, every positional rate is a loss count per matching target-gene
opportunity rather than a count per physical base pair.

## Statistical tests

Between-chromosome heterogeneity is tested with a negative-binomial model:

`decayed count ~ assembly unit + chromosome + offset(log(target gene count))`.

The chromosome term is evaluated by a likelihood-ratio test against the
matching reduced model. Chromosome-specific adjusted rate ratios are compared
with the adjusted grand mean, with Benjamini-Hochberg correction across the 29
chromosomes.

Within-chromosome heterogeneity is tested with:

`decayed count ~ assembly unit + chromosome + five-zone + offset(log(target gene count))`.

The five-zone term is evaluated by a likelihood-ratio test. Zone-specific
adjusted rate ratios use the central zone as the reference, with
Benjamini-Hochberg correction across the five contrasts. Opportunity-based
chi-square tests are retained as descriptive closure tests, including one
five-zone test for each chromosome.

## Production result

The checksum-closed production analysis is
`results/decayed_chromosome_distribution_v1_20260725` in the external data
store. It contains:

- 171,866 article-method `decayed` unit-gene rows;
- 138,599 spatially placed `decayed` rows used in positional numerators;
- 33,267 spatially unlocalized `decayed` rows reported but not forced into a
  position;
- 19,888 placed strict-pseudogenized rows within the primary `decayed` set;
- 971,672 target annotated genes across the 23 registered GFF files.

Chromosome identity has a strong primary effect
(`negative-binomial LRT chi-square = 1211.09`, 28 degrees of freedom,
`P = 2.50e-237`). The largest adjusted all-decayed burdens occur on Chr25
(242.20 per 1,000 genes; rate ratio 1.760), Chr24 (187.87; 1.365), and Chr22
(169.89; 1.234). The smallest occur on Chr16 (77.07; 0.560), Chr11 (87.20;
0.634), and Chr21 (102.64; 0.746).

The five-zone effect is also significant
(`negative-binomial LRT chi-square = 211.59`, 4 degrees of freedom,
`P = 1.21e-44`). Relative to the central zone:

- terminal: adjusted rate ratio 0.910, `BH q = 8.78e-7`;
- subterminal: 1.125, `BH q = 7.23e-10`;
- outer-intermediate: 1.154, `BH q = 9.44e-14`;
- inner-intermediate: 1.075, `BH q = 1.39e-4`.

The strict-pseudogenized sensitivity independently shows terminal depletion
(rate ratio 0.757, `BH q = 1.90e-15`); its other zones do not differ
significantly from the central zone after correction. This agreement supports
the broad positional trend, while the primary inference remains based on all
article-method `decayed` calls.

## Interpretation limits

An ordinary `decayed` residual can be local, displaced on the same chromosome,
or located on another chromosome. It is evidence for an observed homologous
residual, not proof of a particular rearrangement mechanism. The current
analysis does not claim centromere association, transposable-element causation,
or a specific inversion/translocation breakpoint without independent
annotations. It tests reproducible chromosome and within-chromosome
heterogeneity in the article-method decayed signal.

