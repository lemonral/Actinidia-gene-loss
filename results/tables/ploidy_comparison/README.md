# Ploidy-group comparison

`unit_loss_summary.tsv` is the 23-genome input used for Table S12. Each
haplotype, subgenome, or single-genome assembly is one statistical observation:
11 polyploid genomes and 12 diploid genomes. Positive loss is `decayed +
deleted`, and the resolved denominator is `retained + decayed + deleted`.

Regenerate the exact Mann-Whitney and mean-difference permutation tests with:

```bash
python scripts/statistics/compare_ploidy_loss_rates.py \
  --unit-summary results/tables/ploidy_comparison/unit_loss_summary.tsv \
  --output-dir results/tables/ploidy_comparison/generated
```

The two exact tests answer different questions. The Mann-Whitney test compares
the rank distributions of the two groups; the permutation test evaluates the
observed difference in arithmetic means. Both enumerate every possible group
assignment and are therefore deterministic.
