# Publication figure policy

Only the latest production figures are maintained. Historical, diagnostic, and
superseded figures remain outside the publication bundle.

Display labels do not alter internal identifiers or analysis groupings:

- omit `unphased`;
- display the two parental lineages as *A. zhejiangensis* A and
  *A. zhejiangensis* B, without a multiplication sign;
- italicize Latin binomials and keep informative assembly suffixes upright;
- use `(a)`, `(b)`, `(c)`, and `(d)` for panel labels;
- omit panel titles when the axis labels, legend, and caption identify the
  content.

The same display policy also omits the technical suffixes `ActinidiaBase v1`
and `unresolved polyploid unit`. These changes affect labels only; analytical
identifiers and metadata remain unchanged.

The four-panel genome-evolution overview is assembled from completed figures
by `scripts/figures/compose_genome_evolution_overview.py`. The compositor only
trims white margins, rescales panels proportionally, and adds panel labels. It
does not rerun or modify any scientific analysis. The source images are passed
on the command line because the Circos and Ks panels are maintained outside
this repository.

It can be regenerated without rerunning scientific analyses:

```bash
python scripts/figures/compose_genome_evolution_overview.py \
  --circos /path/to/circos.png \
  --orthofinder results/figures/orthofinder_species_profiles/gene_category_composition.png \
  --ks /path/to/ks_kde_mountains3.png \
  --phylogeny results/figures/phylogeny_cafe/primary_phylogeny_cafe.png \
  --output-dir results/figures/genome_evolution_overview
```

Only input basenames, byte counts, and SHA-256 checksums are recorded in the
publication bundle.
