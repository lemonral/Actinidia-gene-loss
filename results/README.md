# Curated results

Only validated small summaries, plot data, and final figures are stored here.
Raw and intermediate outputs remain in the external data store. Every figure
directory contains the plotted table and validation/checksum metadata, and no
file in this tree contains a private runtime path.

Canonical revised outputs currently include:

- `tables/manuscript_method_loss`: 23-unit reproduction of the article's
  SynOrths plus genome-wide tBLASTX `retained/decayed/deleted` rule, followed by
  a separate conservative cause refinement that never rewrites the article
  class;
- `tables/manuscript_method_downstream`: the matching 23-unit article-method
  shared/non-shared tables, 13-lineage complete/partial states, and exact
  terminal/internal/repeated-loss event assignments on the frozen topology;
- `tables/unit_functional_enrichment_article_method`,
  `figures/unit_functional_enrichment_detail`, and
  `figures/unit_functional_enrichment_article_method`: primary GO/KEGG
  enrichment for 23 independent `decayed + deleted` foregrounds. The detailed
  figure names representative terms and encodes contributing lost-gene counts
  and BH significance; the older category-count panel is retained as QC.
  Haplotypes and subgenomes remain separate, with no retained-state
  requirement in other units;
- `tables/species_specific_functional_enrichment` and
  `figures/species_specific_functional_enrichment`: GO/KEGG results for the
  1,167 pure single-terminal species-specific losses, retained as a narrow
  supplementary sensitivity;
- `tables/unit_loss_evidence_scaffold` and
  `figures/loss_evidence_classification`: the 23 separate historical-threshold
  unit counts, per-unit decayed frameshift/stop classes, all-23-unit
  shared/non-shared partition, and exact event placement on a topology-only
  assembly-unit scaffold; no haplotype or subgenome is collapsed;
- `tables/scaffold_functional_enrichment` and
  `figures/scaffold_functional_enrichment`: GO/KEGG enrichment for all 39
  internal-node, unit-group, and unit-terminal event sets on that scaffold,
  using the 33,998 genes resolved across all 23 units as a common background;
  all true zero-significance nodes remain visible;
- `tables/loss_evidence_branch_summary` and
  `figures/tree_branch_loss_evidence`: the earlier 13-lineage aggregation,
  retained as a supplementary comparison rather than the primary unit-level
  event display;
- `figures/manuscript_method_unit_loss`: the article-method main loss-count
  composition for all 23 units;
- `figures/expression_copy_loss_all_decayed`: pooled shared + non-shared
  `decayed`-only relationships with reference leaf expression and reference
  protein-family size for all 23 independent units;
- `tables/nlr_loss_types` and `figures/nlr_loss_types`: complete NLR
  repertoires and 23-unit non-shared NLR loss partitioned into six mutually
  exclusive evidence groups, without species aggregation;
- `tables/nlr_loss_structural_classes` and
  `figures/nlr_loss_structural_classes`: the same 254 non-shared
  `decayed + deleted` NLR calls classified by structural class, a separate
  138-gene shared-loss column, and the complete class-resolved NLR-Annotator
  repertoire for every assembly unit; panel C reports class-specific loss
  percentages using resolved reference-NLR opportunities;
- `tables/loss_mechanism_spatial` and `figures/loss_mechanism_spatial`:
  target-assembly residual coordinates for all positive calls that can be
  placed, split by six evidence groups and harmonized to HY4A
  `Chr01`--`Chr29`;
- `tables/tree_functional_enrichment` and
  `figures/tree_aware_functional_summary`: GO/KEGG results classified by
  complete terminal-branch loss, ancestral internal-branch loss, repeated
  independent loss, shared loss, and partial/homeolog-specific loss in the
  earlier 13-lineage aggregation, retained as a supplementary comparison;
- `figures/assembly_annotation_qc`: 23-unit assembly, annotation, and BUSCO QC;
- `figures/chromosome_similarity_naming`: complete one-to-one Chr01--Chr29
  naming relative to HY4A, with nucleotide and JCVI support shown as QC only;
- `tables/uniform_loss` and `tables/species_loss_uniform`: the stricter 23-unit
  sensitivity taxonomy, frameshift/stop evidence, and 13-lineage aggregation;
- `figures/species_shared_nonshared_loss_uniform`: the matching unified
  13-lineage shared/non-shared evidence figure (287 shared genes);
- `statistics/pseudogenized_positions_uniform`: the primary observed-target-
  locus/end-distance model for strict pseudogenized versus retained loci,
  including all 20,046 pseudogene target coordinates; deleted expected-locus
  midpoints are labelled sensitivity evidence only;
- `figures/pseudogenized_positions_uniform`: the matching 23-unit primary
  observed-target-locus position figure with deleted and uncertain absent;
- `figures/expression_copy_loss_uniform`: the stricter denominator-aware
  expression/copy sensitivity, not the article-comparable main figure;
- `tables/nlr_uniform` and `figures/nlr_uniform`: complete NLR repertoires and
  resolved non-shared reference-NLR loss percentages for all 23 units;
- `tables/functional_enrichment_uniform` and
  `figures/functional_enrichment_uniform`: GO, KEGG KO, and KEGG pathway
  enrichment for the unified positive-loss foreground, with a separate
  strict-pseudogenized-only sensitivity;
- `figures/functional_enrichment_categories_uniform`: category-level counts
  split into GO biological process, molecular function, cellular component,
  KEGG orthology, and KEGG pathway for pooled and all 13 lineage foregrounds;
- `figures/phylogeny_cafe`: the validated TimeTree-secondary-calibrated
  MCMCTree and CAFE5 Base Poisson changes on all non-root branches, with a
  two-tier period/epoch time scale;
- `figures/genome_evolution_overview`: the layout-only four-panel combination
  of the final Circos, OrthoFinder gene-composition, Ks-density, and
  phylogeny/CAFE figures;
- `figures/orthofinder_species_profiles`: the current frozen 17-taxon
  OrthoFinder core/species-specific orthogroup summary and a mutually
  exclusive per-taxon gene-composition panel. This dataset includes
  *Vitis vinifera* and *Rhododendron simsii*, not *R. vialii* or
  *R. delavayi*;
- `tables/assembly_annotation_qc`, `tables/callable_copy_opportunity`, and
  `tables/cafe5`: reviewer-facing QC tables and transparent CAFE5 exclusion
  ledgers;
- `phylogeny/topology_frozen`, `phylogeny/mcmctree_validation`, and
  `phylogeny/toolchain_audit`: accepted topology, posterior summaries, dated
  tree, secondary-calibration audit, and path-free executable/build bindings;
- `statistics/species_pgls_uniform`: the refitted exploratory PGLS and its explicit
  `BLOCKED_DENOMINATOR_AWARE_MODEL_REQUIRED` publication gate.

The article-comparable matrix contains 817,581 unit-gene rows. It classifies
633,957 rows as retained, 171,866 as decayed, 7,961 as deleted, and 3,797 as
outside the historical candidate scope. Its positive numerator is
`decayed + deleted`; the resolved denominator additionally includes retained.
There are 3,616 genes positive in all 23 units and 22,093 positive in at least
one unit. The 13-lineage tree layer has 34,637 genes with complete lineage
state information, while 910 genes with a not-called lineage remain outside
branch inference. Exact patterns comprise 3,280 single terminal-branch, 3,846
single internal-branch, and 5,620 repeated-independent-loss genes; 8,472 genes
have only partial/homeolog-specific lineage loss.

The primary unit-resolved scaffold does not collapse haplotypes or
subgenomes. Among 35,547 genes it classifies 5,250 as one unit-terminal event,
4,032 as one internal-node event, 11,510 as repeated independent events,
13,206 as no loss, and 1,549 as unplaced because at least one unit is
not-called. The resulting 56,602 event–gene memberships span all 39 scaffold
nodes and terminals.

The stricter sensitivity matrix also contains 817,581 rows. It classifies 57,027 rows
as deleted and 20,046 as assembly-sequence-supported pseudogenized, while
uncertain rows remain outside positive numerators and resolved denominators.
The 13-lineage shared-positive set contains 287 genes. The rebuilt NLR bundle
excludes four shared reference NLR genes and reports 210 non-shared reference
NLR genes, 2,487 resolved unit comparisons, and 1,003 positive calls.

The primary functional analysis keeps all 23 assembly units independent. Each
foreground contains every `decayed + deleted` gene in that unit and does not
require retained evidence anywhere else. Its unit-specific background contains
`retained + decayed + deleted` in the same unit, excluding `not_called_loss`.
The 23 foregrounds contain 179,827 unit–gene memberships and retain 6,420
GO/KEGG rows at BH-adjusted `q <= 0.05`; no unit has an all-zero functional
result.

The matching unit-resolved scaffold functional layer tests all 39 node and
terminal event sets against the 33,998 genes resolved in all 23 units. It
retains 905 GO/KEGG rows at BH-adjusted `q <= 0.05`, and nodes without a
significant term remain explicit zero rows.

The supplementary 13-lineage topology-aware analysis tests 31 foregrounds and
retains 1,826 GO/KEGG rows at BH-adjusted `q <= 0.05`. Foregrounds include
exact terminal branches, exact internal branches, repeated independent loss,
the all-23-unit shared set, and partial/homeolog-specific patterns. Complete
lineage loss requires every assembly unit assigned to that biological lineage
to be positive; mixed polyploid states are never promoted to ancestral loss.
These enrichment counts are descriptive and may overlap across biologically
nested categories, for example the all-23-unit set and the Actinidia-stem set.

The narrow pure species-specific sensitivity removes that overlap: each of its
1,167 genes is assigned to exactly one terminal lineage and is resolved
retained in every other lineage. It retains 54 species–term rows at
BH-adjusted `q <= 0.05`, separated into GO biological process, molecular
function, cellular component, KEGG orthology, and KEGG pathway. It is not used
as the primary assembly-unit result. These homology-based labels are
hypothesis-generating and are not experimental functional validation.

For the pooled non-shared foreground, the stricter sensitivity identifies 29 GO,
27 KEGG KO, and 30 KEGG pathway terms at BH-adjusted `q < 0.05`. The strict-
pseudogenized-only sensitivity identifies 18, 21, and 32 terms, respectively.
No GO or KEGG term is significant for the 287-gene shared foreground. These
are homology-based enrichment annotations, not direct experimental functional
validation. Older GO identifiers absent from the checksum-bound ontology are
reported as unresolved and are not silently remapped.
The complete all-tested-term table remains in the external data store; this
repository keeps the significant rows, foreground definitions/gene IDs, and
the checksum-bearing summary needed to audit that larger output.

Figure labels are deliberately shorter than internal analysis-unit names:
*A. zhejiangensis* A/B, *A. rufa*, and *A. macrosperma*. The display-only
qualifier `unphased` and the multiplication sign are omitted consistently.
This presentation rule does not change any sample, lineage, or denominator.

The current expression/reference-family-size figure pools shared and
non-shared genes. Its numerator is `decayed` only and its denominator is
`retained + decayed + deleted` within the matching unit and expression bin or
family-size class; `not_called_loss` is excluded. Expression is the declared
*C. scandens* leaf raw featureCounts sample divided into 14 rank bins.
Reference family size is the original CD-HIT 0.90 reference-protein cluster
size, with classes 1--7 retained because each contains more than 100 reference
genes; it is not a target-species CNV estimate.

The matching NLR analysis excludes 138 of 214 reference NLR
genes because they are positive in all 23 units. The remaining 76 non-shared
reference NLR genes provide 1,738 resolved unit comparisons and 254
`decayed + deleted` calls. Repertoire counts are retained for each assembly
unit, and no species aggregation is used. The 254 calls comprise 14
frameshift-supported, 14 in-frame-stop-supported, 13 combined-disruption, 16
no-qualifying-hit, 0 truncation/partial-alignment candidate, and 197
residual-sequence mechanism-unresolved calls.
The orthogonal structural classification assigns those same 254 calls to
CC-NBARC (28), CC-NBARC-LRR (141), NBARC (15), NBARC-LRR (51), TIR (3),
TIR-LRR (0), TIR-NBARC (4), and TIR-NBARC-LRR (12). Structural class is
inherited from the exact NLR-Annotator call for the corresponding reference
gene; it is not inferred anew from a fragmented target residual.
The current publication figure aligns all 6,034 NLR-Annotator calls across the
23 complete assembly-unit repertoires (185--329 per unit) with 3,428
shared-plus-non-shared unit-level loss calls, partitioned by the same
structural classes. Two target calls in the additional
TIR-CC-NBARC-LRR class are retained explicitly. Its heatmap uses the requested
descriptive target-repertoire comparison, `loss / (annotated + loss)`, rather
than `loss / annotated`; this quantity remains labelled a relative burden
because the two inputs are different gene universes. The accompanying
reference-opportunity tables retain the formal resolved loss percentages.
The branch panel contains 188 exact NLR loss-event placements (45 terminal and
143 internal); multi-assembly species are merged only after requiring every
assigned haplotype or subgenome to be `decayed` or `deleted`.

The mechanism-stratified spatial analysis includes all 179,827 positive
unit-gene rows, of which 139,994 have a target-assembly residual coordinate.
It includes ordinary `decayed` rows when an existing qualifying residual
alignment supplies a coordinate. Among the 150,547 unresolved-residual rows,
69,436 are local to the expected interval, 7,706 are candidate
same-chromosome displacements, 40,094 are candidate interchromosomal
displacements, and 33,311 remain unlocalized. Candidate displacement is not
proof of rearrangement because a best genome-wide hit can be paralogous.
Chromosome-count heterogeneity is evaluated separately for every evidence
group against within-unit chromosome-length opportunities after HY4A label
harmonization. The *C. scandens* reference GFF contains no TE/repeat features,
so no TE-association result is claimed.

Strict pseudogenized evidence remains separate from all article-method
numerators. Its 20,046 unit–gene rows comprise 11,689 frameshift-only, 3,266
in-frame-stop-only, and 5,091 combined calls. This conservative layer supports
mechanistic classification and observed-target spatial analysis but is never
added on top of an already positive article `decayed` or `deleted` call.
For the decayed-only mechanism panel, the exact historical-class intersection
is 11,559 frameshift-only, 3,258 in-frame-stop-only, and 5,071 combined calls.
The remainder of the uniform strict total consists of 157 historically
deleted calls and one call outside the historical candidate scope.

The unified expression/copy figure remains a separately labelled stricter
sensitivity using `deleted + strict pseudogenized` and excluding uncertain
rows. It is not mixed into the article-method main trend.

The earlier `species_loss`, `species_shared_nonshared_loss`, `spatial_loss`,
`uniform_loss_positions`, `expression_copy_loss`,
`expression_copy_loss_manuscript_method`, `nlr_nonshared`, and `species_pgls`
directories are retained only as superseded evidence; they are not the current
production statistics.
