# Repository guide

This guide describes the repository organization and the data flow required to
reproduce the analyses.

## Two sibling roots

The Git repository contains reusable code, configuration, documentation,
tests, curated small tables, and final figures. The external data store contains
downloaded genomes, annotations, raw tool output, temporary files, logs, and
other large or private material.

```text
Actinidia-gene-loss/       version-controlled project
external-data/             user-defined large-data and run store
```

Deleting a derived run below the data store must never delete a source download
or a legacy input. Large legacy genome/GFF assets are soft-linked into the data
store; small protein assets are copied. Every migrated asset is checksum-audited.

## Repository directories

### `config/`

Human-reviewed, machine-readable declarations live here.

- `assemblies.tsv` declares public candidate bundles, identities, accessions,
  scopes, and provisional inclusion states.
- `candidate_assembly_catalog.tsv` is the broader, source-backed discovery
  inventory. It includes DNA-only, sensitivity, duplicate, unavailable, and
  excluded releases and never adds a row to an analysis cohort automatically.
- `downloads.tsv` declares exact asset URLs, destinations, sizes, and publisher
  checksums.
- `legacy_analysis_units.tsv` maps every manuscript-era label to a descriptive
  assembly-unit ID and biological species.
- `analysis_parameters.toml` freezes the project-wide worker ceiling, QC and
  JCVI parameters, loss definitions, spatial rules, and release policy.
- `cohorts.tsv` describes the intended QC, gene-loss, and phylogeny cohorts.
- `phylogeny/` contains taxon, outgroup-asset, and calibration registries.
  `analysis_designs.tsv` separates the primary biological-species,
  assembly-unit diagnostic, representative-swap, and manuscript-reproduction
  analyses. `representation_policy.tsv` freezes polyploid/haplotype selection
  and sensitivity rules before topology inference.
  `public_outgroup_downloads.tsv` freezes exact public outgroup assets and
  checksums; `taxa.tsv` and `fossil_bracketing_taxa.tsv` separately declare
  analysis roles and calibration-bracketing requirements. The current exact
  assets are the *Leea coccinea* Lco1464 genome--GFF pair, the *Catharanthus
  roseus* `GCA_024505715.1` genome--GFF--protein--CDS quartet, and the
  transcriptome-derived *Saurauia tristyla* Angiosperms353 tree-only archive.
  `leea_lco1464_records.tsv` freezes the 140 source contigs and their GFF
  membership, while `leea_lco1464_haplotype_scope.tsv` summarizes the two
  same-individual complete haplotype scopes without inventing chromosome
  labels; see `docs/LEEA_LCO1464_SCOPE.md`.
- `project.env.example` documents local path/proxy variables without containing
  credentials.

Configuration states intent. A generated decision ledger records which
candidates actually pass after QC; code must not silently edit configuration to
make a failed run pass.

### `src/geneloss_repro/`

Importable Python library code lives here. It contains the reusable annotation,
orthology, loss classification, statistics, spatial analysis, taxonomy-label,
plotting, and figure-bundle functions. Scripts should call this library rather
than duplicate biological logic.

### `scripts/`

Command-line entry points are grouped by analysis stage.

- `download/`: verified, resumable public-data acquisition with per-domain
  routing.
- `migration/`: legacy asset migration and path-free legacy QC import.
- `qc/`: chromosome-scope materialization, fail-closed primary-annotation
  standardization, basic assembly/annotation statistics, BUSCO orchestration,
  paired-scope audits, SynOrths preparation/execution, JCVI execution,
  nucleotide-plus-JCVI chromosome-homology assignment, and manifest
  resolution. The assignment command consumes four precomputed score matrices,
  validates exact 29-by-29 closure, and cannot publish final `Chr` labels from
  ambiguous or conflicting evidence. Paired-scope audit manifests and outputs use
  `assembly_unit_id`; their cohort size is derived from the manifest unless an
  explicit exact-count assertion is supplied.
- `gene_loss/`: complete unit matrices, callable-aware species aggregation, and
  shared/non-shared summaries.
- `spatial/`: chromosome-bin, chromosome-end, and optional independently
  supported centromere-distance analysis.
- `downstream/`: expression and copy-number relationship tables.
- `function/`: GO and KEGG foreground preparation, covariate-adjusted tests,
  and supplementary enrichment summaries.
- `statistics/`: deterministic exact ploidy-group and related comparisons.
- `nlr/`: bounded NLR-Annotator execution, validation, repertoire/loss
  reconciliation, and summaries.
- `figures/`: publication renderers. Each renderer writes PNG, PDF, exact plot
  data, caption, validation, and a checksum manifest.
- `phylogeny/`: alignment trimming, back-translation, concatenation, and taxon
  closure checks for the diagnostic and species trees.

An internal identifier may be short for tool compatibility, but every
reader-facing table and figure must use the real biological species and
haplotype/subgenome label.

### `tests/`

Small synthetic fixtures exercise schema validation, edge cases, denominators,
uncertainty, privacy, typography, atomic output, and failure behavior. Tests do
not require or contain manuscript data. Rendering tests use the optional
Matplotlib dependency; non-rendering validation remains testable without it.

### `results/`

Only curated, publication-safe small outputs belong here. A canonical figure
directory contains the figure in PNG and PDF, the exact plot-data TSV, an
English caption, validation, and a checksum manifest. Raw BUSCO, JCVI,
SynOrths, BLAST, and NLR output stays in the external data store.

### `docs/`

- `PIPELINE.md` gives the stage-by-stage workflow.
- `LOSS_CLASSIFICATION.md` defines the primary loss states, denominators, and
  downstream analysis sets.
- `COVARIATE_ADJUSTED_FUNCTIONAL_ENRICHMENT.md` defines the 13-lineage
  species-specific foregrounds, risk-set backgrounds, and covariates.
- `RESULTS_PROVENANCE.md` links reported results to curated evidence bundles.
- `DECISION_RULES.md` explains inclusion, species aggregation, and spatial
  interpretation.
- `ASSEMBLY_SELECTION.md` records candidate-specific evidence and decisions.
- `CANDIDATE_ASSEMBLY_CATALOG.md` explains the public assembly inventory,
  asset-availability states, biological independence, and deduplication rules.
- `DATA_LAYOUT.md` explains storage and proxy boundaries.
- `OUTGROUPS_AND_CALIBRATIONS.md` defines the phylogeny-only outgroup and fossil
  policy.
- `PHYLOGENY_POLICY.md` separates the biological-species tree from the
  assembly-unit diagnostic tree.

## End-to-end data flow

1. **Declare candidates.** Record identities, releases, scopes, URLs, checksums,
   and biological groupings before computation.
2. **Download and verify.** Download only enabled assets. A file advances only
   after declared size/checksum and calculated local SHA-256 reconciliation.
3. **Migrate reproduction inputs.** Copy small legacy assets, link large ones,
   and retain a checksum inventory. These inputs reproduce the submitted
   analysis but do not automatically enter the revised cohort.
4. **Resolve one runtime manifest.** Join declarations, enabled downloads, the
   verified report, and bytes on disk. Downstream code consumes this generated
   manifest rather than ad hoc paths.
5. **Audit scope and quality.** Calculate assembly/annotation statistics and run
   genome plus protein BUSCO for every candidate. Retain unmatched chromosome,
   unplaced, organelle, and annotation-only records in explicit audit tables.
6. **Run the JCVI gate.** Compare every candidate to the same verified
   *Clematoclethra scandens* reference with frozen parameters. Report reference
   and query coverage together. Apply the submitted 50% reference-side gate
   consistently; do not exempt a taxon because it was in the submitted cohort.
7. **Freeze the accepted ledger.** Select the best matched *A. deliciosa*
   bundle, evaluate both *A. eriantha* haplotypes as one individual, and retain
   or exclude *A. rufa* from the primary loss cohort according to the same QC,
   JCVI, SynOrths, and callable-scope rules. All candidates remain in QC tables.
8. **Standardize accepted units.** Split validated multi-haplotype bundles with
   one declared rule across genome, GFF, CDS, and protein. Select primary
   isoforms only after CDS validation, independently reconcile every selected
   CDS/protein sequence by running gffread on the generated selected-only
   gene/mRNA/CDS GFF3, and retain parent-bundle checksums. See
   `PRIMARY_ANNOTATION_STANDARDIZATION.md`.
9. **Call gene loss.** Build candidate losses from syntenic anchors, execute the
   recorded tBLASTX search, and classify each comparison as `retained`,
   `decayed`, `deleted`, or `not_called_loss` under the uniform rules in
   `LOSS_CLASSIFICATION.md`. A qualifying residual homologous sequence supports
   `decayed` but does not by itself establish coding disruption. Frameshift and
   premature-stop evidence identifies a strict pseudogenized subset and never
   rewrites the primary `decayed` or `deleted` class.
10. **Summarize genomes and biological lineages.** Preserve all 23 genomes for
    per-genome comparisons, then combine same-species haplotypes or subgenomes
    for branch-based analyses while retaining *A. zhejiangensis* A and B
    separately. Complete loss requires every constituent genome to be
    `decayed` or `deleted`; mixed states remain partial or homeolog-specific.
11. **Run downstream analyses.** Preserve the declared numerator and denominator
    for every comparison. The primary GO/KEGG analysis uses 13 biological
    lineages. For each focal lineage, genes assigned to tree-node losses on its
    root-to-lineage path are removed. Complete species-specific losses form the
    foreground; the other annotation- and covariate-complete genes in the
    lineage-specific risk set form the comparison background. Logistic score
    tests adjust for four-tissue mean expression and *C. scandens* gene copy
    number (CD-HIT 90% cluster size). Per-genome hypergeometric and strict-
    evidence enrichments remain supplementary analyses.
12. **Render figures.** Use real taxon metadata. Latin binomials are italic;
    haplotype, subgenome, accession, and scope suffixes are upright. Use the
    concise aliases *A. zhejiangensis* A/B, *A. rufa*, and
    *A. macrosperma* without changing internal IDs. No internal stage code is
    reader-facing.
13. **Build phylogenies separately.** Run the declared OrthoFinder 3 ->
    IQ-TREE 2/ASTRAL-Pro -> MCMCTree -> CAFE 5 chain for one
    biological representative per species. Soft-linked *Rhododendron*,
    *Coffea*, and *Vitis* assets are evaluated for the species tree and dating
    but never enter the gene-loss denominator. With no active fossil bracket,
    the completed revised tree uses the four author-approved, checksum-bound
    TimeTree 5 secondary intervals and is not called fossil-calibrated. Keep
    all assembly units in a separate undated diagnostic tree. The archived
    ordinary PGLS is exploratory and is not used for primary inference because
    it does not model unequal resolved denominators.
    Use one declared *Leea coccinea* haplotype for its one species row. The
    transcriptome-only *Saurauia tristyla* Angiosperms353 asset is tree-only and
    is excluded from the OrthoFinder proteome/count matrix, CAFE, and gene-loss
    analyses. The validated CAFE output is Base Poisson; Gamma3 is explicitly
    unavailable after initialization failure and is not retried by deleting
    further families.
14. **Release checks.** Re-run tests single-threaded, verify generated figures
    and workbook formulas visually, scan tracked files for private paths,
    credentials, manuscript files, and bulk data, and inspect the complete
    staged diff before publishing a release.

## Meaning of common states

- `candidate`: downloaded or migrated, but not yet accepted for biological
  inference.
- `accepted`: passed the frozen identity, scope, QC, JCVI, SynOrths, and
  callable-evidence gates for its declared purpose.
- `excluded`: retained in the audit trail but absent from the primary analysis.
- `sensitivity`: intentionally analyzed separately; its denominator must not be
  merged into the primary cohort.
- `uncertain`: evidence is insufficient for a positive or confident-negative
  biological statement. It is not a third loss class.

The word `terminal` is reserved for a phylogenetic tree leaf. Elsewhere use
`assembly unit` or `biological species`.
