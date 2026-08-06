# End-to-end analysis pipeline

This document describes the intended production workflow. A stage is accepted
only when its manifest, exact command, software version, input checksums, log,
and validation report are present. Large inputs and intermediate outputs live
under the external data store; only code and curated small results belong in
Git.

The authoritative parameter registry is `config/analysis_parameters.toml`.
Interpretive inclusion and aggregation rules are explained in
`docs/DECISION_RULES.md`; commands in this document must not silently override
that registry.

## Stage 0: freeze identities and scopes

`config/assemblies.tsv` is the source of truth. Each row uses a readable
machine identifier such as `act_eriantha_hap1_2026`, while the full Latin name,
individual, haplotype/subgenome, ploidy, accession, and assembly scope remain
explicit columns.

Before computation, verify:

- biological species and individual identity;
- whether a file is a full assembly, chromosome subset, or annotation subset;
- matched genome, GFF, CDS, and protein versions;
- sequence-name agreement across the files;
- publisher checksum and local SHA-256;
- which analysis units belong to one biological species.

HAP1/HAP2 and A-F units are technical assembly units, not independent species
replicates.

## Stage 1: acquire and verify public assets

Downloads are declared in `config/downloads.tsv` and stored under directories
that contain real taxon and assembly names.

```bash
python scripts/download/fetch_manifest.py \
  config/downloads.tsv \
  --data-root "$DATA_ROOT" \
  --report "$DATA_ROOT/checksums/download_report.json" \
  --connections 4 \
  --proxy "$MIHOMO_PROXY" \
  --proxy-domain ftp.ncbi.nlm.nih.gov \
  --segmented-proxy-domain ftp.ncbi.nlm.nih.gov \
  --proxy-domain ndownloader.figshare.com \
  --proxy-domain genome.kazusa.or.jp
```

The downloader checks expected bytes and publisher MD5/SHA-256 when supplied,
and always calculates a local SHA-256. A WAF page, login response, or truncated
file cannot pass the declared size/checksum gate. Unavailable files remain
disabled rather than being silently replaced by a different assembly.

The proxy is selected per domain. CNCB/GWH downloads remain direct. Provider
updates and node names are private machine state; the repository records only
the routing policy and verified file checksums.

Rooting and calibration-evaluation outgroup assets are frozen separately so that downloading an
outgroup cannot silently add it to a biological analysis cohort:

```bash
python scripts/download/fetch_manifest.py \
  config/phylogeny/public_outgroup_downloads.tsv \
  --data-root "$DATA_ROOT" \
  --report "$DATA_ROOT/checksums/public_outgroup_download_report.json" \
  --connections 4
```

This command verifies or restores only the already frozen manifest rows. Do
not add another taxon to that manifest: the approved membership and
cardinality rules are in `config/phylogeny/minimal_outgroup_design.tsv`.

This ledger contains the exact *Leea coccinea* Lco1464 v1.0 genome--GFF pair,
the exact *Catharanthus roseus* `GCA_024505715.1_ASM2450571v1`
genome--GFF--protein--CDS quartet, and the published *Saurauia tristyla*
Angiosperms353 tree-only archive. *Leea* contributes at most one
declared representative haplotype per biological-species analysis. The
*Saurauia* archive is never an OrthoFinder proteome/count-matrix, dating, CAFE,
or gene-loss input. Its RNA-seq rescue rows are prohibited and are not
downloaded. Download completion alone changes no inclusion state.

Migrate the frozen manuscript-era inputs separately. The old absolute paths
are read from the private legacy manifest, while Git stores only the public
sample-to-real-name mapping:

```bash
python scripts/migration/migrate_legacy_assets.py \
  --mapping config/legacy_analysis_units.tsv \
  --legacy-manifest /private/path/assembly_qc_manifest.tsv \
  --data-root "$DATA_ROOT" \
  --copy-max-mib 25 \
  --report "$DATA_ROOT/checksums/legacy_asset_migration.json" \
  --resolved-manifest "$DATA_ROOT/checksums/legacy_qc_manifest.tsv"
```

With the default 25 MiB threshold, the legacy protein files are copied while
the much larger genomes and GFF files are soft-linked. All 22 assembly units
receive descriptive IDs and SHA-256 records.

## Stage 2: reconcile full and chromosome scopes

Run `scripts/qc/basic_stats.py` for every genome FASTA and use
`scripts/qc/audit_paired_assembly_scope.py` where both full and chromosome
scopes exist. Report total length, sequence count, N50/L50, N content, GC,
anchored length, anchored fraction, unplaced length, and organelle records.

The paired-scope manifest is assembly-unit based. Its exact first column is
`assembly_unit_id`; every output table and the run metadata retain that same
identifier. The complete non-empty manifest defines the cohort size, so this
audit has no built-in 22-row or manuscript-era assumption. When an independently
frozen ledger provides an exact expected size, add the optional fail-closed
assertion:

```bash
python scripts/qc/audit_paired_assembly_scope.py \
  --manifest "$DATA_ROOT/manifests/paired_assembly_scope.tsv" \
  --output-dir "$DATA_ROOT/qc/paired_assembly_scope" \
  --expected-assembly-unit-count 22
```

Omit `--expected-assembly-unit-count` for a newly declared candidate cohort;
the script then derives the count from the manifest while still rejecting an
empty manifest, duplicate identifiers, and incomplete rows. Official unplaced
status requires a same-accession full deposit that reconciles exactly with its
assembly report. Extra records in a combined or name-only candidate remain
explicitly labelled candidate-only scope and are never promoted to official
unplaced scaffolds.

Annotation scope is assessed separately: count genes, transcripts, proteins,
and CDS; list GFF sequence IDs absent from the FASTA and FASTA sequences absent
from the GFF. Unmatched records are retained in an audit table, not discarded.

The reviewer QC table must distinguish publisher-reported metrics from locally
recomputed values and must identify the exact scope behind each number.

## Stage 3: uniform assembly and annotation QC

BUSCO production QC is restricted to *Actinidia*. For an *Actinidia* assembly
already assessed in the frozen legacy project, import the existing genome and
protein BUSCO tables after exact sample reconciliation and input binding. Reuse
requires the current genome to resolve to the same frozen asset and the current
protein FASTA to have the same SHA-256 checksum. Run BUSCO only for newly
downloaded or genuinely changed *Actinidia* assemblies, including candidates
later excluded from gene-loss analysis. Non-*Actinidia* outgroups do not require
BUSCO for this revision. All newly run rows use the same BUSCO release, lineage
dataset, and date-stamped database copy.

```bash
python scripts/qc/run_busco_batch.py --help
python scripts/qc/collect_busco.py --help
```

After the runs finish, use `scripts/qc/build_qc_publication_tables.py` to join
the exact basic-statistics and BUSCO outputs to an explicit metadata/decision
table. The builder requires exact identifier reconciliation and validates the
BUSCO counts, percentages, software/database signature, and local metric
arithmetic before atomically publishing path-free tables. See
`docs/QC_PUBLICATION_TABLES.md` for the input schema and output contract.

`run_busco_batch.py` detects compression from gzip magic bytes. Gzip genome or
protein FASTA is never passed directly to BUSCO 5.8. Instead, the wrapper
stream-decompresses it in one thread into
`<output-dir>/staged_inputs/<mode>/`, validates UTF-8 FASTA form, and commits a
content-addressed plain FASTA plus deterministic SHA-256 provenance with
atomic file replacement. Reuse rechecks both the compressed source and staged
plain FASTA checksums. Partial compressed-input runs restart only when
`run_input_bindings/` identifies the exact source and staged checksums;
completed-summary reuse also requires BUSCO's recorded input path to match the
content-addressed stage. An older, unbound run directory therefore fails
closed. Review it and use either a new output directory or an explicit
`--force` rerun. Plain FASTA keeps the original direct-input restart and skip
behavior. `--validate-only` performs staging and validation but does not start
BUSCO.

Every future BUSCO invocation launched by `run_busco_batch.py` includes
`--opt-out-run-stats` exactly once. The wrapper validates this privacy contract
both while constructing the command and again at the execution boundary; an
omitted or duplicated flag fails before the BUSCO executable starts. The exact
flag-bearing command is recorded in the `command` column of `batch_status.tsv`;
for an executed job it is also recorded in the first line of both per-sample
logs. For `validated` and `skipped_complete` rows, `command` is a planned
future command and no new BUSCO process ran; the existing per-sample logs are
the authoritative execution record for an older retained run. BUSCO's
`--offline` option only prevents data-file downloads and does **not** disable
its independent anonymous run-statistics request.

A read-only code audit of the installed BUSCO 5.8.2 release found that, without
the opt-out flag, its run-statistics JSON can contain the input byte size and
MD5 checksum, configuration choices, tool versions, timings, and summary
results. The inspected code sends that JSON separately from download handling.
No code path adding raw FASTA sequence content to that JSON was found. This is
a source-code observation for the inspected BUSCO version, not a packet capture
or a guarantee about other releases or the remote service. Runs already active
before this contract was added are not modified or restarted.

The project-wide worker limit is 15. Individual heavy jobs normally remain at
10 workers or fewer, while batch scheduling ensures that the sum of
BUSCO/Augustus/BLAST worker pools never exceeds 15. The QC decision
uses completeness, duplication, fragmentation, missing BUSCOs, assembly scope,
gene count, annotation compatibility, and contamination signals together; it
is not based on one BUSCO percentage alone.

## Stage 4: standardize annotations

For each accepted assembly, extract one reproducibly selected primary isoform
per gene from the matched chromosome-scope genome and GFF. Validate every
transcript before selection; a longer invalid model must not make a gene
disappear when a shorter valid isoform exists. Production extraction is
atomic, records all omitted models, and requires exact independent gffread
agreement for selected IDs, CDS sequences, and proteins. The independent
gffread command receives the chromosome FASTA and the generated selected-only
`sample.primary.gff3` (gene/mRNA/CDS), not the full publisher GFF3. Its CDS and
protein outputs must each have exact selected-ID closure. The full publisher
GFF3 remains the checksummed input to Python validation and selection.

A publisher `gene -> CDS` graph with no transcript rows is rejected by
default. For a release whose graph has been explicitly audited, use the narrow
`--gene-as-transcript` fallback. It preserves each publisher gene ID as one
self-transcript, rejects mixed/multi-parent/undeclared-parent graphs, audits
noncoding pseudogenes without treating them as invalid proteins, and runs the
same exact gffread sequence gate on a normalized top-level mRNA/CDS object.
See `docs/PRIMARY_ANNOTATION_STANDARDIZATION.md` for the complete contract.

```bash
python scripts/qc/extract_primary_annotation.py \
  --genome input.chromosome_scope.genome.fasta \
  --gff input.chromosome_scope.gff3 \
  --sample-id act_eriantha_hap1_2026 \
  --output-dir standardized/act_eriantha_hap1_2026/primary_annotation \
  --gffread /path/to/gffread \
  --require-gffread
```

Then require compatibility with the publisher protein set from the same frozen
release. This second gate requires exact protein ID-set equality. It permits
only terminal-stop normalization and individually audited publisher `X`
positions against canonical derived residues; all other differences reject the
run.

When the publisher FASTA first-token IDs are protein accessions rather than
the selected GFF3 transcript IDs, first run the explicit accession adapter.
It requires complete one-to-one GFF3 transcript/accession/protein closure,
validates the corresponding publisher header fields, subsets to the exact
selected-primary ID set, and proves character-for-character sequence
preservation. Source non-primary isoforms remain in its audit inventory; an
unmapped source record, missing selected record, ambiguity, or output extra
fails atomically.

For an audited publisher `gene -> CDS` graph, the same adapter must receive
the explicit `--gene-as-transcript` option used during primary extraction. It
preserves the gene ID as the self-transcript ID, requires exactly one declared
gene parent and one declared protein accession across all CDS parts, and
rejects mixed transcript graphs, pseudogene CDS, or ambiguity in selected and
non-selected models. When a normal transcript graph or a self-transcript graph
uses the transcript/CDS parent itself as an accession, declare that identity
with `--transcript-accession-source transcript_id` or
`--protein-accession-source cds_parent`; do not reuse attribute-name options to
imply two different semantic fields. Exact ARU, Fuchu, and ActinidiaBase v1
commands are frozen in `docs/PUBLISHER_PRIMARY_PROTEIN_REMAP.md`.

```bash
python scripts/qc/remap_publisher_primary_proteins.py \
  --selected-primary-proteins standardized/sample/primary_annotation/sample.protein.faa \
  --gff standardized/sample/chromosome_scope.gff3.gz \
  --publisher-proteins raw/sample/publisher.protein.faa.gz \
  --sample-id sample \
  --output-dir standardized/sample/publisher_primary_remap
```

See `docs/PUBLISHER_PRIMARY_PROTEIN_REMAP.md`. The remapped FASTA, not the
unfiltered accession-ID FASTA, is then the publisher input to the exact
sequence compatibility gate.

```bash
python scripts/qc/audit_published_protein_compatibility.py \
  --derived-proteins standardized/act_eriantha_hap1_2026/primary_annotation/act_eriantha_hap1_2026.protein.faa \
  --publisher-proteins raw/act_eriantha_hap1_2026/publisher.protein.faa.gz \
  --sample-id act_eriantha_hap1_2026 \
  --output-dir standardized/act_eriantha_hap1_2026/publisher_protein_compatibility
```

Do not send a compatibility failure downstream. First reconcile the exact
genome/GFF/protein release and any audited publisher-primary subset. See
`docs/PUBLISHED_PROTEIN_COMPATIBILITY.md`.

An explicitly configured canonical transcript tag takes precedence among valid
isoforms. Otherwise the tie-break is longest validated spliced CDS, longer
genomic span, then lexical transcript ID. Missing CDS phase and genes with CDS
but no valid coding transcript fail by default. See
`docs/PRIMARY_ANNOTATION_STANDARDIZATION.md` for the full policy, explicit
compatibility overrides, output schemas, and acceptance checklist.

Do not concatenate or relabel haplotypes before QC. For a combined six-haplotype
bundle, first validate the declared header rule and then split genome, GFF, CDS,
and protein with the same mapping. Every derived unit receives its own checksum
and a parent-bundle reference.

### Stage 4b: assign final chromosome homology labels

Keep publisher-scope `PubChr` labels until each 29-chromosome unit has a
complete global one-to-one maximum-nucleotide-similarity assignment to
Hongyang v4.0 HY4A. The production result must contain each `Chr01`--`Chr29`
label exactly once. Absolute coverage, reciprocal-best/separation measures,
HY4P agreement, and independent JCVI anchors are confidence diagnostics rather
than naming blockers. Rebuild the final FASTA/GFF from the frozen provisional
inputs after mapping; never edit them in place or copy one haplotype/subgenome
map to another. Preserve publisher sequence direction, and require exact
CDS/protein sequence closure. See `docs/CHROMOSOME_HOMOLOGY_RENUMBERING.md`.

Once the validated HY4A nucleotide matrix and the independent HY4P/JCVI
diagnostics exist, build the immutable similarity label map with
`scripts/qc/diagnose_chromosome_assignments.py`, then materialize the matched
FASTA/GFF bundle with `scripts/qc/relabel_chromosome_bundle_by_similarity.py`.
The production worker is label-only: every orientation action is `keep`.

The older strict multi-evidence audit can still be run with:

```bash
python scripts/qc/assign_chromosome_homology.py \
  --nucleotide-hy4a "$DATA_ROOT/homology/UNIT/nucleotide_hy4a.tsv" \
  --jcvi-hy4a "$DATA_ROOT/homology/UNIT/jcvi_hy4a.tsv" \
  --nucleotide-hy4p "$DATA_ROOT/homology/UNIT/nucleotide_hy4p.tsv" \
  --jcvi-hy4p "$DATA_ROOT/homology/UNIT/jcvi_hy4p.tsv" \
  --nucleotide-hy4a-provenance "$DATA_ROOT/homology/UNIT/nucleotide_hy4a.provenance.json" \
  --jcvi-hy4a-provenance "$DATA_ROOT/homology/UNIT/jcvi_hy4a.provenance.json" \
  --nucleotide-hy4p-provenance "$DATA_ROOT/homology/UNIT/nucleotide_hy4p.provenance.json" \
  --jcvi-hy4p-provenance "$DATA_ROOT/homology/UNIT/jcvi_hy4p.provenance.json" \
  --parameters config/analysis_parameters.toml \
  --target-asset-registry "$DATA_ROOT/homology/UNIT/target_assets.tsv" \
  --reference-asset-registry config/chromosome_coordinate_references.tsv \
  --reference-chromosome-map-registry config/chromosome_reference_maps.tsv \
  --assembly-unit-id UNIT \
  --target-scope-id UNIT.chromosome_scope_v1 \
  --trusted-repository-commit "$REVIEWED_GIT_COMMIT" \
  --output-dir "$DATA_ROOT/homology/UNIT/assignment_v1"
```

This diagnostic command performs no alignment. It validates exact 29-by-29
identifier and arithmetic closure, solves each matrix independently with a
deterministic global Hungarian assignment, and emits its strict map only on
`PASS_AUTO`. Each sidecar binds the matrix role, unit, scope, target
genome/GFF/protein, exact
HY4A or HY4P assets, frozen chromosome-ID map, generation parameters, and one
upstream report. Nucleotide generation is frozen to minimap2 `2.28-r1209` with
the exact argv
`minimap2 -x asm5 --secondary=no -c --cs=long {reference_fasta} {query_fasta}`;
the target/reference FASTA precedes the query FASTA, and accepted PAF rows must
carry `tp:A:P`, `de:f`, `cg:Z`, and `cs:Z`. Here `de` is minimap2's emitted
gap-compressed divergence estimate. The report must repeat those
identities, use accepted builder workflow version `1.0.0`, and pass the exact
nucleotide-PAF or JCVI-anchor validation set. The strict diagnostic requires a
1.5 top/second ratio and 0.10 normalized margin; assigned JCVI edges additionally
require at least 30 unique anchor pairs, 0.05 reciprocal gene coverage, and
0.05 assigned score. Valid conflicts and low-support evidence are atomically
retained as diagnostic audits; they do not override the author-approved
maximum-similarity production map. Diagnostic publication is protected by an
exclusive lock and no-replace rename, so an existing output is never overwritten. This
local validation command neither connects to a server nor accesses the network.
Its manifest records the reviewed Git commit plus captured policy and registry
hashes; the explicit trust boundary is defined in
`docs/CHROMOSOME_HOMOLOGY_RENUMBERING.md`.

## Stage 5: synteny and callable-coverage gate

The assembly inclusion percentage is based on the JCVI collinearity analysis
against the verified *Clematoclethra* reference. Prepare BED/protein inputs with
`prepare_jcvi_bed.py` and run the same JCVI parameters in both directions.

For each assembly unit, report at least:

- reference genes in anchors / all eligible reference genes;
- target genes in anchors / all eligible target genes;
- reference and target chromosome coverage;
- duplicated anchor rows and one-to-many behavior;
- chromosome-name and scope reconciliation.

One directional percentage is insufficient. An assembly failing the frozen
bidirectional/callable threshold is excluded from the primary gene-loss cohort
but remains in the assembly-QC table and may remain in the phylogeny if protein
QC and orthologue occupancy pass.

SynOrths provides the gene-level orthology evidence used by the loss caller.
An unchanged legacy assembly unit reuses its archived SynOrths table only after
the table passes identifier-column, FASTA, coordinate, denominator, duplicate,
and checksum binding against the exact frozen legacy inputs. A newly downloaded
or changed assembly unit is run from clean standardized inputs. Changing the
genome, annotation, chromosome scope, selected primary proteins, or coordinates
invalidates legacy reuse and requires a fresh run. JCVI and SynOrths answer
related but different questions; neither output is substituted for the other.

## Stage 6: call candidate and positive gene losses

The reference-gene universe is frozen first. For each accepted assembly unit:

1. Normalize raw SynOrths output with `geneloss normalize-synorth`.
2. Identify unmatched reference genes bracketed by local anchors with
   `geneloss call-candidates`.
3. Extract every candidate query without silently dropping missing IDs.
4. Search the target genome with tBLASTX using one recorded schema and threshold
   set.
5. Classify with `geneloss classify-tblastx`.

Production classification uses the callable synteny-aware mode:

- `positive_deleted`: bilateral same-target-chromosome SynOrths anchors define
  a callable local interval and no qualifying translated hit is detected there;
- `uncertain_local_genomic_sequence_detected`: a qualifying local translated
  hit exists, but disruptive-mutation evidence has not been established;
- `uncertain`: the local interval is not callable or required evidence is
  missing.

The manuscript-era genome-wide-hit rule, including its historical
`decayed=pseudogenized` interpretation, is retained only as a labelled legacy
sensitivity. Neither uncertain state is a biological loss class, and both are
excluded from positive-loss rates. The output records a decision reason and
the exact best hit for every candidate.

## Stage 7: aggregate assembly units into species

First build complete per-unit call tables; absence from the positive lists is
called `not_called_loss`, not automatically `retained`.

For species with multiple units, apply the predeclared aggregation rule. For
example, an *A. eriantha* species-level loss requires concordant positive and
callable evidence in both HAP1 and HAP2. A one-haplotype call is reported as
partial/haplotype-specific or uncertain.

After the accepted biological-species cohort is frozen, regenerate:

- the union of all positive losses;
- losses supported in every biological species (shared losses);
- non-shared losses after removing that intersection;
- prevalence across species and across assembly units;
- lineage-restricted and partial-haplotype categories.

No historical shared-loss number or denominator is reused. Shared losses are
reported separately; comparative analyses requested by the reviewer use the
non-shared set.

Species-level aggregation is performed from a complete, callable-aware unit
matrix with `scripts/gene_loss/aggregate_species_loss.py`. The script writes an
auditable four-state species matrix (`positive_complete`, `positive_partial`,
`not_positive`, and `uncertain`) before calculating the shared intersection.

## Stage 8: downstream analyses

All downstream analyses consume the same accepted positive-call master table
and the same species aggregation ledger.

### Expression

Map expression identifiers to the exact annotation version and retain the
arithmetic mean TPM across the four declared *C. scandens* tissues. Divide the
reference genes into 14 rank-first expression bins and calculate a separate
rate for every genome unit and bin. Shared and non-shared calls are included
together. The numerator is article-method `decayed`, and the resolved
denominator is `retained + decayed + deleted`; `not_called_loss` is excluded.
Fit the reported overall regression to the complete unit-by-bin plot table and
retain the plotted data and fit statistics.

### Copy number

In the production analysis, copy number is the size of the *C. scandens*
reference-protein CD-HIT 0.90 similarity cluster. It is therefore a reference
gene-family-size measure, not a target-species CNV call. Keep classes 1--7,
each supported by more than 100 reference genes. For every genome unit and
class, use the same shared-plus-non-shared article-method `decayed` numerator
and resolved denominator used for expression, and retain the fitted overall
relationship.

### Chromosome position

Restrict the primary position analysis to article-method `decayed` loci with
an observed residual-sequence coordinate in the corresponding target genome.
Include shared and non-shared loci together. Exclude deleted calls because
they have no observed target locus, and exclude spatially unlocalized decayed
calls. Harmonize the 29 homologous *Actinidia* chromosome groups to HY4A
`Chr01`--`Chr29` while retaining each genome unit as an independent model row.

Test between-chromosome heterogeneity with a negative-binomial model containing
genome-unit and chromosome effects and the log annotated-gene count as an
offset. Test within-chromosome position after assigning each target gene and
decayed locus to one of five equal zones from the nearest chromosome end to
the centre. The zone model contains genome-unit, chromosome, and zone effects,
again with the matching annotated-gene count as an offset. Treat frameshift,
in-frame-stop, and combined-disruption subsets as mechanistic sensitivities.
Use chromosome-length expectations only for the separate descriptive
placement summaries. Do not infer centromere or transposable-element effects
without independent annotations. See
`docs/DECAYED_CHROMOSOME_POSITION_ANALYSIS.md` for the pooling and
standardization rules.

### NLR repertoire

Run the same NLR annotation workflow on the complete released scope of each
assembly. For every unit report total NLR loci, positive reference-NLR loss
calls, the resolved NLR denominator, and percentage. A chromosome-only total
must be labelled as such and cannot be mixed with a full-assembly total.
The revised primary comparison first removes reference genes in the shared
positive-complete set. Its reference-NLR universe is the unique reference-CDS
sequence IDs called by NLR-Annotator. Each non-shared reference NLR classified
as `retained`, `deleted`, or strict `pseudogenized` contributes one resolved
unit-level denominator opportunity; `deleted` plus strict `pseudogenized`
enter the numerator. Both uncertain states are excluded from numerator and
denominator. Every target repertoire is obtained from its complete declared
29-chromosome analysis unit, while the reference run uses one CDS per
reference gene.

### GO and KEGG enrichment

Use the checksum-bound eggNOG-mapper annotation generated with the
Viridiplantae taxonomic scope. Analyze GO biological process, molecular
function, and cellular component, KEGG orthology, and KEGG pathway membership
as separate annotation systems.

The primary foregrounds are terminal complete-loss events on the accepted
13-lineage topology. Unit-level positive states are article-method `decayed`
or `deleted`. For a multi-unit species, complete loss requires every assigned
haplotype or subgenome to be positive; mixed states remain partial or
homeolog-specific. Remove genes assigned to an ancestral event on the focal
root-to-tip path from that lineage's risk set.

For each lineage and annotation system, fit a null logistic model containing
z-standardized log2(four-tissue mean TPM + 0.1), its squared term,
z-standardized log2(reference CD-HIT 0.90 family size), and its squared term.
Test each functional term with a one-sided efficient score test, requiring at
least five background genes and two terminally lost genes. Apply
Benjamini-Hochberg correction within lineage and annotation system, then refit
significant positive terms to estimate adjusted odds ratios and 95% confidence
intervals. The validated analysis contains 33,998 resolved reference genes,
33,974 genes with complete covariates, 19,192 terminal-event memberships,
32,591 tested lineage-term rows, and 3,646 significant rows.

Representative figure terms are selected deterministically from significant,
converged full-model fits. Rank terms by the number of significant lineages,
the smallest adjusted q value, the mean adjusted log2 odds ratio, and the term
identifier. Reduce GO ancestor-descendant redundancy when significant-lineage
sets overlap strongly. Retain the complete tested table outside this display
subset. Unit-level hypergeometric and topology-scaffold summaries remain
complementary descriptive analyses. See
`docs/COVARIATE_ADJUSTED_FUNCTIONAL_ENRICHMENT.md` for the exact selection
rules and category limits.

## Stage 9: figures and typography

Each final figure is generated from a committed small plot-data table and a
versioned script. Species labels use an italic Latin binomial and an upright
haplotype/subgenome suffix, for example *A. deliciosa* A. Figures must expose
real species/unit names; internal stage codes are not used as reader-facing
labels. Reader-facing labels are deliberately concise: *A. × zhejiangensis*
A/B is displayed as *A. zhejiangensis* A/B, while *A. rufa* and
*A. macrosperma* carry no technical phase suffix. These aliases change
presentation only; internal analysis-unit IDs and species aggregation remain
unchanged.

Required result panels include assembly/annotation QC, bidirectional JCVI
coverage, shared versus non-shared loss, expression, copy number, chromosome
position, NLR total/loss percentage, GO/KEGG enrichment, and the
phylogenomic/dating results.

## Stage 10: phylogeny and dating

Follow `docs/PHYLOGENY_POLICY.md` and
`docs/OUTGROUPS_AND_CALIBRATIONS.md`. The production chain is:

```text
validated representative proteomes
  -> OrthoFinder 3 orthogroups/HOGs and gene families
  -> per-locus MAFFT alignment, exact CDS back-translation, and codon trimming
  -> IQ-TREE 2 locus trees plus partitioned concatenation
  -> ASTRAL-Pro multi-copy coalescent tree
  -> frozen biological-species topology
  -> TimeTree-secondary-calibrated MCMCTree dating under the author-approved scope
  -> CAFE 5 on the matching dated tree and representative-complement counts
```

Build a separate assembly-unit diagnostic tree containing visible HAP,
haplotype, and subgenome suffixes. That tree is not dated and cannot feed CAFE
or PGLS. The primary tree uses one representative complement per biological
species according to `config/phylogeny/representation_policy.tsv`; alternate
polyploid complements are full rerun sensitivities rather than additional
species tips.

Legacy *Rhododendron*, *Coffea*, and *Vitis* assets are soft-linked and
checksum-audited, but selected for the primary dated tree only after
provenance, BUSCO, orthologue occupancy, and calibration-node checks. No fossil
is currently active. A fossil can be activated only when the exact required
bracketing taxa are present in the dated alignment/topology and have passed
asset QC. PGLS then runs only after the exact dated biological-species tree and
callable-aware shared/non-shared species aggregation both pass. Its exploratory
response is the lineage-specific/non-shared positive-loss proportion after
excluding shared positives, partial/uncertain calls, and non-callable genes
from the denominator. The ordinary Gaussian PGLS publication gate remains
blocked until a denominator-aware phylogenetic count model passes. The primary
predictor is `log2_ploidy`. Required sensitivities exclude *A. rufa*, replace
its accepted assembly and rebuild the matching loss/tree inputs, and prune
each biological species from the time tree in turn before refitting.

The exact genome-derived *Leea* and *Catharanthus* candidates follow the same
primary-isoform, compatibility, BUSCO, and occupancy gates as other proteome
tips. Use at most one *Leea* haplotype for the single *L. coccinea* species row,
and only after the predeclared root-stability diagnostic. *Leea* lies outside
Vitoideae, so it cannot activate the *Indovitis* crown-Vitoideae minimum with a
single *Vitis* tip. *Catharanthus* is a one-for-one *Coffea* rooting swap and
cannot bracket an internal Rubiaceae fossil. The
publisher Lco1464 assembly has contigs but no chromosome assignment, so split
its complete h1/h2 record sets by the frozen exact inventory and do not remove
records by length; see `docs/LEEA_LCO1464_SCOPE.md`. The ITS-only *Saurauia*
archive is external context: it is excluded from the OrthoFinder proteome/count
matrix, dating, CAFE family counts, and all gene-loss analyses, and it cannot
bracket a revised dated-tree node.

Dating has three separate manifests. Revised fossil-calibrated MCMCTree remains
blocked until at least one exact bracket passes. The author-approved revised
dated tree has instead completed under transparent TimeTree secondary calibration, with
complete query/version/count/checksum rows in
`config/phylogeny/secondary_timetree_constraints.tsv`. The two posterior chains
and ultrametric tree passed their declared ESS/R-hat and checksum gates. The
five recovered
legacy fixed points remain only in the exact `ape::chronos` reproduction ledger;
that reproduction is currently blocked by the audited 9-tip input versus
17-tip archived-output mismatch. See `config/phylogeny/dating_designs.tsv`.

CAFE5 then used that matching dated tree and one representative-complement
count per lineage. The accepted Base Poisson run analyzed 15,066 families,
estimated lambda as 0.0085698614157905, and reported 1,539 families at
`p < 0.05`. Gamma3 is `UNAVAILABLE_INITIALIZATION_FAILURE`; it is not a
supported production model and no further post hoc family deletion or retry is
allowed.

The sum of simultaneous scientific workers across all locus jobs and tools is
at most 15. Memory and system load are checked before additional jobs start.

## Stage 11: release gate

Run unit tests, manifest/schema validation, checksum reconciliation, privacy
scans, figure rendering checks, and workbook formula checks single-threaded.
Inspect staged Git files for manuscripts, reviewer material, credentials,
absolute server paths, raw genomes, and bulk outputs. Create one final private
local commit after all gates pass; do not push it.
