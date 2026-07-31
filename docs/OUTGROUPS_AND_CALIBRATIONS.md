# Outgroups and fossil-calibrated dating

## Scope separation

Outgroups used to root or date a species tree are not gene-loss comparison
units. The gene-loss denominator contains only the declared *Actinidia*
cohort and its verified *Clematoclethra* reference. *Rhododendron*, *Coffea*,
and *Vitis* must never be added to shared/non-shared loss matrices merely
because their files are available to the phylogeny workflow.

The project retains the legacy outgroup assets by relative soft link so that
the submitted analysis can be reconstructed. Link availability is not proof
of taxonomic identity, accession, annotation compatibility, or suitability
for the new primary tree. Those properties are audited independently.

## Approved minimal-outgroup policy

No additional outgroup taxon is downloaded after this policy freeze. The
machine-readable contract is
`config/phylogeny/minimal_outgroup_design.tsv`. It distinguishes the revised
primary tree, one-for-one rooting sensitivities, exact legacy reproduction,
and global exclusions. Already acquired files are retained for provenance;
retention does not guarantee inclusion.

The revised primary candidate set adds only the following to the focal
*Actinidia*--*Clematoclethra* cohort:

1. exactly one existing *Rhododendron* representative, chosen by frozen QC and
   orthogroup occupancy rather than by a preferred topology;
2. the existing *Coffea arabica* bundle, if its identity, annotation, BUSCO,
   and occupancy gates pass;
3. the existing *Vitis vinifera* bundle, after its accession/release and
   annotation provenance close;
4. zero or one *Leea coccinea* tip, added only if one Lco1464 haplotype passes
   and a predeclared *Vitis*-only versus *Vitis*+*Leea* diagnostic demonstrates
   useful root stability.

*Catharanthus roseus* is not co-included by default. Its already acquired
matched bundle is a one-for-one replacement for *Coffea* in a labelled rooting
sensitivity. This tests whether the distant asterid representative affects the
topology without growing an unconstrained outgroup collection.

This policy minimizes missing orthologues, long compute cascades, and
opportunities to select taxa after seeing a desired tree. More importantly,
it prevents a downloaded taxon from being mistaken for a valid fossil bracket.
Every optional tip must have a declared topology or calibration role before it
can enter an alignment.

## Frozen public outgroup assets

Exact downloadable assets added for rooting and calibration evaluation are declared in
`config/phylogeny/public_outgroup_downloads.tsv`. This file uses the same
fail-closed size/checksum schema as the main download manifest, but membership
in it is only an acquisition decision. A taxon is not accepted into a tree,
OrthoFinder input, or family-count matrix until its declared purpose-specific
QC gates pass.

- *Leea coccinea* uses the exact Lco1464 v1.0 genome and paired gene annotation
  from Zenodo record 13362874. Primary CDS/proteins are derived and audited
  from that genome--GFF pair. The production biological-species tree and
  OrthoFinder species matrix may contain one declared Lco1464 representative
  haplotype only; phased homologues must not be counted as additional species.
  Lco1464 is a phased contig assembly, not a publisher-defined chromosome
  assembly. Its exact 63-record h1 and 77-record h2 scopes and the prohibition
  on length-based chromosome inference are frozen in
  `docs/LEEA_LCO1464_SCOPE.md` and
  `config/phylogeny/leea_lco1464_records.tsv`.
  *Leea* is an optional rooting diagnostic, not a completed *Indovitis*
  calibration bracket. The primary *Indovitis* description places the fossil
  within Vitoideae and favors two internal alternatives
  ([Manchester et al. 2013](https://doi.org/10.3732/ajb.1300008)). *Leea* lies
  outside Vitoideae; therefore *Leea* plus the single current Vitoideae tip,
  *Vitis*, cannot define crown Vitoideae. The fossil remains disabled because a
  second crown-spanning Vitoideae lineage is absent.
- *Catharanthus roseus* uses the matched NCBI
  `GCA_024505715.1_ASM2450571v1` genome, GFF3, protein, and CDS quartet, plus
  the official assembly report and assembly statistics. The assembly report
  defines chromosome, unlocalized, and unplaced sequence scope without name
  heuristics. The
  published protein/CDS files are retained as independent compatibility
  evidence for sequences re-extracted from the exact genome--GFF pair.
  *Catharanthus* can replace *Coffea* in a Gentianales rooting sensitivity, but
  it cannot define a node internal to Rubiaceae. The prior ledger entry that
  treated it as an internal Rubiaceae bracket has been corrected.
- *Saurauia tristyla* occurs in the published trees and concatenated alignment
  from Zenodo record 13835677, but the exact archive audit found zero matching
  records among the Angiosperms353 assembled-sequence files. Its concatenated
  record has only 685 non-missing sites, restricted to the appended ITS region.
  The archive is therefore external topology/age context only: it cannot be
  grafted into the new nuclear alignment and cannot enter OrthoFinder,
  ASTRAL-Pro, CAFE, or any gene-loss numerator or denominator. See
  `docs/SAURAUIA_TREE_ONLY_AUDIT.md`. The two registered NCBI paired-end
  RNA-seq runs are evidence records only and are explicitly not downloaded
  under the approved minimal-outgroup policy; see
  `config/phylogeny/saurauia_nuclear_rescue_candidates.tsv`.

`config/phylogeny/taxa.tsv` declares the allowed analysis roles, while
`config/phylogeny/fossil_bracketing_taxa.tsv` records which calibration each
taxon could bracket. Neither registry activates a fossil automatically.

## Planned tree sets

### Biological-species backbone

This is the primary topology tree. It contains one declared representative
haploid-complement tip per biological species, the verified *Clematoclethra*
lineage, and the smallest justified outgroup set. Haplotypes and subgenomes
from one individual are not independent species tips.

### Assembly-unit diagnostic tree

This tree retains HAP1/HAP2 and polyploid A-F units with explicit suffixes. It
tests phasing, contamination, orthologue occupancy, and topology stability.
It is not used to inflate biological sample size or to define a fossil node.

### Dated sensitivity trees

Alternative outgroup selections and calibration bounds are run separately.
The legacy *Rhododendron*, *Coffea arabica*, and *Vitis vinifera* assets may be
used after provenance and BUSCO checks. No new diploid *Coffea* is acquired.
The already acquired *Catharanthus* bundle supplies the declared one-for-one
rooting swap, and *Leea* supplies the conditional Vitales-root diagnostic.

## Calibration gate

No historical point age is reused without evidence. Every calibration record
must include:

1. the two descendant clades that define the calibrated node;
2. fossil name and specimen or a clearly identified secondary source;
3. minimum and maximum or distributional bound;
4. the biological justification for placement;
5. whether the bound is hard or soft;
6. the sensitivity runs that omit or alter the calibration.

The candidate fossil ledger is stored in
`config/phylogeny/calibrations.tsv`. Candidate rows remain disabled until the
required bracketing taxa are present in the exact alignment and topology used
for dating. A fossil's oldest stratigraphic bound is not treated as a node
maximum: each primary fossil supplies a conservative hard minimum, while a
broad root soft maximum is tested separately.

There are currently **no active fossil calibrations**. A fossil row can receive
`status=active` only after its `activation_gate` is
`passed_exact_bracketing_taxa_and_asset_qc`. A prose requirement alone is not a pass: the run
manifest must name the exact sampled descendant and sister-lineage taxon IDs,
their checksums, and the node they bracket. Removing a bracketing taxon
automatically disables that calibration for the affected sensitivity.

The previously discussed 55 Ma constraint on the *R. delavayi*-*R. simsii*
pair is not part of the recovered exact legacy evidence and is not reused. The
surviving script instead records 41.86 Ma for that pair; it remains legacy-only
under the provenance block below. The existing core-*Rhododendron* tips do not represent crown
*Rhododendron*; a Therorhodion tip such as *R. camtschaticum* is required, and
*R. newburyanum* remains sensitivity-only because its genus-level placement is
uncertain. One *Coffea* tip cannot bracket an internal Rubiaceae node and one
*Vitis* tip cannot bracket crown Vitoideae.

The minimal set does not currently close any fossil bracket. *Parasaurauia*
remains disabled because no compatible nuclear *Saurauia* asset is admitted;
*Paleoenkianthus* and *Rhododendron newburyanum* remain disabled because the
required basal descendants are not acquired. *Indovitis chitaleyae* is a
candidate crown-Vitoideae minimum, not a crown-*Vitis* shortcut; *Leea* cannot
replace the missing second Vitoideae descendant. The Rubiaceae row now has no
numerical bound because the recorded family-wide review does not by itself
freeze an internal fossil placement, and *Catharanthus* is not a Rubiaceae
descendant. Required sampling and each minimal-policy exclusion are listed in
`config/phylogeny/fossil_bracketing_taxa.tsv`.

## Three non-interchangeable dating designs

`config/phylogeny/dating_designs.tsv` freezes three separate products. The
2026-07-20 author decision activates the second design for the revised dated
tree because no fossil bracket is available in the compact taxon set:

- revised primary fossil dating uses MCMCTree only after a candidate fossil has
  exact sampled descendants, asset checksums, and an explicit soft-bound model;
- the user-authorized secondary-TimeTree MCMCTree analysis is active only after rows in
  `config/phylogeny/secondary_timetree_constraints.tsv` records the queried
  descendants, TimeTree version, retrieval time, query URL, contributing
  studies, raw response checksum, and the declared transformation into a bound;
  no undocumented fixed point is allowed, and the result is described as
  TimeTree secondary-calibrated rather than fossil-calibrated;
- the submitted `ape::chronos` analysis remains an exact legacy reproduction,
  never primary evidence.

The secondary design is now complete. MCMCTree produced two 50,001-row
posterior chains; the pooled dated tree passed the declared ESS/R-hat,
constraint, checksum, and ultrametric gates. Its pooled root age is 116.6885 Ma.
This completion does not activate a fossil row and does not change the required
wording: the result is TimeTree secondary-calibrated, never fossil-calibrated.

The server audit recovered five legacy fixed points from the surviving script:
113.97, 41.86, 5.91, 2.29, and 25.33 Ma. Corresponding archived TimeTree values
are retained at full precision in
`config/phylogeny/legacy_chronos_calibrations.tsv`. The script sets
`time_min == time_max` and calls `chronos(..., model="clock",
control=chronos.control())`. Exact artifact hashes and the 9-tip versus 17-tip
mismatch are frozen in `config/phylogeny/legacy_chronos_artifacts.tsv`.

This is evidence of what the legacy directory currently contains, not a claim
that it is rerunnable: the current script input is a later 9-tip tree whereas
the archived dated output contains 17 tips, and the TimeTree version, retrieval
date, contributing studies, original script revision, and `ape` version are
missing. Exact legacy reproduction therefore remains blocked until those
inputs reconcile. The historical points are prohibited in revised production
dating.

The manuscript's separate WGD/Ks scaling reports a 111.4--123.9 Ma
*Vitis vinifera*--*Actinidia chinensis* TimeTree range. It is not one of the
five recovered `chronos` script points and remains a labelled legacy secondary
comparison until its own TimeTree query version, retrieval date, contributing
studies, and exact source asset are recovered.

## Fossil-dating workflow

1. Freeze `config/phylogeny/taxa.tsv`,
   `config/phylogeny/minimal_outgroup_design.tsv`,
   `config/phylogeny/public_outgroup_downloads.tsv`, and all asset checksums.
2. Run genome and protein BUSCO with the same lineage and database version.
3. Apply `config/phylogeny/representation_policy.tsv`, then select one primary
   isoform per gene for each declared biological-species representative.
4. Run OrthoFinder 3 on validated genome-derived representative proteomes and
   retain its complete orthogroup, HOG, gene-tree, and species-tree diagnostics
   with checksums. Do not add or graft the ITS-only *Saurauia* tree asset to
   this proteome/count matrix or its downstream nuclear alignment.
5. Build the strict one-copy concatenation locus set and the validated
   multi-copy ASTRAL-Pro locus set using thresholds declared before production.
6. Align each protein locus with MAFFT, back-translate against exact CDS IDs,
   trim codons deterministically, infer IQ-TREE 2 gene trees, and screen
   failures without silently shrinking denominators.
7. Build a partitioned IQ-TREE 2 concatenation tree and an ASTRAL-Pro
   coalescent tree. Compare support, concordance, quartet score, and
   representative-swap sensitivities.
8. Compare the biological-species topology with the assembly-unit diagnostic
   topology, then freeze the accepted biological-species topology.
9. Audit each proposed fossil against the exact sampled tips. Only rows with
   `activation_gate=passed_exact_bracketing_taxa_and_asset_qc` enter a
   fossil-calibrated MCMCTree analysis. When no fossil row is active, the
   author-approved TimeTree secondary constraints may instead activate a
   separately labelled MCMCTree analysis.
10. Estimate divergence times under the declared soft bounds with at least two
    independent chains; check convergence and effective sample sizes.
11. Repeat the analysis under the declared minimal-outgroup and calibration
    sensitivities and report unstable nodes rather
    than hiding disagreement.
12. Run CAFE 5 only after the dated tree and one representative-complement
    family count per biological species reconcile exactly.

Required diagnostics include a prior-only run, leave-one-constraint-out runs,
alternative placements where relevant, and independent
dating chains with convergence and effective-sample-size checks. TimeTree
intervals are copied only through the checksum-bound constraint ledger and are
never introduced as undocumented fixed points.

The nuclear species tree, plastome tree, gene-family analysis, and gene-loss
comparison tree are separate products even when they share some taxa.

All steps share the project-wide maximum of 10 scientific workers. Concurrent
locus jobs and per-job thread counts are budgeted together; no tool receives a
separate 10-thread allowance.
