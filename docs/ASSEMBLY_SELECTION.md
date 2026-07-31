# Assembly selection and biological aggregation

The source-backed discovery inventory is
`config/candidate_assembly_catalog.tsv`, documented in
`docs/CANDIDATE_ASSEMBLY_CATALOG.md`. That catalog records alternatives and
exclusions; it does not authorize analysis inclusion. Only rows promoted into
`config/assemblies.tsv` after the gates below belong to the executable cohort.

## General gates

An assembly unit can enter the gene-loss candidate cohort only after all of the
following pass:

1. accession, version, biological sample, ploidy, and haplotype identity;
2. genome/GFF/protein scope and sequence-ID closure;
3. contamination and basic assembly statistics;
4. genome and primary-protein BUSCO with the same local lineage and version;
5. bidirectional JCVI/synteny coverage with the reference;
6. clean SynOrths input and output reconciliation;
7. a declared callable reference-gene denominator.

## *Actinidia eriantha*

HAP1 and HAP2 are two haplotypes from one diploid individual. They are run as
separate assembly units. At the biological-species level:

- loss supported by both haplotypes: species-supported loss;
- loss supported by one haplotype only: haplotype-specific absence;
- failure of synteny/callability or conflicting sequence evidence: uncertain.

The two haplotypes must not be counted as two independent species or two
biological replicates.

The production source is the complete six-file Figshare v1 deposit
(`10.6084/m9.figshare.31034008.v1`), which pairs a genome, GFF, and published
protein set for each haplotype. Each deposited genome FASTA contains exactly 29
pseudochromosomes: 636,551,830 bp for HAP1 and 628,273,856 bp for HAP2. Its 29
sequence identifiers match the GFF exactly, with no unmatched or unplaced
records in the public files.

The paper also reports pre-anchoring primary-contig totals of 663.76 Mb and
633.05 Mb. Those paper-level totals are retained as reported metadata, not
silently merged with the smaller deposited chromosome scope. Genome BUSCO and
all downstream comparative results state explicitly which scope they assess.

The NCBI chromosome-only FASTAs are not interchangeable with the Figshare GFFs.
That mixed-source combination produced chromosome-block internal-stop failures
and is quarantined as a failed diagnostic, not used as evidence for
reannotation or as an OrthoFinder input.

## *Actinidia deliciosa*

The replacement must preserve the six-haplotype/174-chromosome biological
design. A diploid assembly or one collapsed representative is not a valid
replacement for A-F.

Two independent candidates are tracked:

- the 2025 Qinmei bundle, whose 174 chromosome records can be split by the
  `_1` to `_6` suffix but whose unplaced FASTA records are not public; and
- the 2026 GWH ADM bundle, which is a different individual/assembly and must
  not be mislabelled as a full Qinmei release.

The exact genome/GFF/CDS/protein partition is validated before six derived rows
are frozen. The candidate with stronger paired scope, BUSCO, annotation closure,
and synteny is used for the main rerun; the other remains a sensitivity result.

Both current splits pass their structural 29-chromosome-per-unit gate. They are
not yet interchangeable:

- Qinmei's chromosome FASTA has no scaffold records, while its GFF contains
  2,031 `scf` sequence IDs and 136,516 feature rows. The corresponding 13,623
  CDS/protein records are retained in the outside-chromosome audit and cannot
  be added to A--F without the missing genome sequences. Publisher proteins in
  A, C, and F each contain one severe internal-stop model; a strict run must
  fail or retain a separately labelled diagnostic omission, never hide it.
- ADM has no unassigned genome or annotation records and no publisher-protein
  internal stops. It is multi-isoform and uses different GWH identifier
  namespaces for transcripts, CDS, and proteins. Publisher compatibility must
  therefore use an exact one-to-one GFF-attribute mapping to the selected
  primary transcripts; direct comparison of the unscreened full sets is
  invalid.

This makes ADM the stronger provisional candidate, not the selected winner.
Selection waits for strict extraction/gffread, mapped publisher compatibility,
per-A--F genome and derived-protein BUSCO, chromosome-homology assignment,
JCVI/callable coverage, and orthogroup occupancy. Whole-bundle BUSCO values are
not substituted for the twelve per-unit assessments.

## *Actinidia rufa*

ARU_r1.0, Fuchu, and the 2024 ActinidiaBase v1 bundle are inventoried. The
exact sequence audit shows that the three downloaded genomes are not exact
mirrors: every pair has zero complete records sharing both length and normalized
sequence SHA-256. ActinidiaBase v1 has 38 records, specifically 29
pseudochromosomes and nine extra contigs. The paper's Supplementary Table S1
reports 100 assembly contigs, but 100 is not the record count of the released
`genome.fa.gz`; both metrics are retained with their contexts. The sampled
ActinidiaBase plant is unnamed in
the publication, so sequence non-identity is resolved while biological-accession
independence remains unresolved. Alternative assemblies are never biological
replicates. If no candidate produces adequate chromosome-matched/callable
coverage, *A. rufa* is excluded from downstream gene-loss comparisons. It may
still be retained in a species tree if protein BUSCO and ortholog occupancy
pass, because phylogenetic inclusion and gene-loss callability are different
questions. Exact evidence and scope maps are documented in
`docs/RUFA_BUNDLE_IDENTITY.md`.

The completed candidate comparison is recorded in
`results/qc/actinidia_rufa_selection/assembly_selection.tsv`. ActinidiaBase v1
is the provisional primary candidate: its exact current-input JCVI coverage is
99.426289% on the *C. scandens* reference side and 98.638443% on the *A. rufa*
side. Chromosome-based inputs contain only its 29 pseudochromosomes; all nine
unplaced contigs and their 1,642 GFF3 feature rows remain in the full-release
QC audit and are excluded from JCVI, spatial analysis, and gene-loss calling.
The exact-bound ARU result is retained as an assembly-swap sensitivity
(98.768210% and 86.126011%, respectively). Fuchu is excluded from the
production candidate ranking because five primary genes map to multiple
publisher proteins; its older JCVI values used a different input set and are
reported only as historical context. Final admission to the gene-loss cohort
still requires clean SynOrths, translated-genome evidence, and a callable
reference-gene denominator. The unnamed ActinidiaBase plant remains one
biological species representative, not an independent replicate of ARU or
Fuchu.
