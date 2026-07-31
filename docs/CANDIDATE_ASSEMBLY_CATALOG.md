# Candidate assembly catalog

`config/candidate_assembly_catalog.tsv` is the public, source-backed inventory
of assembly releases considered during the 2026 revision. It is a discovery and
decision record, not the executable analysis cohort. A row in this catalog never
enters the gene-loss analysis automatically.

The executable cohort remains `config/assemblies.tsv`. Moving a candidate there
requires frozen file versions and checksums, genome/GFF/protein scope closure,
basic statistics, genome and protein BUSCO, bidirectional JCVI coverage, clean
orthology reconciliation, and an explicit callable denominator.

## Classification

- `primary_replacement_candidate`: a complete genome/GFF/protein bundle that is
  eligible for full QC as a possible main-analysis replacement. It is not yet an
  inclusion decision.
- `sensitivity`: a complete secondary bundle, a ploidy-mismatched comparison,
  a population accession, or a repository mirror that must remain separate from
  the primary comparison.
- `unavailable`: the exact runnable bundle is incomplete. DNA-only releases
  remain unavailable for gene-loss calls even when a paper reports an annotation.
- `excluded`: a known hybrid or other non-independent unit that must not be
  treated as a species-level sample.

The four `public_*` fields describe assets verified at the repository snapshot
on 2026-07-17. `accession_version_status` distinguishes a versioned accession
from a GWH accession root whose current file version and checksums still need to
be frozen. `alternate_accessions` and `deduplication_group` prevent repository
mirrors or duplicate submissions from inflating the sample count.

## Biological independence rules

Assembly units and biological samples are deliberately separate:

- The public *A. eriantha* HAP1 and HAP2 files are two haplotypes from the same
  biological accession. They can be processed separately, but they are not two
  species or independent biological replicates. The downloaded releases contain
  29 chromosome-anchored pseudochromosomes per haplotype and omit the unanchored
  primary-contig scope described in the paper.
- Each pair in the *A. eriantha* diversity panel and the *Actinidia* pan-genome
  projects likewise represents two haplotypes of one accession.
- *A. zhejiangensis* A and B are parental haplomes of the same F1 hybrid,
  appropriately written *A. × zhejiangensis*, not independent species. Both are
  excluded from species-level loss comparisons.
- HH01 in the 2026 *A. eriantha* panel is a hybrid accession. Both HH01
  haplotypes are excluded.
- The two NCBI Assembly identifiers listed for *A. chinensis* Guimi No. 2 refer
  to the same 607,462,286-bp submitted assembly. The catalog records one paper-
  linked WGS master and the two identifiers as aliases; the assembly is counted
  once. NCBI does not provide a matched structural annotation for these records,
  so the assembly remains unavailable for this gene-loss workflow.
- ActinidiaBase and GWH labels associated with one 2024 publication are kept in
  one release family unless exact cross-repository evidence proves otherwise.
  Repository labels are not biological replicates. For *A. rufa*, the downloaded
  ActinidiaBase genome is sequence-distinct from ARU and Fuchu; the identity of
  the unnamed ActinidiaBase plant remains unresolved.

## Availability decisions

The 2024 *A. chinensis* graph panel and the 2025 genus super-pan-genome provide
official GWH genome DNA, but GWH does not provide the exact matching annotation
assets used by the publications. These rows therefore fail closed as
`unavailable` for gene-loss analysis until version-matched GFF and protein
files are obtained.

The 2026 *A. eriantha* diversity panel provides DNA, GFF, CDS, and protein files
for the 11 non-hybrid wild accessions and is retained as a sensitivity panel.
Kuimi and ACM4 are annotated tetraploid, 116-chromosome sensitivity candidates.
Neither can replace a diploid *A. chinensis* primary unit without changing the
biological design.

The chromosome-level *A. rufa* MT570001 assembly has been published, but no
public accession or retrievable genome/annotation bundle was found by the
snapshot date. Its paper metrics must not be substituted for sequence files.
ARU and Fuchu remain replacement candidates, conditional on assembly-scope and
bidirectional JCVI gates. The exact downloaded ActinidiaBase v1 genome contains
29 pseudochromosomes plus nine contigs and shares no exact normalized sequence
record with either ARU or Fuchu. It remains a sensitivity candidate because the
publication does not name the sampled Xi'an germplasm-garden plant; sequence
identity and biological-accession identity are deliberately separate. See
`docs/RUFA_BUNDLE_IDENTITY.md`.

For hexaploid *A. deliciosa*, only a candidate that preserves the six A-F
assembly units can replace the original design. Qinmei and ADM are independent
full-bundle candidates and must be evaluated separately; they must never be
combined into one biological accession.

## Public sources

The catalog points to accession-level NCBI or GWH records and to the following
primary publications:

- *A. arguta* M1: [10.1016/j.xplc.2024.100856](https://doi.org/10.1016/j.xplc.2024.100856)
- *A. chinensis* graph genomes: [10.1002/advs.202400322](https://doi.org/10.1002/advs.202400322)
- genus super-pan-genome: [10.1093/hr/uhaf067](https://doi.org/10.1093/hr/uhaf067)
- *A. eriantha* HAP1/HAP2 release: [10.1038/s41597-026-07414-w](https://doi.org/10.1038/s41597-026-07414-w)
- *A. eriantha* diversity panel: [10.1186/s13059-026-04068-0](https://doi.org/10.1186/s13059-026-04068-0)
- ActinidiaBase/GWH assemblies: [10.1186/s12915-024-02002-z](https://doi.org/10.1186/s12915-024-02002-z)
- *A. × zhejiangensis* hybrid evidence: [10.1111/tpj.16336](https://doi.org/10.1111/tpj.16336)
- *A. rufa* MT570001: [10.1007/s44281-026-00103-z](https://doi.org/10.1007/s44281-026-00103-z)

Run the focused catalog validation with:

```bash
python -m unittest tests.config.test_candidate_assembly_catalog
```

The validation checks schema and controlled values, public-safe paths, fail-
closed asset classification, paired-haplotype grouping, hybrid exclusions,
chromosome-subset disclosure, and duplicate-accession handling.
