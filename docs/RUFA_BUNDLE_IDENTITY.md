# *Actinidia rufa* bundle identity and chromosome scope

This audit prevents three downloaded *A. rufa* resources from being counted as
independent assemblies merely because their repositories use different names.
It also prevents the inverse error of calling sequence-distinct resources
mirrors without evidence. The frozen machine-readable results are in
`config/rufa_bundle_scope_identity.tsv` and
`config/rufa_bundle_pairwise_identity.tsv`; exact compressed-asset hashes are
in `config/rufa_bundle_asset_checksums.tsv`.

## Exact downloaded resources

| Unit | Publisher identity | Whole genome | Publisher chromosome scope |
|---|---|---:|---:|
| `act_rufa_aru_r1` | `GCA_030159155.1`, ARU_r1.0, Nakamura-B-male | 29 records; 620,324,227 bp | all 29 records; 620,324,227 bp |
| `act_rufa_fuchu` | `GCA_014362265.1`, A.rufaFuchu_1.0 | 501 records; 647,239,177 bp | records `BJWL01000001.1`--`BJWL01000029.1`; 509,119,555 bp |
| `act_rufa_actinidiabase_v1` | ActinidiaBase 2024 v1 bundle associated with `GWHETLF00000000` | 38 records; 615,891,845 bp | `Chr1`--`Chr29`; 613,679,866 bp |

The publication's Supplementary Table S1 reports **100 contigs** and an
assembly size of 615,891,845 bp. That contig statistic describes the assembly
before or during scaffolding; it is not the record count of the distributed
`genome.fa.gz`. The exact frozen release FASTA contains **38 records**: 29
explicitly named pseudochromosomes and nine additional contigs totaling
2,211,979 bp. Its GFF has feature rows on all 29 chromosomes and eight of the
nine extra contigs. `Contig01298` (121,280 bp) is present in the FASTA but has
no GFF feature row. The paper-reported metric and release-file observations are
kept as separate columns in `config/rufa_bundle_scope_identity.tsv`. The
explicit chromosome/GFF map is
`config/chromosome_maps/act_rufa_actinidiabase_v1.publisher_scope.tsv`; the
nine excluded records, lengths, feature counts, and sequence hashes are frozen
in `config/chromosome_maps/act_rufa_actinidiabase_v1.excluded_records.tsv`.

`PubChr01`--`PubChr29` preserve publisher numbering. They do not assert
homology to the same-numbered *A. chinensis* chromosome; final `Chr01`--`Chr29`
labels require the separate chromosome-homology and JCVI procedure.

## Header-independent sequence test

Every FASTA record was read from the frozen gzip file, whitespace was removed,
and nucleotide characters were converted to uppercase. SHA-256 was calculated
for each resulting sequence. The assembly signature is SHA-256 over the sorted
multiset of lines `length<TAB>record_sha256<LF>`, so it is independent of FASTA
headers, record order, line wrapping, and nucleotide case while retaining
duplicate records.

All three pairwise comparisons have zero records with both the same length and
the same normalized sequence SHA-256. Therefore the downloaded ActinidiaBase
v1 genome is not an exact sequence mirror or alias of either ARU_r1.0 or Fuchu;
ARU_r1.0 and Fuchu are also sequence-distinct. This exact test does not measure
relatedness or synteny and does not replace JCVI.

## Biological identity remains a separate question

ARU_r1.0 identifies its biological material as Nakamura-B-male, and Fuchu is a
named accession. The ActinidiaBase publication
([10.1186/s12915-024-02002-z](https://doi.org/10.1186/s12915-024-02002-z))
states only that mature plants were sampled from a germplasm garden in Xi'an;
the downloaded FASTA/GFF does not provide a plant accession name. Consequently:

- sequence-level non-identity with ARU and Fuchu is resolved;
- the ActinidiaBase plant-level accession is recorded as unresolved;
- ActinidiaBase must not be treated as an independent biological replicate
  until passport or BioSample metadata identifies the plant;
- the ActinidiaBase and GWH labels for the 2024 publication remain one release
  family, but exact cross-repository file identity has not been tested because
  only the ActinidiaBase bundle was included in this three-bundle audit.

Assembly QC may evaluate all three sequence-distinct bundles. Species-level
trees, gene-family counts, gene-loss calls, and PGLS must select one declared
representative and must never treat alternative assemblies as replicate
species observations.
