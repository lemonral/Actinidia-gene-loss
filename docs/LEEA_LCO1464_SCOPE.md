# *Leea coccinea* Lco1464 scope

## Evidence boundary

The exact public bundle is Zenodo record
[13362874](https://doi.org/10.5281/zenodo.13362874):
`Lco1464_v1.0.fasta` (1,104,911,269 file bytes; publisher MD5
`2ea922e9006b2053c291c534d1405cef`) and
`Lco1464_v1.0_gene_annotation_ver1.0.gff3` (102,108,962 bytes;
publisher MD5 `ed8f7699487a95e750490bfa3cf1514a`). The associated
[publication](https://doi.org/10.1038/s41467-025-61387-9) reports one sampled
accession, *L. coccinea* 1464, one HiFi library from that accession, and
hifiasm assembly of the HiFi reads. The two `h1` and `h2` components are
therefore phased assemblies of the same biological individual, not two
species replicates.

Supplementary Data 3 and 5 spell the epithet `cuccinea`, whereas the article
title, main text, dataset context, and accession label use *Leea coccinea*.
The project treats the supplementary spelling as a typographical variant,
normalizes the reader-facing taxon to *L. coccinea*, and preserves `Lco1464`
as the accession/assembly identity.

The publication's Supplementary Data 5 reports 63 h1 sequences totaling
554,723,631 bp and 77 h2 sequences totaling 550,183,858 bp. These values
reconcile exactly to the 140 records and 1,104,907,489 sequence bases in the
downloaded FASTA. The paper calls them sequences, gives no chromosome-to-
sequence assignment, and describes no Hi-C scaffolding for Lco1464. Every
record is named as a hifiasm contig (`Lco1464_v1.0_h1tg...` or
`Lco1464_v1.0_h2tg...`), not as a chromosome. Consequently, this project must
not infer chromosomes from record order, record length, the `l`/`c` suffix,
or a size threshold. In particular, the fact that each haplotype has twelve
records at least 20 Mb long is not evidence that those twelve records are a
complete chromosome set.

## Frozen record inventory

`config/phylogeny/leea_lco1464_records.tsv` is the exact record-level audit of
the downloaded FASTA/GFF pair. It freezes every source sequence ID, its
haplotype, sequence length, GFF feature and gene counts, and whether the GFF
contains any row on that sequence. Source and canonical IDs are deliberately
identical: assigning `Chr` labels would create unsupported biological claims.

The exact three-column inputs for the scope materializer are
`config/chromosome_maps/leea_coccinea_lco1464_h1.full_haplotype.tsv` and
`config/chromosome_maps/leea_coccinea_lco1464_h2.full_haplotype.tsv`. Despite
the shared `chromosome_maps` directory name, these files retain all publisher
contigs and do not claim that any contig is a chromosome. A blank `gff_seqid`
is permitted only for a registry row declared to have zero GFF features; an
unexpected feature on such a sequence makes materialization fail.

The complete GFF contains no sequence ID absent from the FASTA and no feature
end beyond its sequence length. Per-haplotype annotation coverage is:

| Haplotype | Sequences | Assembly bp | Sequences with GFF rows | bp on sequences with GFF rows | Fraction of assembly bp | Gene loci |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| h1 | 63 | 554,723,631 | 46 | 553,353,975 | 99.7531% | 30,024 |
| h2 | 77 | 550,183,858 | 60 | 548,536,384 | 99.7006% | 25,109 |

The seventeen records without GFF features in each haplotype remain part of
that haplotype's genome scope. They are not deleted merely because they lack
gene models. They contribute to whole-haplotype assembly statistics and
genome-mode BUSCO, but naturally contribute no model to protein extraction.

## Materialization and analysis policy

1. Split the combined source FASTA and GFF by exact record membership in the
   frozen inventory. Do not use an inferred numeric cutoff. The h1 bundle must
   contain exactly 63 FASTA records and all GFF rows on the 46 annotated h1
   records; h2 must contain exactly 77 FASTA records and all GFF rows on the
   60 annotated h2 records. Preserve source sequence IDs.
2. Run basic assembly statistics and genome BUSCO on every sequence in each
   complete haplotype. Run strict primary-transcript extraction independently
   on each paired haplotype genome/GFF and run protein BUSCO on each resulting
   proteome. A BUSCO result from the combined h1+h2 diploid file cannot choose
   a representative because duplicated homologues inflate duplication.
3. The production biological-species OrthoFinder matrix, species tree,
   MCMCTree input, and CAFE count matrix may contain exactly one Lco1464
   haplotype. H1 is only the provisional preference because it has fewer
   contigs and more annotated loci; it is not accepted until per-haplotype
   genome/protein BUSCO, extraction integrity, contamination checks, and
   orthologue occupancy pass. Do not select a haplotype because it yields a
   preferred topology.
4. Retain both h1 and h2 with upright suffixes in the assembly-unit diagnostic
   tree. If h1 is primary, rebuild the full species analysis with h2 as the
   predeclared representative-swap sensitivity. Never count h1 and h2 as two
   biological species or sum their family counts for CAFE.
5. Lco1464 is a fossil-bracketing phylogeny outgroup, not an *Actinidia*
   gene-loss comparison unit. It does not enter JCVI chromosome coverage,
   SynOrths loss calling, spatial loss analyses, shared/non-shared loss
   denominators, or PGLS.

This full-haplotype policy is an explicit exception to chromosome-only
processing. Removing smaller Lco1464 contigs would be an undocumented
length-based filter, whereas retaining them cannot create extra species tips
and preserves the publisher assembly for genome and annotation QC.
