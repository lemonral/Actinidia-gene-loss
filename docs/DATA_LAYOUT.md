# Data layout and provenance policy

## Why data are outside Git

Assemblies, GFF/GTF files, proteins, CDS, raw reads, BLAST output, BUSCO run
directories, OrthoFinder results, and NLR-Annotator work directories are too
large for a normal Git repository. They may also have redistribution terms
that differ from those of this code.

The server layout is:

```text
../actinidia_gene_loss_data/
├── downloads/       # immutable downloaded archives and API metadata
├── raw/             # extracted, source-faithful files or legacy soft links
├── standardized/    # renamed and validated analysis inputs
├── checksums/       # source and locally calculated checksums
├── work/            # disposable tool working directories
├── results_large/   # large accepted outputs
└── logs/            # download and compute logs
```

`raw` is append-only. A standardized input must record its source file,
source checksum, transformation command, output checksum, record count, and
sequence-ID contract.

## Legacy links

Previously downloaded outgroup data are linked under
`raw/phylogeny_outgroups/legacy_linked_pending_provenance`. The name is
intentional: a working symlink proves file availability, not taxonomic
identity, accession, annotation compatibility, or publication suitability.

The 22 manuscript-era *Actinidia* assembly units are declared in
`config/legacy_analysis_units.tsv`. Run
`scripts/migration/migrate_legacy_assets.py` against the private historical
manifest to copy assets at or below the declared threshold and soft-link
larger files under `legacy_linked/<real_assembly_unit_id>/`. The generated
checksum report and compatibility QC manifest live in this external data
store. They contain runtime paths and are never committed.

## Optional download proxy

The downloader deliberately removes inherited `HTTP_PROXY`, `HTTPS_PROXY`,
and related variables so that routing is explicit and reproducible. Foreign
repositories may be declared with repeatable `--proxy-domain` arguments and a
user-local Mihomo endpoint. Stable NCBI FTP HTTPS paths may additionally use
`--segmented-proxy-domain` to preserve aria2 resume metadata and bounded
multi-connection transfer. Figshare remains curl-routed because its signed
redirect must be refreshed through the original URL. The signed S3 transfer
also remains proxied when direct S3 throughput is poor.

CNCB/GWH and other domestic hosts remain direct simply by omitting them from
the proxy-domain list. No credentials, subscription URL, or node name is
stored in the repository, and no machine-wide proxy setting is changed.

## Git boundary

Only manifests, code, tests, small summary tables, plot data, and final figures
belong in Git. Absolute server paths, manuscripts, supplementary workbooks,
reviewer responses, private checksums, credentials, and unpublished raw data do
not belong in the repository.
