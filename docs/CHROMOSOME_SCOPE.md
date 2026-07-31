# Chromosome-scope materialization

`scripts/qc/materialize_chromosome_scope.py` creates the exact genome and GFF3
scope used by chromosome-based QC, JCVI, orthology, phylogeny, and gene-loss
analysis. It does not edit or delete the downloaded source files.

## Required map

Chromosome selection and renaming must be reviewed before this step. Provide a
three-column, headered TSV in which every retained genome sequence, retained
GFF3 sequence ID, and canonical name occurs once:

```text
genome_seqid	gff_seqid	canonical_seqid
CEVJBR010000001.1	Chr01	Chr01
CEVJBR010000002.1	Chr02	Chr02
```

The separate genome and GFF columns are required because public plant bundles
may use assembly accessions in FASTA while using chromosome labels in GFF3.
Correspondence must come from reviewed assembly metadata and/or validated
chromosome homology; the program never pairs records by order. The old
`source_seqid/canonical_seqid` two-column form remains accepted only for
bundles in which the FASTA and GFF3 sequence IDs are exactly identical.

The explicit three-column form also supports a deliberately empty `gff_seqid`
for a retained genome record that has been audited and is known to have zero
GFF3 feature rows:

```text
genome_seqid	gff_seqid	canonical_seqid
Lco1464_v1.0_h1tg000001l	Lco1464_v1.0_h1tg000001l	Lco1464_v1.0_h1tg000001l
Lco1464_v1.0_h1tg000021l		Lco1464_v1.0_h1tg000021l
```

This is a featureless-record declaration, not a generic missing-data marker.
The genome record remains in the materialized FASTA and receives a retained
audit row with a blank GFF ID and a feature count of zero. If the source GFF3
contains a feature on that genome sequence ID, publication fails. Empty GFF
IDs are prohibited in the legacy two-column form, and a featureless genome ID
may not collide with another row's mapped GFF ID.

The map is the complete inclusion list. Genome records absent from it are
excluded from downstream chromosome-scope assets, but their IDs, lengths, and
annotation feature counts remain in the audit. A mapped genome ID missing from
the genome or a non-empty mapped GFF ID without a feature row is an error.

## Run

```bash
python scripts/qc/materialize_chromosome_scope.py \
  --genome data/raw/sample.genome.fa.gz \
  --gff data/raw/sample.annotation.gff3.gz \
  --seqid-map config/chromosome_maps/sample.tsv \
  --output-dir data/standardized/sample_chromosome_scope \
  --prefix sample
```

The command is single-process and starts no worker threads. Plain and gzip
inputs are recognized by file content. The output genome and GFF3 are written
as reproducible gzip files. The output directory must not already exist.

## Published files

- `<prefix>.genome.fa.gz`: retained FASTA records with canonical first tokens.
- `<prefix>.annotation.gff3.gz`: matching retained features with canonical
  column-1 sequence IDs; retained `##sequence-region` directives are renamed.
- `audit/sequence_scope.tsv`: every input genome sequence, retained/excluded
  status, genome/GFF/canonical IDs, base-pair length, and matched GFF3 feature
  count.
- `audit/feature_counts.tsv`: retained/excluded counts by genome ID, GFF ID,
  canonical ID, and feature type.
- `audit/validation.json`: policy, input identities, reconciliation counts, and
  explicit PASS checks.
- `checksums.tsv`: SHA-256 and byte size for all inputs and published payloads.

Nothing is published unless the entire input passes. Fatal conditions include
duplicate FASTA or map IDs, a non-bijective map, mapped genome IDs or non-empty
mapped GFF IDs missing from their inputs, a feature row on a genome ID declared
featureless, an unmapped GFF3 sequence ID without an identical genome ID,
coordinates outside the explicitly associated genome sequence, broken or
cross-scope `ID`/`Parent` links, conflicting repeated feature IDs, and embedded
GFF3 FASTA. Compatible repeated IDs for multipart features, such as multiple
CDS rows belonging to the same CDS object, are accepted.

Protein and CDS sequences should be extracted from the validated materialized
genome/GFF3 pair in the next pipeline stage. This keeps genome, annotation,
derived proteins, JCVI inputs, and gene-loss coordinates in one explicit
chromosome scope.
