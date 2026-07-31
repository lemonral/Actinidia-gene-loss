# *Actinidia deliciosa* A--F bundle partitioning

## Purpose

The Qinmei 2025 and ADM 2026 downloads are each one hexaploid bundle, not one
analysis unit. The gene-loss, BUSCO, JCVI, expression, copy-number, and spatial
workflows require six visible assembly units (A--F). The splitter therefore
uses the same partition for the genome, GFF3, CDS, and protein set and refuses
to publish a usable manifest when their identifiers do not close exactly.

The two bundles are independent individuals and are never combined. Each is
split and evaluated separately; QC and JCVI determine which complete A--F set
enters the main analysis. The other set remains a sensitivity candidate.

## Explicit mapping contract

`config/deliciosa_polyploid_partitions.tsv` contains six rows per parent
bundle. Every row freezes the source token, displayed A--F label, biological
individual, derived `assembly_unit_id`, expected chromosome count, and whether
known unplaced annotation records are permitted.

- Qinmei uses chromosome suffix `_1` through `_6`. Chromosome and annotation
  identifiers must independently imply the same token.
- ADM uses the genome header fields `OriSeqID=ChrNNhN` and `Chromosome ANN`.
  When both fields are present they must agree. The resulting accession-to-unit
  map assigns GFF3 sequence IDs and annotation `Position=` accessions exactly.
  The `AchdmhNcNN...` record identifier and `OriSeqID` provide independent
  cross-checks. ADM CDS accessions (`GWHT`) and protein accessions (`GWHP`)
  occupy different namespaces, so equality is tested through reciprocal
  header links (`Protein=` and `mRNA=`) and their exact GFF3
  `Accession`/`Protein_Accession` values. Their raw record IDs are not expected
  to be identical.

These rules are source-specific by design. A changed publisher header is an
error requiring inspection and an explicit mapping update, not a reason to
guess a partition.

## Run

First create the verified parent manifest with
`scripts/qc/resolve_asset_manifest.py`, requiring all four roles. Then run one
bundle at a time:

```bash
python scripts/qc/split_deliciosa_polyploid.py \
  --resolved-manifest /path/to/resolved_public_candidates.tsv \
  --bundle-id act_deliciosa_qinmei_2025 \
  --mapping config/deliciosa_polyploid_partitions.tsv \
  --output-dir /path/to/standardized/deliciosa_qinmei_2025_A_to_F

python scripts/qc/split_deliciosa_polyploid.py \
  --resolved-manifest /path/to/resolved_public_candidates.tsv \
  --bundle-id act_deliciosa_adm_2026 \
  --mapping config/deliciosa_polyploid_partitions.tsv \
  --output-dir /path/to/standardized/deliciosa_adm_2026_A_to_F
```

The output directory must not exist. Input paths may be relative to the parent
manifest. The program recomputes and checks every declared input SHA-256 before
creating a staging directory. It is single-process and streams sequence data.

## Validation and outputs

On a passing run, `resolved_assembly_units.tsv` contains exactly six portable
rows. Its asset paths are relative to the manifest and are accepted directly
by the basic-statistics and BUSCO wrappers. All rows remain QC candidates;
`include_gene_loss=false` until assembly selection is complete.

The `audit` directory contains:

- `checksums.tsv`: verified parent hashes plus every derived and retained-file
  SHA-256;
- `fasta_records.tsv`: assignment, length, and normalized-sequence SHA-256 for
  every genome, CDS, and protein record;
- `partition_summary.tsv`: A--F chromosome, feature, CDS, protein, and exact
  GFF3-identifier counts;
- `id_join_exceptions.tsv`: every exact CDS/protein/GFF3 set difference;
- `unassigned_records.tsv`: each retained unmatched or unplaced sequence/seqid;
- `unassigned/*.gz`: the corresponding source records or GFF3 rows, retained
  verbatim at record level rather than silently discarded;
- `validation_issues.tsv` and `run_metadata.json`: the publication gate.

For Qinmei, recognized `scf` annotations are allowed only in the retained
unassigned scope. They are not attached to A--F because the matching scaffold
genome FASTA is absent from the public 174-chromosome download. An unexpected
identifier, an ADM unplaced record, conflicting partition evidence, a missing
chromosome, coordinate overflow, duplicate FASTA ID, CDS/protein difference,
or non-exact GFF3 join blocks the run.

A blocked run is installed for inspection with status `BLOCKED` and no
`resolved_assembly_units.tsv`; the command exits 2. A malformed input or
checksum mismatch installs nothing and exits 1. This distinction prevents
partial partitions from entering downstream QC while preserving evidence
needed to correct a source mapping.
