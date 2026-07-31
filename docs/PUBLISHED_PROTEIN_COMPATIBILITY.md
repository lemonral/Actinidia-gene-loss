# Published-protein compatibility gate

This gate verifies that proteins reconstructed from a chromosome-scope genome
and GFF3 are compatible with the protein FASTA distributed for the same
publisher release. Its purpose is to detect a mismatched genome/GFF/protein
bundle before OrthoFinder, BUSCO-on-derived-proteins, JCVI, or SynOrths uses
that bundle. It does not align proteins, repair gene models, choose isoforms,
or reannotate an assembly.

Run the gate after `extract_primary_annotation.py` has passed its independent
gffread check. The derived input is the resulting `sample.protein.faa`; the
publisher input must be frozen from the same release as the genome and GFF3.
If the publisher FASTA contains multiple isoforms but the derived set contains
one primary transcript per gene, first create and checksum an explicitly
audited publisher-primary subset using the frozen primary-transcript ledger.
When the GFF3 and publisher FASTA expose explicit transcript/protein accession
fields, use `scripts/qc/remap_publisher_primary_proteins.py` and require its
publication gate to pass; see `docs/PUBLISHER_PRIMARY_PROTEIN_REMAP.md`.
That remapper also has a narrow, explicit gene-as-transcript mode for audited
publisher `gene -> CDS` graphs and explicit self-accession sources for bundles
whose publisher protein ID is the transcript/CDS parent. Those modes do not
weaken this compatibility gate: the remapped selected ID set must still be
exact, and every sequence must pass the same residue-level comparison below.
The compatibility program itself never filters either input.

## Exact identifier contract

The first whitespace-delimited token in each FASTA header is the protein ID.
IDs must be unique within each file, and the two ID sets must be exactly equal.
A missing derived ID, a publisher-only isoform, or a duplicated ID rejects the
complete run. Header descriptions are not compared.

Protein sequences are uppercase and are not case-normalized. Embedded
whitespace, empty sequences, unsupported residue symbols, and internal or
repeated stop codons are rejected.

## Only two permitted sequence differences

After ID reconciliation, proteins are compared position by position. The only
permitted differences are:

1. At most one terminal `*` may be present or absent in either sequence. The
   comparison removes that terminal marker and audits whether the two inputs
   differed in terminal-stop representation. Any internal or repeated `*`
   fails.
2. A publisher `X` may match one unambiguous canonical residue
   (`ACDEFGHIKLMNPQRSTVWY`) in the derived sequence. This wildcard is one-way:
   a derived `X` cannot match a concrete publisher residue, and publisher `X`
   cannot rescue a derived ambiguous `B`, `J`, `O`, `U`, or `Z`. Every accepted
   wildcard position is reported using one-based coordinates. Equal `X`/`X`
   positions are exact matches and are not counted as wildcard use.

Lengths must be equal after terminal-stop removal. Every other residue or
length difference rejects the whole run. There is no similarity threshold,
alignment, identifier truncation, or best-hit rescue.

## Command

```bash
python scripts/qc/audit_published_protein_compatibility.py \
  --derived-proteins "$DATA_ROOT/standardized/act_eriantha_hap1_2026/primary_annotation/act_eriantha_hap1_2026.protein.faa" \
  --publisher-proteins "$DATA_ROOT/raw/act_eriantha_hap1_2026/publisher.protein.faa.gz" \
  --sample-id act_eriantha_hap1_2026 \
  --output-dir "$DATA_ROOT/standardized/act_eriantha_hap1_2026/publisher_protein_compatibility"
```

Plain or gzip-compressed FASTA is accepted; gzip is detected from file content.
The wrapper uses one process, starts no worker pool, and refuses an existing
output path.

## Atomic output contract

The command stages the complete directory beside the requested output and
renames it into place only after every ID and sequence passes. On failure it
exits nonzero, removes the staging directory, and publishes no partial audit.
The PASS directory contains:

- `sample.published_protein_compatibility.summary.tsv`: one row with the exact
  ID-set gate, record counts, accepted normalization counts, and input
  SHA-256 values;
- `sample.published_protein_compatibility.records.tsv`: one row per protein,
  raw/normalized lengths, terminal-stop flags, publisher-X wildcard counts and
  one-based positions, mismatch counts, and final status;
- `run_manifest.json`: workflow version, explicit policies, path-free input
  basenames/checksums, execution declaration, and publication gate;
- `checksums.tsv`: byte counts and SHA-256 values for every published file
  except itself.

No output contains an absolute input or output path. The program records only
input basenames and checksums.

## Acceptance and interpretation

A bundle is accepted only when the command exits zero, both `status` and
`publication_gate` are `PASS`, `exact_ID_set` is `true`, both record counts are
equal, the nonpermitted mismatch count is zero, and every checksum verifies.
Publisher-X wildcard counts must be retained in QC reporting; they must not be
described as exact residue identity.

The summary distinguishes raw exact records (`exact_record_count`) from records
that become exact after the explicitly permitted terminal-stop normalization
(`normalized_exact_record_count`). The latter still excludes every record that
uses a publisher-X wildcard.

A failure is evidence that the declared files are not sequence-compatible
under this strict contract. The first response is to verify release identity,
assembly scope, chromosome map, and publisher-primary ID selection. Do not
silently omit failed genes or infer that a blanket de novo reannotation is
required. Any separately justified annotation repair must receive a new
analysis-unit identifier and repeat all QC gates.
