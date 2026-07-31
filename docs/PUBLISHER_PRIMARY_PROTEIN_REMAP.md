# Publisher-primary protein ID remapping

Some matched publisher bundles use different identifier namespaces in their
GFF3 and protein FASTA. The ADM *Actinidia deliciosa* release is one example:
the GFF3 transcript ID is an `Achdm...t1` value, while the protein FASTA
first-token ID is a `GWHP...` accession. The publisher FASTA header and GFF3
provide an explicit bridge:

```text
GFF3 mRNA: ID=Achdmh1c01g00001.t1;Accession=GWHTJJEQ000001.1
GFF3 CDS:  Parent=Achdmh1c01g00001.t1;Protein_Accession=GWHPJJEQ000001.1
protein:   >GWHPJJEQ000001.1  mRNA=GWHTJJEQ000001.1  OriID=Achdmh1c01g00001.t1
```

`scripts/qc/remap_publisher_primary_proteins.py` is the only production
adapter for that situation. It runs after primary-isoform selection and before
`audit_published_protein_compatibility.py`. It is an identifier remapper and
subsetter, not a sequence comparison, best-hit matcher, or annotation repair.

## Fail-closed mapping contract

The selected-primary protein FASTA defines the exact selected transcript-ID
set through its unique first-token IDs. Before subsetting anything, the tool
validates the complete source bundle:

1. every accepted GFF3 transcript/self-gene ID and declared transcript
   accession is unique;
2. every coding model has exactly one publisher accession across all CDS
   parts, from the configured case-sensitive attribute or explicit sole-parent
   source;
3. every protein accession maps to exactly one transcript/self-gene model;
4. the complete GFF3 coding-protein accession set equals the complete
   publisher FASTA first-token ID set, with no missing or extra record;
5. in `metadata` mode, each publisher header `OriID` equals its GFF3 transcript
   ID and each header `mRNA` equals its GFF3 transcript `Accession`; in
   `first_token` mode, the publisher first token is the protein accession and
   the complete one-to-one transcript bridge must come from declared GFF3
   attributes;
6. every selected-primary transcript has exactly one validated publisher
   protein mapping;
7. the remapped output ID set and order exactly equal the selected-primary
   input, with no duplicate, missing, or extra output record;
8. every output sequence is character-for-character equal to its source
   publisher sequence.

Publisher proteins belonging to valid but non-selected isoforms are expected.
They are not "extra" source records because they close exactly to GFF3 coding
transcripts. They are excluded from the primary subset and retained in the
source inventory. A publisher record with no GFF3 coding accession, or a GFF3
coding accession with no publisher record, is an actual closure error and
rejects the complete run.

No identifier is truncated, inferred from a prefix, or rescued by similarity.
No protein residue is uppercased, trimmed, stop-normalized, aligned, or
replaced. Internal/repeated stops, lowercase or unsupported symbols, malformed
headers, and duplicate identifiers fail before publication.

## ADM command

The defaults exactly name the audited ADM fields (`ID`, `Accession`, `Parent`,
`Protein_Accession`, `OriID`, and `mRNA`):

```bash
python scripts/qc/remap_publisher_primary_proteins.py \
  --selected-primary-proteins "$DATA_ROOT/standardized/deliciosa_adm_2026_A_to_F/assembly_units/act_deliciosa_adm_2026_A/primary_annotation_strict/act_deliciosa_adm_2026_A.protein.faa" \
  --gff "$DATA_ROOT/standardized/deliciosa_adm_2026_A_to_F/assembly_units/act_deliciosa_adm_2026_A/act_deliciosa_adm_2026_A.gff3.gz" \
  --publisher-proteins "$DATA_ROOT/standardized/deliciosa_adm_2026_A_to_F/assembly_units/act_deliciosa_adm_2026_A/act_deliciosa_adm_2026_A.protein.fa.gz" \
  --sample-id act_deliciosa_adm_2026_A \
  --output-dir "$DATA_ROOT/standardized/deliciosa_adm_2026_A_to_F/assembly_units/act_deliciosa_adm_2026_A/publisher_primary_remap"
```

Plain or gzip-compressed inputs are accepted; gzip is detected from file
content. For a release with the same explicit mapping model but different
attribute names, freeze the schema in the command with
`--transcript-id-attribute`, `--transcript-accession-attribute`,
`--cds-parent-attribute`, `--protein-accession-attribute`,
`--publisher-transcript-key`, and `--publisher-mrna-accession-key`. Do not use
those options to guess a relationship that the publisher did not declare.

## Explicit gene-as-transcript graph

Some matched publisher annotations contain `gene` and `CDS` rows but no
accepted `mRNA` or `transcript` row. This is a graph representation choice, not
permission to guess a transcript ID. Run primary extraction with its separate
`--gene-as-transcript` gate first. The remapper has the same explicit option
and otherwise rejects a gene-to-CDS graph.

In this mode the complete remap GFF3 must contain zero accepted transcript
rows. Each declared gene ID is preserved as its self-transcript ID. Every CDS
must name exactly one declared gene parent; CDS attached to a declared
`pseudogene` are rejected. All parts of a coding gene must declare the same
single publisher protein accession. Missing, comma-multiple, or conflicting
accessions reject the run. The full coding gene set is checked, including
non-selected genes, so ambiguity outside the selected chromosome-primary set
cannot pass silently.

There are two explicit accession-source choices; neither performs inference:

- `--transcript-accession-source transcript_id` records the unchanged
  transcript/self-gene ID when no separate transcript-accession field exists.
- the protein accession normally comes from the case-sensitive attribute
  named by `--protein-accession-attribute`. If the publisher protein first
  token is exactly the sole CDS `Parent`, use
  `--protein-accession-source cds_parent`.

`--gene-as-transcript` refuses a GFF3 that contains even one accepted
transcript row. It therefore cannot rewrite or partially mix a normal
gene-transcript-CDS graph.

### *A. rufa* ARU_r1.0 command

The Kazusa ARU GFF3 uses the gene ID as the sole CDS `Parent`, and the same ID
is the publisher protein FASTA first token:

```bash
python scripts/qc/remap_publisher_primary_proteins.py \
  --selected-primary-proteins "$DATA_ROOT/standardized/act_rufa_aru_r1_publisher_scope/primary_annotation_strict/act_rufa_aru_r1.protein.faa" \
  --gff "$DATA_ROOT/downloads/Actinidia_rufa/ARU_r1.0/ARU1.0.genes.gff.gz" \
  --publisher-proteins "$DATA_ROOT/downloads/Actinidia_rufa/ARU_r1.0/ARU1.0.proteins.fasta.gz" \
  --sample-id act_rufa_aru_r1 \
  --output-dir "$DATA_ROOT/standardized/act_rufa_aru_r1_publisher_scope/publisher_primary_remap" \
  --gene-as-transcript \
  --transcript-accession-source transcript_id \
  --protein-accession-source cds_parent \
  --publisher-header-mode first_token
```

### *A. rufa* Fuchu command

The NCBI Fuchu GFF3 attaches CDS directly to gene IDs and declares the
publisher accession in the case-sensitive CDS `protein_id` attribute:

```bash
python scripts/qc/remap_publisher_primary_proteins.py \
  --selected-primary-proteins "$DATA_ROOT/standardized/act_rufa_fuchu_publisher_scope/primary_annotation_strict/act_rufa_fuchu.protein.faa" \
  --gff "$DATA_ROOT/downloads/Actinidia_rufa/Fuchu/GCA_014362265.1_A.rufaFuchu_1.0_genomic.gff.gz" \
  --publisher-proteins "$DATA_ROOT/downloads/Actinidia_rufa/Fuchu/GCA_014362265.1_A.rufaFuchu_1.0_protein.faa.gz" \
  --sample-id act_rufa_fuchu \
  --output-dir "$DATA_ROOT/standardized/act_rufa_fuchu_publisher_scope/publisher_primary_remap" \
  --gene-as-transcript \
  --transcript-accession-source transcript_id \
  --protein-accession-attribute protein_id \
  --publisher-header-mode first_token
```

Use the complete matched source GFF3 in the remap, not the chromosome-filtered
GFF3, because exact source closure must also audit publisher proteins belonging
to Fuchu unplaced records. Only the already selected chromosome-primary IDs
are written to the remapped output.

The frozen Fuchu bundle currently fails this exact command before publication.
Five direct-gene models contain conflicting `protein_id` values across their
CDS rows: `gene-Acr_00g0080410` (unplaced), `gene-Acr_01g0014530`,
`gene-Acr_07g0010490`, `gene-Acr_13g0008400`, and
`gene-Acr_14g0001000`. The latter four are in chromosome scope. A
gene-as-transcript concatenation can produce one synthetic sequence for such a
gene, but it cannot establish a one-to-one mapping to two publisher proteins.
There is deliberately no option to drop or choose one accession. Fuchu remains
blocked at the publisher-compatibility gate unless a separately versioned,
biologically justified model-resolution stage is implemented and re-audited.

## NCBI first-token command

Matched NCBI bundles commonly put the protein accession in the publisher
FASTA first token and declare the transcript/protein bridge only in GFF3, for
example `ID=rna-...;orig_transcript_id=...` on `mRNA` rows and
`Parent=rna-...;protein_id=KAI...` on `CDS` rows. Use the explicit
`first_token` mode for this schema:

```bash
python scripts/qc/remap_publisher_primary_proteins.py \
  --selected-primary-proteins primary_annotation_strict/sample.protein.faa \
  --gff matched_full_source_annotation.gff.gz \
  --publisher-proteins matched_full_source_protein.faa.gz \
  --sample-id sample \
  --output-dir publisher_primary_remap \
  --transcript-accession-attribute orig_transcript_id \
  --protein-accession-attribute protein_id \
  --publisher-header-mode first_token
```

The full matched source GFF3 and publisher protein FASTA are used for source
accession closure even when the selected-primary set was derived from a
chromosome-only scope. In this mode the two publisher-header mapping columns
are empty and `publisher_header_mapping_check` is `not_applicable`; exact GFF3
protein-accession closure, one-to-one mapping, selected-ID closure, and
character-for-character sequence preservation remain mandatory. This mode is
not a license to infer accessions from descriptions or sequence similarity.

## Explicit self-accessions in a normal transcript graph

The ActinidiaBase v1 *A. rufa* GFF3 has normal `gene -> mRNA -> CDS` rows, but
its mRNA rows have no separate transcript accession and its publisher protein
first token equals the mRNA ID/CDS `Parent`. Keep normal transcript mode and
declare those identities explicitly:

```bash
python scripts/qc/remap_publisher_primary_proteins.py \
  --selected-primary-proteins "$DATA_ROOT/standardized/act_rufa_actinidiabase_v1_publisher_scope/primary_annotation_strict/act_rufa_actinidiabase_v1.protein.faa" \
  --gff "$DATA_ROOT/downloads/Actinidia_rufa/ActinidiaBase_v1/Actinidia_rufa_v1_gene_model.gff3.gz" \
  --publisher-proteins "$DATA_ROOT/downloads/Actinidia_rufa/ActinidiaBase_v1/protein.fa.gz" \
  --sample-id act_rufa_actinidiabase_v1 \
  --output-dir "$DATA_ROOT/standardized/act_rufa_actinidiabase_v1_publisher_scope/publisher_primary_remap" \
  --transcript-accession-source transcript_id \
  --protein-accession-source cds_parent \
  --publisher-header-mode first_token
```

This is not `--gene-as-transcript`: accepted mRNA rows must be present, and the
normal transcript graph is retained. The explicit source modes avoid treating
one GFF3 attribute as two independent fields and preserve the default
attribute-based contract for existing commands.

With the frozen files audited on 2026-07-17, ARU closed exactly at 63,947
source/selected records and all 63,947 remapped sequences passed the separate
publisher-compatibility gate exactly. ActinidiaBase v1 closed at 47,228 source
records, retained 47,005 chromosome-primary records, audited 223 full-source
records as excluded nonprimary, and all 47,005 remapped sequences passed the
separate gate exactly. These counts are evidence for those frozen inputs only;
the published output manifests and checksums remain the authoritative record.

The wrapper uses one process and no worker pool. It stages a complete output
directory beside the requested destination, rechecks that all three inputs
were unchanged, and atomically renames the directory only after every gate
passes. Existing output paths are never overwritten.

## Output contract

For sample `sample`, a PASS directory contains:

- `sample.publisher_primary.remapped.protein.faa`: publisher sequences with
  selected transcript IDs as first-token identifiers;
- `sample.publisher_primary.mapping.tsv`: one selected transcript per row,
  both publisher accessions, source/output sequence hashes, and PASS status;
- `sample.publisher_protein.source_inventory.tsv`: every source publisher
  protein labelled `SELECTED_AND_REMAPPED` or `EXCLUDED_NONPRIMARY`;
- `sample.publisher_primary.summary.tsv`: graph mode, source
  gene/transcript/model counts, source/selected/output counts, and the exact
  closure, one-to-one, and sequence-preservation gates;
- `input_checksums.tsv`: path-free input roles, basenames, sizes, and SHA-256;
- `run_manifest.json`: path-free schema, policies, counts, and execution
  declaration;
- `checksums.tsv`: byte count and SHA-256 for every published output except
  itself.

Accept the remap only when the command exits zero, `status` and
`publication_gate` are `PASS`, both closure columns, `one_to_one_mapping`, and
`sequence_preservation` are `true`, output count equals selected count, and
all checksums verify.

## Required next gate

The remap proves identifier provenance and sequence preservation only. It does
not assert that publisher and genome-derived proteins agree. Feed the remapped
FASTA directly to the separate compatibility gate:

```bash
python scripts/qc/audit_published_protein_compatibility.py \
  --derived-proteins primary_annotation_strict/sample.protein.faa \
  --publisher-proteins publisher_primary_remap/sample.publisher_primary.remapped.protein.faa \
  --sample-id sample \
  --output-dir publisher_protein_compatibility
```

If either stage fails, review the exact release, chromosome scope, GFF3 graph,
and publisher header schema. Do not silently drop the failed record and do not
replace an explicit mapping with a similarity search.
