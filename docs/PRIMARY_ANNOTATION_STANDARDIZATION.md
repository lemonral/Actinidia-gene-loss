# Primary annotation standardization

This stage derives one internally consistent CDS/protein set from each matched
chromosome-scope genome and GFF3 pair. It does not predict new genes and does
not modify or replace the publisher files. Raw, full-scope, and excluded
scaffold records remain in the external data store; the standardized bundle is
a derived analysis input.

Use `scripts/qc/extract_primary_annotation.py`, not the older convenience
extractor, for production. The wrapper uses one process, starts no worker pool,
stages all files in a sibling temporary directory, and renames that directory
into place only after every validation succeeds. It never overwrites an
existing output directory.

## Required inputs

Run this stage only after `materialize_chromosome_scope.py` has produced a
reviewed chromosome FASTA/GFF3 pair. The inputs must satisfy all of the
following conditions:

- FASTA sequence IDs are unique and every GFF3 `seqid` is present in the FASTA;
- intervals are within sequence bounds;
- transcript rows have one `ID` and exactly one gene `Parent`;
- CDS rows have one or more declared transcript parents;
- transcript and CDS sequence IDs, strands, and intervals agree;
- CDS phase is `0`, `1`, or `2` by default;
- gene and transcript IDs are globally unambiguous.

Plain and gzip-compressed FASTA/GFF3 are accepted. Compression is detected by
file content rather than filename suffix.

### Explicit gene-as-transcript compatibility mode

Some publisher GFF3 files, including candidate *A. rufa* annotations, encode a
`gene -> CDS` graph and contain no `mRNA` or `transcript` rows. The strict
default still rejects those files. After the matched genome/GFF3 release and
graph shape have been audited, add `--gene-as-transcript` to enable the narrow
compatibility mode.

The option is a fallback, not a general graph repair. It synthesizes exactly
one in-memory self-transcript per declared gene only when the complete input
contains zero accepted transcript-feature rows. The self-transcript keeps the
publisher gene ID unchanged as its transcript, CDS/protein FASTA, coordinate,
and audit-table ID. A CDS must have exactly one `Parent`, and that parent must
identify a declared gene-level row. The run rejects undeclared parents,
multiple parents, CDS attached directly to a declared `pseudogene`, and mixed
files containing both declared transcript models and direct gene-parent CDS.
Transcript-level `--canonical-tag` rules are not meaningful for synthesized
models and are therefore rejected in this mode.

A gene or pseudogene without CDS is still represented in the audits. It is
reported as `no_CDS_transcript` / `noncoding_or_no_CDS`, is not counted as an
invalid coding gene, and does not produce a protein. Thus a noncoding
pseudogene does not block a bundle that also contains valid coding genes.

Because a GFF3 ID belongs to one global namespace, a gene row and its
self-transcript cannot both use the unchanged publisher ID in the normalized
comparison GFF3. In compatibility mode, `sample.primary.gff3` therefore
contains a top-level `mRNA` with the publisher gene ID plus its selected CDS
rows, and deliberately omits the duplicate-ID gene row. gffread is run on this
minimal `mRNA/CDS` object and must return the exact preserved selected-ID set
and the same sequences as the Python extractor. The manifest records graph
mode `gene_as_transcript`, zero source transcript rows, the number of
synthesized self-transcripts, and comparison scope
`selected_primary_GFF3_top_level_mRNA_CDS_only_gene_as_transcript`.

```bash
python scripts/qc/extract_primary_annotation.py \
  --genome chromosome_scope.genome.fa.gz \
  --gff chromosome_scope.annotation.gff3.gz \
  --sample-id act_rufa_aru_r1 \
  --output-dir primary_annotation \
  --gene-as-transcript \
  --gffread /path/to/gffread \
  --require-gffread
```

## Transcript validation and selection

Every transcript is validated before any isoform is selected. A coding
candidate is valid only when its strand-aware CDS can be reconstructed, the
GFF3 phase chain is consistent, the CDS length is divisible by three, the
translation is non-empty, and no internal stop codon is present. A terminal
stop codon is retained in the CDS and removed from the protein. Ambiguous
codons are translated as `X` and reported in `QC_flags`; absence of an initial
methionine or terminal stop is also reported rather than used to invent a gene
model.

For each gene, selection is deterministic:

1. If ordered `--canonical-tag` rules were configured, use valid transcripts
   matching the first rule that has a valid match for that gene.
2. Among the remaining eligible transcripts, choose the longest validated
   spliced CDS.
3. Break a length tie by the longer genomic transcript span.
4. Break a second tie by lexical transcript ID.

A rule can test presence (`--canonical-tag canonical`) or an exact/comma-token
value (`--canonical-tag tag=canonical`). Repeat the option to declare a lower
priority rule. Canonical tags are preferences, not exemptions from sequence
validation: an invalid tagged transcript cannot displace a valid untagged
transcript.

The strict default rejects a run if any gene has CDS rows but no valid coding
transcript. After inspecting and documenting those genes, an explicitly
approved sensitivity or publisher-compatibility run may use
`--invalid-coding-gene-policy omit`; every omission remains in the gene and
transcript audit tables. Likewise, missing CDS phase is rejected by default.
`--missing-phase-policy zero` is available only as an explicit, audited legacy
compatibility choice.

## Independent gffread gate

The Python implementation first validates the complete chromosome-scope GFF3,
selects one coding transcript per gene, and writes
`sample.primary.gff3`. Normally that normalized file contains only the selected
gene, mRNA, and CDS rows. The explicit gene-as-transcript mode instead contains
the selected top-level self-mRNA and CDS rows, as described above, so the
publisher ID is not changed or duplicated. When gffread is found, the wrapper
runs it on the original chromosome FASTA plus this selected-primary GFF3; it
does **not** pass the full publisher GFF3 to gffread.

This comparison object is deliberate. Unselected isoforms and non-CDS features
such as publisher `start_codon`/`stop_codon` rows cannot extend or otherwise
change a selected transcript during the independent extraction, and gffread
does not allocate/write the complete publisher annotation when only selected
models are being checked. The full source GFF3 is still the frozen input to the
Python validation and primary-selection audit; it is neither ignored nor
modified.

The gffread CDS and protein ID sets must each exactly equal the selected
transcript ID set, with no missing or extra record. CDS agreement is exact.
Protein agreement is exact after removal of a terminal `*` only. Any ID-set or
sequence difference rejects the complete run and publishes nothing. The
manifest records the selected-primary GFF3 basename, byte count, SHA-256,
scope `selected_primary_GFF3_gene_mRNA_CDS_only`, and the explicit fact that
the full source GFF3 was not passed to gffread. The selected GFF3 is also
checked for modification during the external command.

Production runs must use `--require-gffread` and record an explicit executable
when it is not already on `PATH`. A development run without gffread is allowed,
but its comparison table and manifest state
`NOT_RUN_GFFREAD_NOT_AVAILABLE`; an explicitly disabled check states
`NOT_RUN_EXPLICITLY_DISABLED`. In both cases the manifest publication gate is
`BLOCKED_GFFREAD_NOT_RUN`.

The current production environment freezes official gffread 0.12.9. Its
Linux x86-64 release archive SHA-256 is
`3effcbb3ccf5a76305887df0d39feb69160aed6091e42f173bb0ebf422f18d17`.
The executable is stored outside Git with the other versioned tools; only the
version and checksum belong in the repository.

```bash
python scripts/qc/extract_primary_annotation.py \
  --genome "$DATA_ROOT/standardized/act_eriantha_hap1/chromosome_scope.genome.fa.gz" \
  --gff "$DATA_ROOT/standardized/act_eriantha_hap1/chromosome_scope.annotation.gff3.gz" \
  --sample-id act_eriantha_hap1_2026 \
  --output-dir "$DATA_ROOT/standardized/act_eriantha_hap1/primary_annotation" \
  --gffread /path/to/gffread \
  --require-gffread
```

## Publisher-protein compatibility gate

Independent gffread agreement proves that the Python and gffread extractors
interpret the selected CDS models against the declared genome consistently;
it does not prove that the full annotation belongs to the same release as the
publisher protein FASTA. After the primary bundle passes gffread, run
`scripts/qc/audit_published_protein_compatibility.py` in a separate atomic
output directory. Production requires exact first-token protein ID-set
equality and residue agreement after only terminal-stop normalization and an
explicitly audited, one-way publisher-`X` wildcard against a canonical derived
amino acid. Every other difference fails the complete compatibility run.

See `docs/PUBLISHED_PROTEIN_COMPATIBILITY.md` for the exact sequence policy,
command, output schemas, and failure interpretation. A failed compatibility
gate must trigger a release/scope/asset-pair review; it must not be bypassed by
silently omitting affected genes.

Do not add a `--canonical-tag` merely because another assembly uses one. Rules
are assembly-specific input metadata and must be frozen in the analysis-unit
manifest before production.

## Output contract

For sample ID `sample`, the atomic directory contains:

- `sample.protein.faa`: one terminal-stop-free protein per selected transcript;
- `sample.cds.fa`: the matching strand- and phase-aware CDS;
- `sample.primary.gff3`: a minimal normalized gene/mRNA/CDS annotation, or the
  explicitly labelled top-level self-mRNA/CDS comparison object in
  gene-as-transcript mode;
- `sample.primary_isoforms.tsv`: gene-to-selected-transcript map, lengths,
  selection rule, and QC flags;
- `sample.transcript_audit.tsv`: every source transcript, validation result,
  canonical match, and disposition;
- `sample.gene_audit.tsv`: every source gene ID, valid/invalid transcript
  counts, selection, or omission reason;
- `sample.coords.tsv`: headered coordinates for downstream joins;
- `sample.coords`: headerless five-column SynOrths compatibility file;
- `sample.gffread_comparison.tsv`: two PASS rows per selected transcript in a
  publication-ready run; `run_manifest.json` binds that comparison to the
  checksum of `sample.primary.gff3` and records exact selected-ID closure;
- `run_manifest.json`: path-free input checksums, policy, counts, tool status,
  and one-process execution declaration;
- `checksums.tsv`: checksums of all published files except itself.

FASTA IDs are transcript IDs. Downstream gene-level joins must use
`sample.primary_isoforms.tsv`; they must not infer a gene ID by truncating a
transcript name.

## Acceptance checklist

A bundle is accepted only when:

- the command exits zero and `run_manifest.json` has `status: PASS`;
- `publication_gate` and gffread comparison status are both `PASS` for production;
- the separate published-protein compatibility directory has
  `publication_gate: PASS`, exact ID-set equality, and zero nonpermitted
  sequence mismatches;
- `selected_genes`, protein records, CDS records, primary-isoform rows, and
  coordinate rows are equal;
- the invalid coding-gene count agrees with the declared policy;
- all files listed in `checksums.tsv` verify;
- BUSCO is rerun on this derived primary protein set, while the publisher
  protein BUSCO remains separately reported as annotation-source QC.

This standardized protein set is used by OrthoFinder, phylogeny, gene-family,
JCVI protein matching, and SynOrths. It is not evidence that all assemblies
were uniformly reannotated; the Methods and supplementary table must state
that matched published annotations were normalized to one validated primary
isoform per gene.
