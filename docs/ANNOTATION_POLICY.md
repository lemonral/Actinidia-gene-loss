# Annotation policy for cross-genome comparisons

## Decision

The primary rebuild does **not** uniformly predict every genome de novo. It
uses the exact publisher genome and matched GFF/protein release, restricts the
annotation to the reviewed chromosome scope, and normalizes that annotation to
one validated coding transcript per gene. This is annotation standardization,
not reannotation.

Blanket reannotation would require the same repeat library, RNA-seq and protein
evidence, gene predictor versions, training data, filtering, and manual curation
for every species. Without equivalent evidence, a new uniform-looking pipeline
can introduce species-specific false absences and change gene-family sizes more
than the biological signal. It would also be substantially more expensive.

## Required gates for every assembly unit

1. Freeze one matched genome--GFF--publisher-protein release and its checksums.
2. Audit full versus chromosome scope; retain excluded scaffold/unplaced
   records and their genes in an exclusion table.
3. Materialize the chromosome FASTA and GFF without changing feature chains or
   coordinates.
4. Select one primary coding transcript per gene only after validating every
   candidate model.
5. Independently re-extract CDS and proteins with gffread from the generated
   selected-primary gene/mRNA/CDS GFF3 and require exact selected-ID and
   sequence agreement. The complete publisher GFF3 remains the audited Python
   selection input but is not passed to this selected-model comparison.
6. Map publisher protein identifiers to selected transcript identifiers only
   through explicit one-to-one GFF attributes, then require publisher-protein
   compatibility.
7. Run BUSCO separately on the chromosome genome and on the standardized
   derived primary protein set. Retain publisher-protein BUSCO as a separate
   source-annotation metric.
8. Record gene/model omissions, internal stops, unsupported phases, ID mapping
   exceptions, and chromosome-scope exclusions instead of silently dropping
   them.

The standardized derived protein set is the only input used for OrthoFinder,
JCVI protein matching, SynOrths, phylogeny, and comparable gene-family counts.

## When targeted reannotation is justified

An assembly is considered for targeted repair only when all of the following
hold:

- assembly identity and chromosome scope pass;
- genome BUSCO and sequence integrity are good;
- the matched annotation or derived proteins fail materially;
- no better matched publisher annotation exists;
- the failure can plausibly be corrected with adequate species-specific
  transcript/protein evidence.

Such a repair receives a new assembly-unit identifier, a separate provenance
record, and all QC, chromosome-homology, orthology, JCVI, and gene-loss gates
are rerun. It is never patched silently into a published annotation.

## Reporting language

Methods and supplementary tables should state that matched published
annotations were standardized to one validated primary isoform per gene. They
must not claim that all genomes were uniformly reannotated. Genome BUSCO,
publisher-protein BUSCO, derived-primary-protein BUSCO, gene count, transcript
count, chromosome-scope fraction, excluded unplaced/scaffold counts, and every
analysis decision are reported as distinct fields.
