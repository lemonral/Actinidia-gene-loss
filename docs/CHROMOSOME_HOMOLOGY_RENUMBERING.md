# Chromosome homology and final numbering

## Why publisher labels are provisional

`PubChr01` through `PubChr29` mean only that a record belongs to the reviewed
publisher chromosome scope. Equal publisher numbers across assemblies are not
assumed to be homologous. Every haplotype or subgenome is mapped independently;
an A--F or HAP1/HAP2 mapping is never copied from one unit to another.

The historical coordinate reference was the Hongyang v4.0 HY4A haplome, not a
generic `GWHDODM` A/B pair. Its exact matched genome, GFF, CDS, and protein
checksums are frozen in `config/chromosome_coordinate_references.tsv`. HY4P is
an independent confirmation reference. The old PAF files are retained as
cross-check evidence only: the old relabelling scripts did not validate
bijections, overlapping PAF intervals, GFF closure, or failure states.

The recovered manuscript-era script ran `minimap2 -x asm5 -t 10` with
*A. chinensis* as the minimap2 target and each other *Actinidia* assembly as
the query. A later script summed PAF matching bases for every query/target
chromosome pair. This correctly identifies the historical intent, but the
original materialization step confused PAF query and target identifiers: all
archived `reordered.fasta` products are zero-byte files. It also did not update
GFF coordinates or chromosome orientation. Therefore the PAF files and Excel
summary are historical evidence, not production-renamed genomes.

## Evidence workflow

For every reviewed 29-chromosome assembly unit, run four whole-chromosome
alignments: target to HY4A, HY4A to target, target to HY4P, and HY4P to target.
The frozen minimap2 version is `2.28-r1209`, and the exact path-independent argv
template is
`minimap2 -x asm5 --secondary=no -c --cs=long {reference_fasta} {query_fasta}`.
The target/reference FASTA is the first positional argument and the sequence
being queried is the second. Instantiate that template in each direction
without appending, dropping, or reordering argv. Every retained PAF alignment
must carry `tp:A:P`, `de:f`, `cg:Z`, and `cs:Z`; retain only those primary rows
with MAPQ at least 20, alignment block length at least 10 kb, and minimap2's
gap-compressed divergence tag `de` at most 0.15.
Coverage is the union of query and target intervals; overlapping PAF rows are
never summed as independent bases.

Build a nucleotide score matrix containing bidirectional covered bases and
fractions, matching bases, weighted divergence, orientation, row and column
reciprocal-best flags, and top-versus-second evidence. Independently build a
JCVI gene-anchor matrix against the same HY4A/HY4P references. The JCVI score
uses the harmonic mean of unique anchored-gene coverage in the two directions,
not raw anchor-row counts. Its frozen generation policy is LAST protein
alignment, `cscore=0.7`, tandem `Nmax=10`, maximum gene distance `20`, minimum
anchor block size `4`, and coverage calculated from raw JCVI anchors.

The assignment program consumes the matrices; it never launches an alignment.
Each matrix is long-form and must contain the complete Cartesian product of
exactly 29 query IDs and 29 reference IDs (841 data rows). Four separate files
are required: nucleotide/HY4A, JCVI/HY4A, nucleotide/HY4P, and JCVI/HY4P.
The exact column order is part of matrix schema `1.0.0`. Each matrix also has a
different required provenance sidecar and checksum-bound upstream `PASS`
report; a matrix without both documents is not an assignment input.

Nucleotide columns are:

```text
query_chromosome
reference_chromosome
canonical_chromosome
score
query_covered_bp
query_length_bp
query_coverage
reference_covered_bp
reference_length_bp
reference_coverage
reciprocal_coverage
matching_bases
weighted_divergence
orientation
```

`query_coverage` and `reference_coverage` must reproduce their interval-union
covered-bp fractions. `reciprocal_coverage` is their minimum. The assignment
score is
`harmonic_mean(query_coverage, reference_coverage) * (1 - weighted_divergence)`.
The frozen arithmetic tolerance is `1e-9`; matrix producers should therefore
write at least 12 significant digits rather than publication-rounded values.
The accepted orientation values are `+`, `-`, `mixed`, and `none`. Orientation
is reported as a diagnostic only. Under the author-approved production rule,
publisher sequence direction is preserved for every chromosome: no FASTA
record is reverse-complemented and no GFF coordinate or strand is flipped.
The earlier direction-harmonization implementation is retained only as
superseded diagnostic code and is not invoked by the production relabelling
queue.
An unaligned cell remains present with both covered-bp values and
`matching_bases` equal to zero, both coverage fractions and `score` equal to
zero, `weighted_divergence` equal to one, and orientation `none`.

JCVI columns are:

```text
query_chromosome
reference_chromosome
canonical_chromosome
score
query_anchored_genes
query_eligible_genes
query_gene_coverage
reference_anchored_genes
reference_eligible_genes
reference_gene_coverage
unique_anchor_pairs
```

The two gene-coverage fractions must reproduce the unique anchored-gene
numerators and eligible-gene denominators. `score` is their harmonic mean.
Each reference ID maps to exactly one canonical label and each matrix must map
bijectively onto `Chr01` through `Chr29`. HY4A and HY4P may use reference IDs
such as `Chr01A` and `Chr01P`; `canonical_chromosome` supplies their reviewed
common label.

Solve the nucleotide and JCVI matrices separately with a global one-to-one
assignment. Production `Chr01` through `Chr29` names use the HY4A nucleotide
global optimum and require only complete 29-by-29 scope plus a unique
bijection. Absolute support, HY4P, and JCVI are reported as confidence/QC and
do not block a label. The production naming policy is frozen in
`[chromosome_naming]` in `config/analysis_parameters.toml`.

The older strict four-matrix diagnostic evaluates all of the following:

- exact 29 by 29 scope and a bijection;
- nucleotide and JCVI assignments agree;
- every assigned pair is row- and column-reciprocal best;
- top/second ratio, normalized margin, and reciprocal nucleotide coverage pass;
- each assigned JCVI edge has at least 30 unique anchor pairs, reciprocal
  gene coverage at least 0.05, and assigned score at least 0.05;
- HY4A and HY4P imply the same canonical chromosome label;
- every sidecar, upstream report, matrix, target asset registry, reference asset
  registry, and reference chromosome-map registry checksum reconciles.

The strict diagnostic reciprocal nucleotide-coverage floor is 0.05. This
release-specific floor was frozen after exact primary-alignment interval-union auditing showed
0.081--0.193 coverage on the 29 assigned *A. rufa* chromosomes. It is accepted
only together with 29/29 agreement between HY4A and HY4P and 29/29 independent
agreement between nucleotide and JCVI assignments. These values classify
confidence in the revised production workflow; they do not determine whether
a unique maximum-similarity label is emitted.

For every assigned edge, the program evaluates both its query row and its
reference column. The edge must be the unique best score in both directions.
The top/second ratio and `(top - second) / top` normalized margin must pass in
both directions. The frozen ratio is 1.5 and the frozen normalized margin is
0.10. Both values remain explicit audit gates even though the ratio threshold
already implies a larger margin for positive finite scores. A zero second
score does not waive the margin, unique-best, or positive-top requirements.
All four matrices are solved independently; neither one evidence type nor HY4A
can seed another assignment.

### Fail-closed provenance contract

The target asset registry binds one `assembly_unit_id` and one
`target_scope_id` to the exact verified chromosome-scope genome, GFF, and
protein bytes. `config/chromosome_coordinate_references.tsv` binds the exact
HY4A and HY4P genome, GFF, protein, and CDS assets and their biological roles.
`config/chromosome_reference_maps.tsv` independently freezes every HY4A/HY4P
reference chromosome ID and its canonical `Chr01`--`Chr29` label. Matrix files
cannot redefine that mapping. The exact SHA-256 identities of both reference
registries are pinned in `[chromosome_homology]`; supplying a rewritten
registry under another path therefore fails before any matrix is parsed.

Each matrix sidecar must bind all of the following by basename, byte count, and
SHA-256 where applicable:

- its exact matrix role (`nucleotide_hy4a`, `jcvi_hy4a`,
  `nucleotide_hy4p`, or `jcvi_hy4p`), reference slot, reference ID, reference
  role, reference-map ID, target unit, and target scope;
- the matrix bytes, target registry and target genome/GFF/protein bindings;
- the reference-asset registry and exact HY4A or HY4P asset bindings;
- the frozen chromosome-map registry and the matrix-generation parameters;
- for nucleotide matrices, minimap2 `2.28-r1209` and the exact argv template
  above, in addition to the primary/MAPQ/block/`de`/coverage policy;
- one adjacent upstream validation report whose own checksum is frozen in the
  sidecar.

The upstream report must say `PASS`, repeat those identities and checksums, and
use the accepted matrix-builder workflow version `1.0.0`, and contain the exact
all-true check set for its matrix kind. Nucleotide reports
cover bidirectional PAF identity, primary-only filtering, MAPQ, block-length,
`de`, interval-union denominators, matrix arithmetic, and checksum closure.
JCVI reports cover bidirectional anchor inputs, BED/protein identity, eligible
gene denominators, unique anchors, matrix arithmetic, and checksum closure.
Copied HY4A/HY4P matrices, swapped A/P roles, a sidecar for another assembly
unit or scope, and a consistently but incorrectly relabelled four-matrix set
therefore fail before assignment.

Run the older strict diagnostic only after all four score matrices have their
own upstream checksums and validation reports:

```bash
python scripts/qc/assign_chromosome_homology.py \
  --nucleotide-hy4a "$DATA_ROOT/homology/UNIT/nucleotide_hy4a.tsv" \
  --jcvi-hy4a "$DATA_ROOT/homology/UNIT/jcvi_hy4a.tsv" \
  --nucleotide-hy4p "$DATA_ROOT/homology/UNIT/nucleotide_hy4p.tsv" \
  --jcvi-hy4p "$DATA_ROOT/homology/UNIT/jcvi_hy4p.tsv" \
  --nucleotide-hy4a-provenance "$DATA_ROOT/homology/UNIT/nucleotide_hy4a.provenance.json" \
  --jcvi-hy4a-provenance "$DATA_ROOT/homology/UNIT/jcvi_hy4a.provenance.json" \
  --nucleotide-hy4p-provenance "$DATA_ROOT/homology/UNIT/nucleotide_hy4p.provenance.json" \
  --jcvi-hy4p-provenance "$DATA_ROOT/homology/UNIT/jcvi_hy4p.provenance.json" \
  --parameters config/analysis_parameters.toml \
  --target-asset-registry "$DATA_ROOT/homology/UNIT/target_assets.tsv" \
  --reference-asset-registry config/chromosome_coordinate_references.tsv \
  --reference-chromosome-map-registry config/chromosome_reference_maps.tsv \
  --assembly-unit-id UNIT \
  --target-scope-id UNIT.chromosome_scope_v1 \
  --trusted-repository-commit "$REVIEWED_GIT_COMMIT" \
  --output-dir "$DATA_ROOT/homology/UNIT/assignment_v1"
```

An output directory is new and immutable. Malformed schemas, missing matrix
cells, inconsistent lengths/denominators, non-finite values, arithmetic
mismatches, changes detected while an input snapshot is being captured, or
reference-map conflicts abort before any output is published. Every input is
opened once as a nonblocking, non-symlink, bounded regular file; hashing and
parsing use those same captured bytes. A later mutation of an unrelated input
pathname cannot alter the already captured result: the output manifest retains
the exact captured hash, and a deliberately changed input requires a new run.
Valid but biologically insufficient matrices produce one atomic
failure-audit directory and a non-zero CLI exit. That directory does not
contain a final chromosome map. The command obtains an exclusive publication
lock and uses a no-replace rename; it never overwrites an existing output path,
including a dangling symlink. It reads local files only: it does not contact an
analysis server, access the network, or launch scientific computation.

Production materialization uses
`scripts/qc/relabel_chromosome_bundle_by_similarity.py`. It writes a new genome
and GFF rather than editing publisher files, changes each FASTA/GFF sequence ID
by the accepted bijection, and assigns `keep` to every orientation action.
Feature coordinates, strands, IDs, and attributes otherwise remain unchanged.
Primary CDS and protein sequences must close exactly against the relabelled
genome/GFF, and the copied CDS/protein files must remain byte-identical.

## Failure states

Malformed scope, a non-bijective HY4A maximum-similarity result, or failed
genome/GFF/CDS/protein closure stops production numbering. Low absolute
support, a nucleotide/JCVI conflict, or HY4A/HY4P disagreement remains visible
in the diagnostic output and confidence flag but does not override a complete
unique HY4A label map under the approved naming policy.

Machine-readable states from the older strict diagnostic are more specific:

- `PASS_AUTO`: the complete unit passes and a 29-row final map is emitted;
- `CONFLICT_NUCLEOTIDE_JCVI`: at least one reference-specific assignment
  differs between evidence types;
- `HY4A_HY4P_DISAGREEMENT`: both evidence types agree within each reference,
  but the inferred canonical labels disagree between references;
- `AMBIGUOUS_RECIPROCAL_BEST`: an assigned edge is not the unique row and
  column best;
- `AMBIGUOUS_SEPARATION`: a row or column ratio/margin fails;
- `LOW_RECIPROCAL_COVERAGE`: an assigned nucleotide edge fails the frozen
  minimum in HY4A or HY4P;
- `LOW_JCVI_ABSOLUTE_SUPPORT`: an assigned JCVI edge has fewer than 30 unique
  anchor pairs, reciprocal gene coverage below 0.05, or score below 0.05;
- `NON_BIJECTIVE`: the combined canonical result is not one-to-one.

Multiple states are retained in `failure_states`; the summary `status` uses a
fixed priority only to provide one sortable label. A valid failure audit is
not a reviewed structural exception. Structural exceptions remain a separate
human review process and cannot be manufactured by this command.
Rows whose local evidence passes inside a strict-diagnostic failed unit are
labelled `NOT_PUBLISHED_UNIT_FAILURE`, never `PASS_AUTO`. This describes that
diagnostic bundle only, not the production maximum-similarity naming decision.

A reviewed exception must record the reason, chosen reference chromosome,
reviewer and date, and every input checksum. Any changed input invalidates it.

### Trust boundary

The trusted inputs are the reviewed Git commit named on the command line (the
policy, schemas, and assignment code), the frozen one-unit target registry, the
two checked-in reference registries, and the upstream matrix builders that issue
the checksum-bound `PASS` reports. `run_manifest.json` records the full commit ID
and the captured policy and target/reference-registry hashes. This local gate
detects accidental changes, partial tampering, copied roles, and inconsistent
substitutions. It does not claim cryptographic authenticity against an adversary
who can rewrite the reviewed code and every trust root coherently. A signed
release manifest may be added at a future publication boundary.

## Publication outputs

Each unit produces nucleotide and JCVI score matrices, a chromosome-assignment
table, validation JSON, checksums, and a final reviewed map. Only `PASS_AUTO`
or checksum-valid `PASS_REVIEWED_EXCEPTION` maps can materialize final genome
and GFF files. Final files are rebuilt from frozen publisher-scope inputs, not
edited in place, and must reproduce identical CDS and protein sequences after
renaming and any reviewed whole-chromosome reverse complement.

The assignment bundle contains normalized copies of all four input matrices,
four per-matrix assignment-evidence tables, the combined unit table, a summary,
validation JSON, a run manifest, and SHA-256 checksums. The manifest records
the validated matrix-sidecar/report and registry bindings, but contains only
input basenames, byte counts, and checksums; absolute runtime paths are never
published. `PASS_AUTO` additionally contains
`UNIT.final_chromosome_map.tsv`. A failure bundle deliberately omits it.

This chromosome-homology JCVI analysis is distinct from the bidirectional
*Clematoclethra scandens* JCVI analysis used for gene-loss inclusion and
callable coverage.
