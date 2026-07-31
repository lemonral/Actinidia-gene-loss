# Reviewed chromosome-scope maps

These tables are explicit publisher-scope mappings used to create provisional
chromosome-only genome/GFF3 pairs. They are not, by themselves, evidence that
the publisher chromosome number is homologous to the same-numbered
`A. chinensis` chromosome.

The `PubChr01`--`PubChr29` labels deliberately preserve that distinction. A
second reviewed mapping based on bidirectional chromosome similarity and JCVI
anchors will assign final `Chr01`--`Chr29` analysis labels. The final mapping
must be one-to-one or record a documented exception; sequence order is never
used as evidence.

Mapping evidence for the files currently present here is the chromosome label
embedded in the publisher FASTA header, reconciled to the matched GFF3
sequence ID. The HAP1/HAP2 maps apply only to the exact six-file Figshare v1
bundle; the quarantined NCBI-FASTA/Figshare-GFF diagnostic is not a valid input.
For Fuchu, only the 29 records explicitly labelled chromosome
1--29 are retained. Records labelled `unknown` are excluded but remain in the
materializer audit. The exact ActinidiaBase v1 download contains 38 records,
not 100 contigs: `Chr1`--`Chr29` plus nine extra contigs. Its explicit 29-row
map retains the 29 publisher pseudochromosomes; the nine extras remain in the
materializer audit. One extra record, `Contig01298`, has no GFF feature row.
See `docs/RUFA_BUNDLE_IDENTITY.md` for the frozen scope and sequence-identity
evidence.

`catharanthus_roseus_asm2450571.publisher_scope.tsv` is an outgroup scope map,
not an *Actinidia* homology map. Its eight linkage-group rows come from the
official NCBI assembly report for `GCA_024505715.1`. The 313 unplaced primary
scaffolds and three non-nuclear chloroplast scaffolds remain in the source
bundle and the materialization audit, but are excluded from the
chromosome-scope phylogeny input. Publisher linkage-group labels are retained;
they are not renamed against the 29-chromosome *A. chinensis* registry.
