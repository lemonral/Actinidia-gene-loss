# Audit of the manuscript-era chromosome numbering

The recovered legacy workflow used `Actinidia_chinensis.fasta` as the
minimap2 target and each other *Actinidia* assembly as the query:

```text
minimap2 -x asm5 -t 10 Actinidia_chinensis.fasta TARGET_SPECIES.fasta
```

Its second-stage summary grouped PAF rows by query chromosome and Chinese
kiwifruit target chromosome, summed matching bases, and retained the largest
sum. This documents the original biological intention: transfer chromosome
numbers from *A. chinensis* by whole-genome similarity.

The archived implementation is not used to publish renamed genomes. In PAF,
column 1 is the assembly supplied as the minimap2 query and column 6 is the
minimap2 target. The old materializer treated these roles inconsistently when
looking up target records, and every archived `reordered.fasta` output is
zero bytes. It also selected local evidence without a global one-to-one gate,
did not reconcile overlapping PAF intervals, did not confirm with an
independent *A. chinensis* haplome or JCVI anchors, did not update GFF, and did
not harmonize chromosome direction.

The rebuild preserves the old PAF and Excel files as historical evidence but
repeats chromosome assignment for every new or changed 29-chromosome unit. It
uses exact Hongyang v4.0 HY4A assets as the coordinate and orientation
reference, HY4P as an independent label confirmation, bidirectional minimap2
and JCVI evidence, and a global one-to-one assignment. A separate verified
materializer then renames chromosomes and, where the HY4A orientation evidence
is strong and unambiguous, reverse-complements the whole chromosome while
transforming every GFF coordinate and strand. CDS and protein sequences must
remain exactly identical after that transformation.
