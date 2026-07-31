# Reviewer response: assembly and annotation quality

## Reviewer comment (Line 131)

> The study uses multiple genome assemblies for comparative analyses (e.g., gene copy number, orthologues, paralogues, and functional classification). A clear assessment of the quality of these assemblies is necessary. Please provide a summary in the supplementary material (table or figure), including metrics such as genome size, BUSCO completeness (for both genome and gene sets), and number of annotated genes.

## Response

We agree and have added a uniform assembly and annotation quality summary for all 23 *Actinidia* assembly units used in the revised chromosome-based analyses (Supplementary Table S27 and the accompanying QC figure). The table reports assembly identity and analysis scope, genome size, sequence count, N content, sequence N50/L50, annotated gene and transcript counts, the protein set used, and BUSCO completeness for both the genome and the analyzed protein set.

All comparable BUSCO values were generated or recovered from exact input-bound runs using BUSCO v5.8.2 and `embryophyta_odb10` (dataset date 2024-01-08; 1,614 orthologues). Previously completed *Actinidia* BUSCO results were reused only when the current genome resolved to the same source file and the current protein FASTA had the same SHA-256 checksum. BUSCO was rerun for newly downloaded or changed *Actinidia* assemblies. Outgroup BUSCO values were not mixed into this table.

Assembly-wide QC and chromosome-dependent analysis scope are now reported separately. For example, the selected *A. rufa* release contains 38 nuclear sequences (29 pseudomolecules plus nine unplaced contigs) and was assessed as the complete release. The nine unplaced contigs remain documented in the QC audit, whereas JCVI, chromosome-position, and gene-loss analyses use the explicitly declared 29-pseudomolecule subset. The publisher annotation contains 47,228 genes; the chromosome analysis contains 47,005 valid primary coding genes, with six invalid coding models omitted under the frozen annotation criteria. Thus, unplaced sequences were not silently discarded from assembly QC and were not assigned artificial chromosome positions.

The revised dataset also treats the six *A. deliciosa* units (A-F), four *A. arguta* units (A-D), two *A. eriantha* haplotypes (HAP1/HAP2), and two *A. × zhejiangensis* parental-lineage units (A/B) separately where unit-level gene-loss evidence is required. *A. eriantha* HAP1 and HAP2 are paired haplotypes from one individual and are not described as independent species replicates. Primary coding sets for newly processed assemblies use the longest valid spliced CDS per gene and require frame, internal-stop, and genome/GFF/CDS/protein sequence closure. We did not perform blanket reannotation because this would introduce an additional annotation-pipeline effect across otherwise publisher-supported gene sets.

These additions make assembly quality, annotation scope, and the distinction between full-release QC and chromosome-only comparative analyses explicit and auditable.

The downstream revision uses these same 23 declared chromosome analysis units
without treating haplotypes or subgenomes as independent species replicates.
Chromosome labels were harmonized to Hongyang v4 HY4A by a global one-to-one
maximum-nucleotide-similarity assignment, while publisher sequence direction
was preserved. The relabelled genome and GFF were accepted only after exact
CDS/protein sequence closure. Gene-loss denominators, position analyses,
expression/copy-number summaries, and NLR summaries therefore all refer to the
same frozen unit identities and chromosome names.

For the nine newly harmonized units, the chromosome-naming audit reports a
complete unique `Chr01`--`Chr29` bijection in every unit (261 labels total).
The primary name is the global one-to-one maximum-nucleotide-similarity match
to HY4A. Absolute nucleotide coverage and the independent HY4P/JCVI results are
shown as confidence/QC rather than used to suppress a unique name. This is the
author-approved naming rule; no chromosome was reverse-complemented merely to
match HY4A direction. The corresponding figure and plot-data table retain the
lower-support calls visibly instead of hiding them.

All revised claims are linked to their path-free validation and checksum
artifacts in the private evidence index. The index also preserves negative
publication gates: ordinary Gaussian PGLS remains exploratory pending a
denominator-aware phylogenetic count model, Gamma3 CAFE support is not claimed,
and no centromere analysis is reported without independent centromere
intervals.

## Proposed Methods text

Assembly statistics were calculated from the exact FASTA files used in the revised analysis. Genome and protein completeness were evaluated with BUSCO v5.8.2 against `embryophyta_odb10` (dataset date 2024-01-08; n = 1,614). Previously completed BUSCO results were reused only after exact genome-path and protein-checksum identity checks; newly downloaded or changed *Actinidia* assemblies were evaluated de novo. For chromosome-dependent analyses, publisher-defined pseudomolecules were retained, while unplaced sequences were preserved in assembly-wide QC and audit tables but excluded from chromosome-position calculations. For newly standardized annotations, one longest valid spliced coding transcript was selected per gene and required exact genome/GFF/CDS/protein compatibility.

Chromosome names in newly processed 29-pseudomolecule units were assigned by
the global one-to-one maximum-nucleotide-similarity mapping to Hongyang v4
HY4A, requiring one occurrence of every label from `Chr01` to `Chr29`.
Independent HY4P and JCVI assignments and absolute alignment support were
retained as QC diagnostics. Publisher sequence direction was preserved.
Genome and GFF sequence IDs were changed by the same bijection, and the
relabelled bundle was accepted only after exact CDS and protein sequence
closure. The article-comparable main trend uses one rule for all old and new
units: an exact SynOrths anchor is retained; among historical missing-gene
candidates, a genome-wide tBLASTX hit with identity at least 50%, bit score at
least 50, and e-value below `1e-5` is decayed, without a length minimum; a
candidate without such a hit is deleted. The main loss numerator is decayed
plus deleted, and rows outside the historical candidate scope are not called.
A separate conservative evidence layer records callable local deletion and
high-quality Miniprot frameshift/in-frame-stop support. This layer refines
mechanism but never rewrites or double-counts the article classes.

## Proposed supplementary legends

**Supplementary Table S27. Assembly and annotation quality of the 23
*Actinidia* chromosome-analysis units.** The table reports the declared
biological species and assembly-unit scope, genome size, sequence count, N
content, sequence N50/L50, publisher gene and transcript counts, analyzed
primary coding-gene count, and exact-bound genome/protein BUSCO summaries.
Genome-wide release scope and chromosome-analysis scope are distinguished;
unplaced sequences are documented but do not receive artificial chromosome
positions. No outgroup BUSCO values are included.

**Assembly and annotation QC figure.** Panels summarize genome size and
contiguity, analyzed coding-gene counts, and exact-bound genome/protein BUSCO
completeness for the same 23 units. Unit suffixes denote haplotypes,
subgenomes, or assembly releases and are not independent species replicates.
Missing or inapplicable values remain explicit rather than being imputed.

**NLR repertoire and classified non-shared reference-NLR loss figure.** The
left panel reports complete NLR repertoires and six mutually exclusive
loss-evidence groups for each of the 23 chromosome-analysis units. The right
panel reports `decayed + deleted` as a percentage of resolved `retained +
decayed + deleted` comparisons after excluding 138 reference NLR genes
positive in all 23 units. The remaining 76 non-shared reference NLR genes
contribute 1,738 resolved comparisons and 254 positive calls: 14
frameshift-supported, 14 in-frame-stop-supported, 13 combined-disruption, 16
no-qualifying-hit, no truncation/partial-alignment candidate, and 197
residual-sequence mechanism-unresolved calls. Not-called rows enter neither
numerator nor denominator, and no species aggregation is used. Latin
binomials are italicized, whereas haplotype, subgenome, parental-lineage, and
assembly-release suffixes are upright.

**Mechanism-stratified spatial distribution of gene-loss evidence.** Target
assembly residual coordinates are divided among six loss-evidence groups after
all 23 assemblies are harmonized to HY4A `Chr01`--`Chr29`. Chromosome panels
compare observed counts with within-unit chromosome-length opportunities.
Candidate same-chromosome and interchromosomal displacements report the best
existing genome-wide residual alignment and are not interpreted as confirmed
inversions or translocations. Unlocalized calls are not given synthetic
coordinates. The audited *C. scandens* GFF lacks repeat/transposon features;
therefore no TE-association panel is shown.
