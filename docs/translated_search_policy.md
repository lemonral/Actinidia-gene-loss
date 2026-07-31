# Translated genome-search policy

The final historical result files were empirically reconciled to a `tblastx`
post-filter of percent identity at least 50, bit score at least 50, and E-value
below `1e-5`, with no alignment-length minimum. The earlier 30% identity and
100-aa draft rule did not reproduce those final files.

The rebuilt analysis retains that rule as a labelled historical reproduction.
Its primary classification is more conservative. An unobserved reference gene
is callable only when left and right SynOrths anchors occur within 20 reference
genes, map unambiguously to the same target chromosome, delimit no more than
5 Mb after 10-kb padding, and the interval is at least 80% A/C/G/T. A callable
locus with no qualifying local translated hit is a positive deletion. A local
translated hit is reported as uncertain genomic sequence, not as a pseudogene,
unless an independent disruptive-mutation test supports pseudogenization.
