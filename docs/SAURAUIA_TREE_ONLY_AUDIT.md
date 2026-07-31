# *Saurauia tristyla* tree-only archive audit

## Result

The exact `Angiosperms353.zip` file from Zenodo record 13835677 is useful as
published topology and age context, but it is **not a compatible source of
Angiosperms353 loci for the new nuclear species tree**.

The checksum-verified archive was audited for the exact label
`Actinidiaceae_Saurauia_tristyla`. The label occurs in five members:

- the published concatenated alignment;
- the published dated tree;
- the published molecular-branch-length tree;
- the ITS tree; and
- the published ASTRAL tree.

It occurs in zero files under `Angiosperms353/assembled_seqs/`. In the 7,698 bp
published concatenated alignment, the terminal has 685 non-missing sites
(8.90%), all between positions 6,874 and 7,629. This terminal is therefore
supported by the appended ITS region rather than the archive's 353 nuclear
target loci.

Consequences are fail-closed:

1. do not add this archive to OrthoFinder, the strict single-copy alignment,
   ASTRAL-Pro, CAFE, or a gene-loss denominator;
2. do not graft the published *S. tristyla* tip onto a newly inferred nuclear
   tree;
3. keep the *Parasaurauia* calibration disabled until a compatible
   *Saurauia* nuclear-locus or genome/proteome asset is obtained and passes the
   same identifier, occupancy, and topology gates; and
4. retain the archive only as external published topology/age context.

NCBI currently indexes no nuclear genome assembly for this species, but it
does index two paired-end RNA-seq runs (`SRR11994221` and `SRR28027655`). They
are recorded as not-yet-downloaded rescue candidates in
`config/phylogeny/saurauia_nuclear_rescue_candidates.tsv`. If the
Actinidiaceae fossil is retained, the economical rescue is to recover only the
same frozen nuclear target-locus set from one run, validate orthology and
occupancy, and use the second biological sample as a sensitivity. A full de
novo transcriptome is not required merely to create a tree tip.

## Reproduction

Run the audit against the private bulk-data root:

```bash
python scripts/phylogeny/audit_tree_only_archive.py \
  --archive "$GENELOSS_DATA/downloads/phylogeny_tree_only/Saurauia_tristyla/Angiosperms353.zip" \
  --taxon-label Actinidiaceae_Saurauia_tristyla \
  --assembled-prefix Angiosperms353/assembled_seqs/ \
  --concat-member Angiosperms353/alignments/molecular_branch_length_estimation_alignment/CONCAT.fasta \
  --output "$GENELOSS_DATA/qc/phylogeny/saurauia_tristyla_tree_only_audit.json"
```

The expected archive identity is 37,809,948 bytes, MD5
`0f667a1877774bccc639ae204b6b8f33`, and SHA-256
`27f2dd659d5b46411fddc2646865686eef324d86082da0bc3783c397fe5d7ae3`.
