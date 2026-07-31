from __future__ import annotations

from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from scripts.phylogeny.audit_tree_only_archive import audit_archive


def test_reports_tree_only_taxon_without_assembled_locus(tmp_path: Path) -> None:
    archive = tmp_path / "tree_only.zip"
    concat_member = "bundle/alignments/CONCAT.fasta"
    assembled_prefix = "bundle/assembled/"
    with ZipFile(archive, "w", ZIP_DEFLATED) as bundle:
        bundle.writestr(
            concat_member,
            ">Other_taxon\nACGTACGT\n>Target_taxon\n----ACN?\n",
        )
        bundle.writestr(
            f"{assembled_prefix}locus1.fasta",
            ">Other_taxon\nACGT\n",
        )
        bundle.writestr("bundle/published.tre", "(Target_taxon,Other_taxon);\n")

    report = audit_archive(
        archive=archive,
        taxon_label="Target_taxon",
        assembled_prefix=assembled_prefix,
        concat_member=concat_member,
    )

    assert report["assembled_exact_header_match_count"] == 0
    assert report["matching_members"] == [
        concat_member,
        "bundle/published.tre",
    ]
    assert report["concatenated_alignment"] == {
        "member": concat_member,
        "record_found": True,
        "alignment_length": 8,
        "nonmissing_sites": 2,
        "nonmissing_fraction": 0.25,
        "first_nonmissing_1based": 5,
        "last_nonmissing_1based": 6,
        "nonmissing_blocks_1based_inclusive": [[5, 6]],
        "missing_characters": "-.?NX",
    }
