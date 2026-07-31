"""Tests for metadata-driven Matplotlib taxon labels."""

from __future__ import annotations

import unittest

from geneloss_repro.labels import (
    TaxonLabelError,
    format_downstream_taxon_label,
    format_taxon_label,
    format_taxon_label_from_metadata,
)


class TaxonLabelTest(unittest.TestCase):
    def test_all_deliciosa_and_zhejiangensis_subgenomes_are_metadata_driven(self) -> None:
        units = [
            *(('Actinidia deliciosa', suffix) for suffix in "ABCDEF"),
            *(('Actinidia zhejiangensis', suffix) for suffix in "AB"),
        ]
        labels = [
            format_taxon_label_from_metadata(
                {
                    "biological_species": species,
                    "haplotype_or_subgenome": suffix,
                },
                suffix_fields=("haplotype_or_subgenome",),
            )
            for species, suffix in units
        ]

        self.assertEqual(labels[0], r"$\mathit{Actinidia\ deliciosa}$ $\mathrm{A}$")
        self.assertEqual(labels[5], r"$\mathit{Actinidia\ deliciosa}$ $\mathrm{F}$")
        self.assertEqual(labels[6], r"$\mathit{Actinidia\ zhejiangensis}$ $\mathrm{A}$")
        self.assertEqual(labels[7], r"$\mathit{Actinidia\ zhejiangensis}$ $\mathrm{B}$")

    def test_haplotypes_accession_and_scope_are_explicitly_upright(self) -> None:
        label = format_taxon_label_from_metadata(
            {
                "biological_species": "Actinidia eriantha",
                "haplotype_or_subgenome": "HAP1",
                "accession": "GCA_030345205.1",
                "assembly_scope": "chromosome_partition",
            }
        )
        self.assertEqual(
            label,
            r"$\mathit{Actinidia\ eriantha}$ $\mathrm{HAP1}$ | "
            r"$\mathrm{GCA\_030345205.1}$ | $\mathrm{chromosome\_partition}$",
        )
        self.assertNotIn(r"\mathit{HAP1", label)

        hap2 = format_taxon_label(
            "Actinidia eriantha",
            ["HAP2"],
            abbreviate_genus=True,
        )
        self.assertEqual(hap2, r"$\mathit{A.\ eriantha}$ $\mathrm{HAP2}$")

    def test_unsuffixed_species_omits_empty_metadata(self) -> None:
        label = format_taxon_label_from_metadata(
            {
                "biological_species": "Actinidia rufa",
                "haplotype_or_subgenome": "",
            }
        )
        self.assertEqual(label, r"$\mathit{Actinidia\ rufa}$")

    def test_author_approved_downstream_labels_are_concise(self) -> None:
        self.assertEqual(
            format_downstream_taxon_label(
                "Actinidia x zhejiangensis parental lineage A",
                abbreviate_genus=True,
            ),
            r"$\mathit{A.\ zhejiangensis}$ $\mathrm{A}$",
        )
        self.assertEqual(
            format_downstream_taxon_label(
                "Actinidia x zhejiangensis parental lineage B",
                abbreviate_genus=True,
            ),
            r"$\mathit{A.\ zhejiangensis}$ $\mathrm{B}$",
        )
        self.assertEqual(
            format_downstream_taxon_label(
                "Actinidia rufa", ["ActinidiaBase v1"], abbreviate_genus=True
            ),
            r"$\mathit{A.\ rufa}$",
        )
        self.assertEqual(
            format_downstream_taxon_label(
                "Actinidia macrosperma", abbreviate_genus=True
            ),
            r"$\mathit{A.\ macrosperma}$",
        )

    def test_hybrid_marker_is_omitted_and_lineage_suffix_is_upright(self) -> None:
        label = format_taxon_label(
            "Actinidia x zhejiangensis", ["A"], abbreviate_genus=True
        )
        self.assertEqual(
            label,
            r"$\mathit{A.\ zhejiangensis}$ $\mathrm{A}$",
        )
        self.assertNotIn(r"\mathit{A}$", label)

        parental_lineage = format_taxon_label(
            "Actinidia x zhejiangensis parental lineage A",
            ["A"],
            abbreviate_genus=True,
        )
        self.assertEqual(
            parental_lineage,
            r"$\mathit{A.\ zhejiangensis}$ "
            r"$\mathrm{parental\ lineage\ A}$",
        )
        self.assertEqual(parental_lineage.count(r"\mathrm{A}"), 0)

    def test_unphased_suffix_is_omitted_from_publication_labels(self) -> None:
        self.assertEqual(
            format_taxon_label(
                "Actinidia rufa", ["unphased"], abbreviate_genus=True
            ),
            r"$\mathit{A.\ rufa}$",
        )
        self.assertEqual(
            format_taxon_label(
                "Actinidia macrosperma",
                ["unresolved polyploid unit"],
                abbreviate_genus=True,
            ),
            r"$\mathit{A.\ macrosperma}$",
        )

    def test_suffix_columns_and_order_are_selected_by_metadata(self) -> None:
        row = {
            "species_name": "Clematoclethra scandens",
            "assembly_unit": "reference",
            "accession": "ASM1",
        }
        label = format_taxon_label_from_metadata(
            row,
            species_field="species_name",
            suffix_fields=("accession", "assembly_unit"),
            separator=", ",
        )
        self.assertEqual(
            label,
            r"$\mathit{Clematoclethra\ scandens}$ $\mathrm{ASM1}$, $\mathrm{reference}$",
        )

    def test_malformed_biological_species_fails_clearly(self) -> None:
        malformed = [
            "Actinidia",
            "actinidia eriantha",
            "Actinidia Eriantha",
            "Actinidia eriantha HAP1",
            "Actinidia ×",
            "A. eriantha",
            "",
        ]
        for value in malformed:
            with self.subTest(value=value):
                with self.assertRaisesRegex(TaxonLabelError, "full two-word Latin binomial"):
                    format_taxon_label(value)

        with self.assertRaisesRegex(TaxonLabelError, "must be a string"):
            format_taxon_label(None)  # type: ignore[arg-type]

    def test_missing_species_and_unsafe_suffix_fail_clearly(self) -> None:
        with self.assertRaisesRegex(TaxonLabelError, "missing required species field"):
            format_taxon_label_from_metadata({"haplotype_or_subgenome": "A"})
        with self.assertRaisesRegex(TaxonLabelError, "unsafe"):
            format_taxon_label("Actinidia deliciosa", ["A$B"])
        with self.assertRaisesRegex(TaxonLabelError, "iterable of values"):
            format_taxon_label("Actinidia deliciosa", "HAP1")
        with self.assertRaisesRegex(TaxonLabelError, "sequence of column names"):
            format_taxon_label_from_metadata(
                {"biological_species": "Actinidia deliciosa"},
                suffix_fields="accession",  # type: ignore[arg-type]
            )


if __name__ == "__main__":
    unittest.main()
