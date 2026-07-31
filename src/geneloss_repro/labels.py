"""Metadata-driven taxon labels for Matplotlib figures.

The functions in this module return Matplotlib MathText strings.  The Latin
binomial is explicitly italic, while assembly-unit qualifiers are explicitly
upright.  This keeps labels typographically correct even when a caller changes
the surrounding tick-label font style.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from typing import Any


DEFAULT_SUFFIX_FIELDS = (
    "haplotype_or_subgenome",
    "accession",
    "assembly_scope",
)

_BINOMIAL_PATTERN = re.compile(
    r"^(?P<genus>[A-Z][a-z]+) (?:(?P<hybrid>x|×) )?"
    r"(?P<epithet>[a-z]+(?:-[a-z]+)*)"
    r"(?: parental lineage (?P<parental_lineage>[A-Za-z0-9][A-Za-z0-9._+-]*))?$"
)
_SAFE_SUFFIX_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 ._:/()+-]*$")
_OMITTED_PUBLICATION_SUFFIXES = {
    "unphased",
    "actinidiabase v1",
    "unresolved polyploid unit",
    "unresolved_polyploid_unit",
}


class TaxonLabelError(ValueError):
    """Raised when taxon metadata cannot produce an unambiguous label."""


def _parse_biological_species(
    biological_species: str,
) -> tuple[str, str, bool, str | None]:
    if not isinstance(biological_species, str):
        raise TaxonLabelError(
            "biological_species must be a string containing a full Latin binomial"
        )
    normalized = " ".join(biological_species.split())
    match = _BINOMIAL_PATTERN.fullmatch(normalized)
    if match is None:
        raise TaxonLabelError(
            "biological_species must be a full two-word Latin binomial, optionally "
            "with an upright hybrid marker x/× or an upright 'parental lineage' "
            "qualifier, and have a capitalized genus and lowercase species epithet; received "
            f"{biological_species!r}"
        )
    return (
        match.group("genus"),
        match.group("epithet"),
        match.group("hybrid") is not None,
        match.group("parental_lineage"),
    )


def _upright_mathtext(value: str | None, *, field_name: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise TaxonLabelError(
            f"{field_name} must be a string or None; received {type(value).__name__}"
        )
    normalized = " ".join(value.split())
    if not normalized:
        return None
    if _SAFE_SUFFIX_PATTERN.fullmatch(normalized) is None:
        raise TaxonLabelError(
            f"{field_name} contains characters that are unsafe in a Matplotlib "
            f"taxon label: {value!r}"
        )
    escaped = normalized.replace("_", r"\_").replace(" ", r"\ ")
    return rf"$\mathrm{{{escaped}}}$"


def format_taxon_label(
    biological_species: str,
    suffixes: Iterable[str | None] = (),
    *,
    abbreviate_genus: bool = False,
    separator: str = " | ",
) -> str:
    """Return a Matplotlib label with an italic binomial and upright suffixes.

    Parameters
    ----------
    biological_species
        A full Latin binomial, for example ``"Actinidia deliciosa"``.
    suffixes
        Optional assembly-unit qualifiers such as subgenome, haplotype,
        accession, or assembly scope. Empty values and ``None`` are omitted.
    abbreviate_genus
        Render the genus as its initial (for example, ``A. deliciosa``) while
        retaining the full binomial as the validated input.
    separator
        Plain-text separator placed between non-empty upright suffixes.

    Notes
    -----
    The returned value can be passed directly to ``Axes.set_xticklabels``,
    ``Axes.set_yticklabels``, or any Matplotlib text method with MathText
    enabled (the Matplotlib default).
    """

    genus, epithet, _hybrid, parental_lineage = _parse_biological_species(
        biological_species
    )
    displayed_genus = f"{genus[0]}." if abbreviate_genus else genus
    species_label = rf"$\mathit{{{displayed_genus}\ {epithet}}}$"

    if isinstance(suffixes, (str, bytes)):
        raise TaxonLabelError(
            "suffixes must be an iterable of values, not a single string"
        )
    if (
        not isinstance(separator, str)
        or "\n" in separator
        or "\r" in separator
        or "$" in separator
    ):
        raise TaxonLabelError(
            "separator must be a single-line plain-text string without '$'"
        )

    upright_suffixes = []
    if parental_lineage is not None:
        rendered = _upright_mathtext(
            f"parental lineage {parental_lineage}",
            field_name="biological_species parental lineage",
        )
        assert rendered is not None
        upright_suffixes.append(rendered)
    for index, suffix in enumerate(suffixes):
        if (
            isinstance(suffix, str)
            and " ".join(suffix.split()).casefold()
            in _OMITTED_PUBLICATION_SUFFIXES
        ):
            continue
        if (
            parental_lineage is not None
            and isinstance(suffix, str)
            and " ".join(suffix.split()) == parental_lineage
        ):
            continue
        rendered = _upright_mathtext(suffix, field_name=f"suffixes[{index}]")
        if rendered is not None:
            upright_suffixes.append(rendered)
    if not upright_suffixes:
        return species_label
    return f"{species_label} {separator.join(upright_suffixes)}"


def format_taxon_label_from_metadata(
    metadata: Mapping[str, Any],
    *,
    species_field: str = "biological_species",
    suffix_fields: Sequence[str] = DEFAULT_SUFFIX_FIELDS,
    abbreviate_genus: bool = False,
    separator: str = " | ",
) -> str:
    """Build a taxon label from one manifest or result-table row.

    Missing suffix columns are treated as empty so the same formatter can be
    used for species-level and assembly-unit-level tables.  The species column
    is required.  ``suffix_fields`` controls which metadata dimensions appear
    and in what order; no species or assembly identifier is hard-coded here.
    """

    if not isinstance(metadata, Mapping):
        raise TaxonLabelError("metadata must be a mapping such as a TSV row")
    if species_field not in metadata:
        raise TaxonLabelError(
            f"metadata is missing required species field {species_field!r}"
        )
    if not isinstance(suffix_fields, Sequence) or isinstance(
        suffix_fields, (str, bytes)
    ):
        raise TaxonLabelError(
            "suffix_fields must be a sequence of column names, not a string"
        )
    invalid_fields = [field for field in suffix_fields if not isinstance(field, str) or not field]
    if invalid_fields:
        raise TaxonLabelError("suffix_fields must contain only non-empty strings")

    return format_taxon_label(
        metadata[species_field],
        (metadata.get(field) for field in suffix_fields),
        abbreviate_genus=abbreviate_genus,
        separator=separator,
    )


def format_downstream_taxon_label(
    biological_species: str,
    suffixes: Iterable[str | None] = (),
    *,
    abbreviate_genus: bool = False,
    separator: str = " | ",
) -> str:
    """Return the concise, author-approved label used in downstream figures.

    Internal biological-lineage names and assembly identifiers remain unchanged.
    Only the publication display is simplified. Hybrid markers and the
    ``unphased`` assembly qualifier are omitted, while the two
    *A. zhejiangensis* parental lineages remain distinguishable by upright
    ``A``/``B`` suffixes. Other informative assembly-unit suffixes are retained.
    """

    genus, epithet, _hybrid, parental_lineage = _parse_biological_species(
        biological_species
    )
    normalized = " ".join(biological_species.split())
    if parental_lineage is not None:
        base_species = f"{genus} {epithet}"
        return format_taxon_label(
            base_species,
            (parental_lineage,),
            abbreviate_genus=abbreviate_genus,
            separator=separator,
        )
    return format_taxon_label(
        normalized,
        suffixes,
        abbreviate_genus=abbreviate_genus,
        separator=separator,
    )


def format_downstream_taxon_label_from_metadata(
    metadata: Mapping[str, Any],
    *,
    species_field: str = "biological_species",
    suffix_fields: Sequence[str] = DEFAULT_SUFFIX_FIELDS,
    abbreviate_genus: bool = False,
    separator: str = " | ",
) -> str:
    """Build a concise downstream label from a manifest or result-table row."""

    if not isinstance(metadata, Mapping):
        raise TaxonLabelError("metadata must be a mapping such as a TSV row")
    if species_field not in metadata:
        raise TaxonLabelError(
            f"metadata is missing required species field {species_field!r}"
        )
    if not isinstance(suffix_fields, Sequence) or isinstance(
        suffix_fields, (str, bytes)
    ):
        raise TaxonLabelError(
            "suffix_fields must be a sequence of column names, not a string"
        )
    invalid_fields = [
        field for field in suffix_fields if not isinstance(field, str) or not field
    ]
    if invalid_fields:
        raise TaxonLabelError("suffix_fields must contain only non-empty strings")
    return format_downstream_taxon_label(
        metadata[species_field],
        (metadata.get(field) for field in suffix_fields),
        abbreviate_genus=abbreviate_genus,
        separator=separator,
    )
