"""Primary-isoform sequence extraction for standardized loss-analysis inputs."""

from __future__ import annotations

from pathlib import Path

from .gff import Transcript, collect_transcripts, load_fasta, select_longest_cds_isoform
from .io_utils import SchemaError, write_tsv


_COMPLEMENT = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")

# Plant nuclear genes use the standard genetic code.  Unknown/ambiguous codons
# become X; this is recorded in QC instead of silently discarded.
_STANDARD_CODE = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L", "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*", "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M", "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V", "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _cds_phase(feature_phase: str, transcript_id: str) -> int:
    """Return a validated GFF3 CDS phase.

    GFF3 requires a phase of 0, 1, or 2 for CDS rows.  A small number of
    publisher files use ``.``; treating that value as zero preserves the old
    compatibility behavior, while the extraction QC still records malformed
    models through the independent publisher-protein comparison.
    """

    phase = 0 if feature_phase in {".", ""} else int(feature_phase)
    if phase not in {0, 1, 2}:
        raise SchemaError(f"{transcript_id}: invalid CDS phase {feature_phase!r}")
    return phase


def build_spliced_cds(transcript: Transcript, genome: dict[str, str]) -> str:
    """Extract a CDS in transcript orientation and validate its GFF3 phases.

    Phase is *not* a per-exon trimming instruction.  For an internal CDS row,
    its leading ``phase`` bases complete the codon started in the preceding
    row and therefore remain in the spliced sequence.  Only a non-zero phase
    on the first 5'-most CDS row is removed, because those bases precede the
    first complete codon in a partial model.  Every later phase is checked
    against the preceding row using the Sequence Ontology GFF3 recurrence.
    """
    if transcript.chrom not in genome:
        raise SchemaError(f"{transcript.transcript_id}: chromosome {transcript.chrom!r} not in genome FASTA")
    ordered = sorted(transcript.cds_features, key=lambda item: item.start, reverse=transcript.strand == "-")
    pieces: list[str] = []
    initial_phase: int | None = None
    previous_phase: int | None = None
    previous_length: int | None = None
    for feature_index, feature in enumerate(ordered):
        if feature.sequence_id not in genome:
            raise SchemaError(
                f"{transcript.transcript_id}: CDS sequence {feature.sequence_id!r} not in genome FASTA"
            )
        sequence = genome[feature.sequence_id][feature.start - 1:feature.end]
        if len(sequence) != feature.end - feature.start + 1:
            raise SchemaError(
                f"{transcript.transcript_id}: GFF CDS interval {feature.sequence_id}:{feature.start}-{feature.end} "
                "exceeds FASTA sequence length"
            )
        if transcript.strand == "-":
            sequence = reverse_complement(sequence)
        phase = _cds_phase(feature.phase, transcript.transcript_id)
        if feature_index == 0:
            initial_phase = phase
            if phase >= len(sequence):
                raise SchemaError(f"{transcript.transcript_id}: initial CDS phase removes entire segment")
        else:
            assert previous_phase is not None and previous_length is not None
            expected_phase = (3 - ((previous_length - previous_phase) % 3)) % 3
            if phase != expected_phase:
                raise SchemaError(
                    f"{transcript.transcript_id}: inconsistent CDS phase chain; "
                    f"expected {expected_phase}, found {phase} at "
                    f"{feature.sequence_id}:{feature.start}-{feature.end}"
                )
        pieces.append(sequence)
        previous_phase = phase
        previous_length = len(sequence)
    assert initial_phase is not None
    return "".join(pieces)[initial_phase:].upper()


def translate_standard(cds: str) -> str:
    usable = len(cds) - (len(cds) % 3)
    return "".join(_STANDARD_CODE.get(cds[index:index + 3], "X") for index in range(0, usable, 3))


def write_fasta(path: str | Path, records: list[tuple[str, str]], line_width: int = 60) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), line_width):
                handle.write(f"{sequence[start:start + line_width]}\n")


def extract_annotation(
    genome_fasta: str | Path,
    gff3: str | Path,
    output_dir: str | Path,
    sample_id: str,
) -> dict[str, Path]:
    """Create protein, CDS, coordinates, primary-isoform map and QC tables.

    Files retain transcript IDs in the FASTA/legacy coord file because SynOrths
    and the archived result files use transcript IDs.  ``isoform_map.tsv`` is
    therefore mandatory for translating back to gene-level identifiers.
    """
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    genome = load_fasta(genome_fasta)
    transcripts, audit = collect_transcripts(gff3)
    selected = select_longest_cds_isoform(transcripts)
    if not selected:
        raise SchemaError(f"{gff3}: no transcript with usable CDS was found")

    proteins: list[tuple[str, str]] = []
    cdss: list[tuple[str, str]] = []
    isoform_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = list(audit)
    coordinates: list[dict[str, object]] = []
    for transcript in selected:
        try:
            cds = build_spliced_cds(transcript, genome)
        except SchemaError as exc:
            qc_rows.append({
                "record_type": "transcript", "record_id": transcript.transcript_id,
                "reason": str(exc), "line_number": "",
            })
            continue
        protein = translate_standard(cds)
        status = "ok"
        if len(cds) % 3:
            status = "cds_length_not_divisible_by_3_trimmed_for_translation"
        if "*" in protein[:-1]:
            status = f"{status};internal_stop" if status != "ok" else "internal_stop"
        if "X" in protein:
            status = f"{status};ambiguous_codon" if status != "ok" else "ambiguous_codon"
        proteins.append((transcript.transcript_id, protein.rstrip("*")))
        cdss.append((transcript.transcript_id, cds))
        coordinates.append({
            "transcript_id": transcript.transcript_id,
            "gene_id": transcript.gene_id,
            "chromosome": transcript.chrom,
            "start": transcript.start,
            "end": transcript.end,
            "strand": transcript.strand,
        })
        isoform_rows.append({
            "sample_id": sample_id,
            "gene_id": transcript.gene_id,
            "transcript_id": transcript.transcript_id,
            "chromosome": transcript.chrom,
            "start": transcript.start,
            "end": transcript.end,
            "strand": transcript.strand,
            "spliced_cds_length": len(cds),
            "genomic_span": transcript.genomic_span,
            "protein_length": len(protein.rstrip("*")),
            "selection_rule": "longest_spliced_CDS_then_genomic_span_then_transcript_id",
            "qc_status": status,
        })

    if not proteins:
        raise SchemaError("all selected transcripts failed sequence extraction; inspect extraction_qc.tsv")
    protein_path = output / f"{sample_id}.protein.faa"
    cds_path = output / f"{sample_id}.cds.fa"
    coords_path = output / f"{sample_id}.coords.tsv"
    legacy_coords_path = output / f"{sample_id}.coords"
    isoform_path = output / f"{sample_id}.isoform_map.tsv"
    qc_path = output / f"{sample_id}.extraction_qc.tsv"
    write_fasta(protein_path, proteins)
    write_fasta(cds_path, cdss)
    write_tsv(coords_path, coordinates, ["transcript_id", "gene_id", "chromosome", "start", "end", "strand"])
    # SynOrths V1.5 expects five headerless columns.  Do not use this as the
    # canonical metadata; it is a compatibility artifact only.
    with open(legacy_coords_path, "w", encoding="utf-8") as handle:
        for row in coordinates:
            handle.write(
                f"{row['transcript_id']}\t{row['chromosome']}\t{row['start']}\t{row['end']}\t{row['strand']}\n"
            )
    write_tsv(
        isoform_path, isoform_rows,
        ["sample_id", "gene_id", "transcript_id", "chromosome", "start", "end", "strand", "spliced_cds_length", "genomic_span", "protein_length", "selection_rule", "qc_status"],
    )
    write_tsv(qc_path, qc_rows, ["record_type", "record_id", "reason", "line_number"])
    return {
        "protein": protein_path, "cds": cds_path, "coords": coords_path,
        "legacy_coords": legacy_coords_path, "isoform_map": isoform_path, "qc": qc_path,
    }
