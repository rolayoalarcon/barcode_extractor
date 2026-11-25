import sys
from Bio import SeqIO
from pathlib import Path
# Prefer the third-party `regex` module for fuzzy matching; fall back to
# stdlib `re` (exact only) so the script still runs without it.
try:
    import regex as _re
    HAVE_FUZZY = True
except ImportError:  # pragma: no cover
    import re as _re
    HAVE_FUZZY = False


def build_pattern(upstream, downstream, min_len, max_len, max_errors):
    """Compile a pattern that captures the barcode between the two flanks.

    The barcode group is bounded to [min_len, max_len] of nucleotide
    characters. With the `regex` module and max_errors > 0, each flank is
    matched allowing up to `max_errors` substitutions.
    """
    barcode = f"(?P<bc>[ACGTacgt]{{{min_len},{max_len}}})"

    if HAVE_FUZZY and max_errors > 0:
        up = f"(?:{upstream}){{s<={max_errors}}}"
        down = f"(?:{downstream}){{s<={max_errors}}}"
    else:
        if max_errors > 0 and not HAVE_FUZZY:
            sys.stderr.write(
                "WARNING: `regex` module not installed; falling back to "
                "EXACT flank matching (--max-errors ignored).\n"
            )
        up, down = upstream, downstream

    flags = _re.IGNORECASE
    if HAVE_FUZZY:
        # Without BESTMATCH, `regex` returns the LEFTMOST match within the error
        # budget, which can be a worse alignment that shifts the flank and pulls
        # flank bases into the barcode. BESTMATCH returns the fewest-error match.
        # NB: it must be a compile flag; `pattern.search(s, BESTMATCH)` would
        # pass it as the `pos` argument (= 4096) and search from that offset.
        flags |= _re.BESTMATCH
    return _re.compile(up + barcode + down, flags)


def find_barcode(seq, pattern):
    """Search forward strand then reverse complement.

    Returns (barcode, strand, start, end) or (None, None, None, None).

    `start`/`end` are 0-based, end-exclusive coordinates on the ORIGINAL read
    (as written in the FASTA), so `str(seq[start:end])` is the barcode region
    on the forward strand. For '-' hits the reported barcode is the reverse
    complement of that region (i.e. read in the orientation it was matched).
    """
    read_len = len(seq)
    candidates = (("+", seq), ("-", seq.reverse_complement()))
    for strand, oriented in candidates:
        s = str(oriented)
        m = pattern.search(s)  # BESTMATCH is set at compile time (see build_pattern)
        if m:
            # Span of the barcode group itself, in the searched orientation
            s_start, s_end = m.span("bc")
            if strand == "+":
                start, end = s_start, s_end
            else:
                # Map reverse-complement coordinates back to the original read:
                # position i in revcomp <-> position (L - 1 - i) in the read,
                # so the half-open interval [a, b) becomes [L - b, L - a).
                start, end = read_len - s_end, read_len - s_start
            return m.group("bc").upper(), strand, start, end
    return None, None, None, None


def extract_barcodes(input_fasta: Path,
                     output_prefix: str,
                     upstream_seq: str,
                     downstream_seq: str,
                     min_len: int = 10,
                     max_len: int = 30,
                     max_errors: int = 2) -> Path:
    """
    Extract the barcode between `upstream_seq` and `downstream_seq` from each
    sequence in the input FASTA, writing a TSV report and a FASTA of barcodes.

    TSV columns:
        read_id, strand, barcode, bc_start, bc_end
    where bc_start/bc_end are 0-based, end-exclusive positions on the original
    read (Python slice / BED convention; barcode length = bc_end - bc_start).

    Args:
        input_fasta (Path): Path to the input FASTA file.
        output_prefix (str): prefix for the output FASTA and tsv file where extracted barcodes will be saved.
        upstream_seq (str): The upstream sequence to match.
        downstream_seq (str): The downstream sequence to match.
        min_len (int): The minimum length of the barcode.
        max_len (int): The maximum length of the barcode.
        max_errors (int): The maximum number of substitutions allowed in each flank.

    Returns:
        Path: Path to the output TSV file.
    """

    # Build the regex pattern for barcode extraction
    pattern = build_pattern(upstream_seq, downstream_seq, min_len, max_len, max_errors)

    tsv_path = input_fasta.parent / f"{output_prefix}_barcodes.tsv"
    fasta_path = input_fasta.parent / f"{output_prefix}_barcodes.fasta"

    n_total = 0
    n_found = 0
    with open(tsv_path, "w") as out, open(fasta_path, "w") as fa_out:
        out.write("read_id\tstrand\tbarcode\tbc_start\tbc_end\n")
        for rec in SeqIO.parse(input_fasta, "fasta"):
            n_total += 1
            bc, strand, start, end = find_barcode(rec.seq, pattern)
            if bc is None:
                # No anchor found on either strand; record as NA so you can
                # audit your recovery rate rather than silently dropping reads.
                out.write(f"{rec.id}\tNA\tNA\tNA\tNA\n")
                continue
            n_found += 1
            out.write(f"{rec.id}\t{strand}\t{bc}\t{start}\t{end}\n")
            fa_out.write(f">{rec.id} strand={strand} start={start} end={end}\n{bc}\n")

    print(f"Barcodes found in {n_found}/{n_total} reads -> {tsv_path}")
    return tsv_path