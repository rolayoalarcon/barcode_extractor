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
    matched allowing up to `max_errors` substitutions/indels.
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
 
    return _re.compile(up + barcode + down, _re.IGNORECASE)

def find_barcode(seq, pattern):
    """Search forward strand then reverse complement.
 
    Returns (barcode, strand) or (None, None) if no match on either strand.
    """
    candidates = (("+", seq), ("-", seq.reverse_complement()))
    for strand, oriented in candidates:
        s = str(oriented)
        # BESTMATCH (regex module) returns the match with the fewest errors,
        # which matters when fuzzy matching could find several alignments.
        if HAVE_FUZZY:
            #m = pattern.search(s, _re.BESTMATCH)
            m = pattern.search(s)
        else:
            m = pattern.search(s)
        if m:
            return m.group("bc").upper(), strand
    return None, None


def extract_barcodes(input_fasta: Path, 
                     output_prefix: str, 
                     upstream_seq: str,
                     downstream_seq: str,
                     min_len: int = 10,
                     max_len: int = 30,
                     max_errors: int = 2) -> Path:
    """
    Extracts the first `barcode_length` bases from each sequence in the input FASTA file
    and writes them to the output FASTA file.

    Args:
        input_fasta (Path): Path to the input FASTA file.
        output_prefix (str): prefix for the output FASTA and tsv file where extracted barcodes will be saved.
        upstream_seq (str): The upstream sequence to match.
        downstream_seq (str): The downstream sequence to match.
        min_len (int): The minimum length of the barcode.
        max_len (int): The maximum length of the barcode.
        max_errors (int): The maximum number of errors allowed in the flanks.

    Returns:
        Path: Path to the output FASTA file.
    """

    # Build the regex pattern for barcode extraction
    pattern = build_pattern(upstream_seq, downstream_seq, min_len, max_len, max_errors)

    # Determine the output directory and create it if it doesn't exist
    OUTPUT_DIRECTORY = input_fasta.parent

    out = open(input_fasta.parent / f"{output_prefix}_barcodes.tsv", "w")
    fa_out = open(input_fasta.parent / f"{output_prefix}_barcodes.fasta", "w")

    n_total = 0
    n_found = 0
    try:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            n_total += 1
            bc, strand = find_barcode(rec.seq, pattern)
            if bc is None:
                # No anchor found on either strand; record as NA so you can
                # audit your recovery rate rather than silently dropping reads.
                out.write(f"{rec.id}\tNA\tNA\n")
                continue
            n_found += 1
            out.write(f"{rec.id}\t{strand}\t{bc}\n")
            if fa_out:
                fa_out.write(f">{rec.id} strand={strand}\n{bc}\n")
    finally:
        if out is not sys.stdout:
            out.close()
        if fa_out:
            fa_out.close()

    