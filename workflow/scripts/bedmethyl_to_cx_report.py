#!/usr/bin/env python3
"""Convert modkit bedMethyl format to Bismark CX_report format.

Determines methylation context (CG/CHG/CHH) from the reference genome.

bedMethyl columns (modkit output):
    0: chrom
    1: start (0-based)
    2: end
    3: base (C)
    4: score
    5: strand (. = combined, + or -)
    6-8: thickStart, thickEnd, itemRgb
    9: coverage
    10: percent methylated
    11: N_mod (modified count)
    12: N_canonical (unmodified count)

CX_report format (Bismark):
    0: chrom
    1: position (1-based)
    2: strand (+/-)
    3: count methylated
    4: count unmethylated
    5: context (CG/CHG/CHH)
    6: trinucleotide sequence

Usage:
    python bedmethyl_to_cx_report.py input.bed.gz reference.fa output.CX_report.txt.gz
"""

import sys
import gzip
import pysam
from contextlib import ExitStack


def get_context_and_trinuc_from_seq(chrom_seq, chrom_len, pos, strand):
    """Determine methylation context and trinucleotide from cached chromosome sequence.

    Args:
        chrom_seq: uppercase chromosome sequence string
        chrom_len: length of chromosome
        pos: 0-based position of the C
        strand: '+' or '-'

    Returns:
        (context, trinucleotide) tuple
    """
    if strand == '+':
        # Forward strand: look at C and next 2 bases
        end = min(pos + 3, chrom_len)
        seq = chrom_seq[pos:end]
        if len(seq) < 2:
            return ("CHH", seq + "N" * (3 - len(seq)))

        trinuc = seq[:3] if len(seq) >= 3 else seq + "N" * (3 - len(seq))

        # Determine context
        if len(seq) >= 2 and seq[1] == 'G':
            context = "CG"
        elif len(seq) >= 3 and seq[2] == 'G':
            context = "CHG"
        else:
            context = "CHH"
    else:
        # Reverse strand: C is actually a G on forward strand
        # Look at 2 bases before and the G position
        start = max(0, pos - 2)
        seq = chrom_seq[start:pos + 1]

        # Reverse complement
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        seq_rc = ''.join(comp.get(b, 'N') for b in reversed(seq))

        if len(seq_rc) < 3:
            seq_rc = "N" * (3 - len(seq_rc)) + seq_rc
        trinuc = seq_rc[:3]

        # Determine context (same logic as forward, but on rev comp)
        if len(trinuc) >= 2 and trinuc[1] == 'G':
            context = "CG"
        elif len(trinuc) >= 3 and trinuc[2] == 'G':
            context = "CHG"
        else:
            context = "CHH"

    return (context, trinuc)


def convert_bedmethyl_to_cx_report(input_file, reference_file, output_file):
    """Convert bedMethyl to CX_report format with context from reference.

    Optimized to load each chromosome sequence into memory once, avoiding
    per-position file I/O overhead from repeated fasta.fetch() calls.
    """

    # Open reference genome
    fasta = pysam.FastaFile(reference_file)

    # Cache for current chromosome sequence (avoids millions of fetch calls)
    current_chrom = None
    chrom_seq = None
    chrom_len = 0

    with ExitStack() as stack:
        # Handle gzipped or plain input
        if input_file.endswith('.gz'):
            infile = stack.enter_context(gzip.open(input_file, 'rt'))
        else:
            infile = stack.enter_context(open(input_file, 'r'))

        # Handle stdout or file output
        if output_file == '/dev/stdout':
            outfile = sys.stdout
        elif output_file.endswith('.gz'):
            outfile = stack.enter_context(gzip.open(output_file, 'wt'))
        else:
            outfile = stack.enter_context(open(output_file, 'w'))

        for line in infile:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 13:
                continue

            chrom = fields[0]
            start = int(fields[1])
            strand = fields[5]
            n_mod = int(fields[11])
            n_canonical = int(fields[12])

            # Load chromosome sequence if we've moved to a new chromosome
            if chrom != current_chrom:
                try:
                    chrom_seq = fasta.fetch(chrom).upper()
                    chrom_len = len(chrom_seq)
                    current_chrom = chrom
                    print(f"Loaded chromosome {chrom} ({chrom_len:,} bp)", file=sys.stderr)
                except KeyError:
                    print(f"Warning: chromosome {chrom} not found in reference", file=sys.stderr)
                    chrom_seq = ""
                    chrom_len = 0
                    current_chrom = chrom

            # Handle combined strand data (strand = '.')
            # For combined CpG, report as + strand
            if strand == '.':
                strand = '+'

            # Get context and trinucleotide from cached sequence
            if chrom_len == 0:
                context, trinuc = ("CNN", "CNN")
            else:
                context, trinuc = get_context_and_trinuc_from_seq(chrom_seq, chrom_len, start, strand)

            # Convert to 1-based position for CX_report
            position = start + 1

            # Write CX_report format
            outfile.write(f"{chrom}\t{position}\t{strand}\t{n_mod}\t{n_canonical}\t{context}\t{trinuc}\n")

    fasta.close()


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} input.bed.gz reference.fa output.CX_report.txt.gz")
        sys.exit(1)

    input_file = sys.argv[1]
    reference_file = sys.argv[2]
    output_file = sys.argv[3]

    convert_bedmethyl_to_cx_report(input_file, reference_file, output_file)


if __name__ == "__main__":
    main()
