#!/usr/bin/env python3
"""Convert modkit bedMethyl format to Bismark CX_report format for DMRcaller.

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
    python bedmethyl_to_cx_report.py input.bed.gz output.CX_report.txt.gz context
"""

import sys
import gzip
from contextlib import ExitStack


def convert_bedmethyl_to_cx_report(input_file, output_file, context):
    """Convert bedMethyl to CX_report format."""

    # Generate placeholder trinucleotide based on context
    if context == "CG":
        trinuc = "CGA"
    elif context == "CHG":
        trinuc = "CAG"
    elif context == "CHH":
        trinuc = "CAA"
    else:
        trinuc = "CNN"

    with ExitStack() as stack:
        # Handle gzipped or plain input
        if input_file.endswith('.gz'):
            infile = stack.enter_context(gzip.open(input_file, 'rt'))
        else:
            infile = stack.enter_context(open(input_file, 'r'))

        # Handle gzipped or plain output
        if output_file.endswith('.gz'):
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

            # Convert to 1-based position
            position = start + 1

            # Handle combined strand data (strand = '.')
            # For combined CpG, report as + strand
            if strand == '.':
                strand = '+'

            # Write CX_report format
            outfile.write(f"{chrom}\t{position}\t{strand}\t{n_mod}\t{n_canonical}\t{context}\t{trinuc}\n")


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} input.bed.gz output.CX_report.txt.gz context")
        print("  context: CG, CHG, or CHH")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    context = sys.argv[3].upper()

    if context not in ["CG", "CHG", "CHH"]:
        print(f"Error: context must be CG, CHG, or CHH (got '{context}')")
        sys.exit(1)

    convert_bedmethyl_to_cx_report(input_file, output_file, context)


if __name__ == "__main__":
    main()
