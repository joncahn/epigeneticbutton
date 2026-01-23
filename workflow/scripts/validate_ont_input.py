#!/usr/bin/env python3
"""
Validation script for ONT methylation inputs.

Validates:
- modBAM files: checks for MM/ML tags and read count > 0
- bedMethyl files: validates column count and format

Usage:
    python validate_ont_input.py <input_type> <input_file> [<chrom_sizes>]

Arguments:
    input_type: Either "modBAM" or "bedMethyl"
    input_file: Path to the input file
    chrom_sizes: Path to chrom.sizes file (optional, for reference validation)

Exit codes:
    0: Validation passed
    1: Validation failed
"""

import sys
import subprocess
import gzip
from pathlib import Path


def validate_modbam(bam_path: str, chrom_sizes: str = None) -> tuple[bool, str]:
    """
    Validate a modBAM file for MM/ML tags.

    Args:
        bam_path: Path to the BAM file
        chrom_sizes: Optional path to chrom.sizes for reference validation

    Returns:
        Tuple of (is_valid, message)
    """
    bam_path = Path(bam_path)

    if not bam_path.exists():
        return False, f"File not found: {bam_path}"

    # Check for MM/ML tags in first 1000 reads
    try:
        result = subprocess.run(
            ["samtools", "view", str(bam_path)],
            capture_output=True,
            text=True,
            timeout=120
        )
        if result.returncode != 0:
            return False, f"Failed to read BAM file: {result.stderr}"

        lines = result.stdout.strip().split('\n')
        if not lines or lines == ['']:
            return False, "BAM file contains no reads"

        # Check for MM tag in reads (sample first 1000)
        mm_found = False
        ml_found = False
        read_count = 0

        for line in lines[:1000]:
            if not line or line.startswith('@'):
                continue
            read_count += 1
            if '\tMM:' in line or '\tMM:Z:' in line:
                mm_found = True
            if '\tML:' in line or '\tML:B:' in line:
                ml_found = True
            if mm_found and ml_found:
                break

        if read_count == 0:
            return False, "BAM file contains no aligned reads"

        if not mm_found:
            return False, "No MM (methylation) tags found in BAM file. Ensure this is a modBAM with base modifications."

        if not ml_found:
            return False, "No ML (methylation likelihood) tags found in BAM file. Ensure this is a modBAM with base modifications."

    except subprocess.TimeoutExpired:
        return False, "Timeout while reading BAM file"
    except Exception as e:
        return False, f"Error validating BAM file: {str(e)}"

    # Check reference alignment if chrom.sizes provided
    if chrom_sizes:
        try:
            # Get BAM header chromosomes
            result = subprocess.run(
                ["samtools", "view", "-H", str(bam_path)],
                capture_output=True,
                text=True
            )
            if result.returncode != 0:
                return False, f"Failed to read BAM header: {result.stderr}"

            bam_chroms = set()
            for line in result.stdout.split('\n'):
                if line.startswith('@SQ'):
                    parts = line.split('\t')
                    for part in parts:
                        if part.startswith('SN:'):
                            bam_chroms.add(part[3:])

            # Get reference chromosomes
            ref_chroms = set()
            with open(chrom_sizes) as f:
                for line in f:
                    if line.strip():
                        ref_chroms.add(line.split('\t')[0])

            # Check overlap
            overlap = bam_chroms & ref_chroms
            if not overlap:
                return False, f"No matching chromosomes between BAM and reference. BAM has: {list(bam_chroms)[:5]}..., Reference has: {list(ref_chroms)[:5]}..."

            coverage = len(overlap) / len(ref_chroms) * 100
            if coverage < 50:
                return False, f"Only {coverage:.1f}% of reference chromosomes found in BAM. Realignment may be needed."

        except Exception as e:
            return False, f"Error checking reference alignment: {str(e)}"

    return True, f"Valid modBAM with MM/ML tags ({read_count}+ reads checked)"


def validate_bedmethyl(bedmethyl_path: str, chrom_sizes: str = None) -> tuple[bool, str]:
    """
    Validate a bedMethyl file format.

    bedMethyl format (11 columns):
    1. chrom
    2. start (0-based)
    3. end
    4. name (modified base code)
    5. score
    6. strand
    7. start (thick)
    8. end (thick)
    9. color
    10. coverage (Nvalid)
    11. percent modified

    Args:
        bedmethyl_path: Path to the bedMethyl file
        chrom_sizes: Optional path to chrom.sizes for coordinate validation

    Returns:
        Tuple of (is_valid, message)
    """
    bedmethyl_path = Path(bedmethyl_path)

    if not bedmethyl_path.exists():
        return False, f"File not found: {bedmethyl_path}"

    # Load reference chromosomes if provided
    ref_chroms = {}
    if chrom_sizes:
        try:
            with open(chrom_sizes) as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        ref_chroms[parts[0]] = int(parts[1])
        except Exception as e:
            return False, f"Error reading chrom.sizes: {str(e)}"

    # Open file (handle gzip)
    try:
        if str(bedmethyl_path).endswith('.gz'):
            fh = gzip.open(bedmethyl_path, 'rt')
        else:
            fh = open(bedmethyl_path, 'r')
    except Exception as e:
        return False, f"Error opening file: {str(e)}"

    try:
        line_count = 0
        valid_lines = 0
        chroms_found = set()

        for line in fh:
            line_count += 1
            if line_count > 10000:  # Sample first 10000 lines
                break

            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track'):
                continue

            fields = line.split('\t')

            # Check column count (bedMethyl should have at least 10-11 columns)
            if len(fields) < 10:
                return False, f"Line {line_count}: Expected at least 10 columns, got {len(fields)}"

            chrom = fields[0]
            chroms_found.add(chrom)

            # Validate coordinates
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                return False, f"Line {line_count}: Invalid coordinates (start={fields[1]}, end={fields[2]})"

            if start < 0:
                return False, f"Line {line_count}: Start position cannot be negative ({start})"

            if end <= start:
                return False, f"Line {line_count}: End ({end}) must be greater than start ({start})"

            # Validate strand
            if fields[5] not in ['+', '-', '.']:
                return False, f"Line {line_count}: Invalid strand '{fields[5]}' (expected +, -, or .)"

            # Validate coverage (column 10) and percent (column 11 if present)
            try:
                coverage = int(fields[9])
                if coverage < 0:
                    return False, f"Line {line_count}: Coverage cannot be negative ({coverage})"
            except ValueError:
                return False, f"Line {line_count}: Invalid coverage value '{fields[9]}'"

            if len(fields) >= 11:
                try:
                    percent = float(fields[10])
                    if percent < 0 or percent > 100:
                        return False, f"Line {line_count}: Percent modified should be 0-100, got {percent}"
                except ValueError:
                    return False, f"Line {line_count}: Invalid percent modified value '{fields[10]}'"

            # Check chromosome against reference
            if ref_chroms and chrom not in ref_chroms:
                # Only warn, don't fail - there might be valid reasons
                pass
            elif ref_chroms and end > ref_chroms.get(chrom, float('inf')):
                return False, f"Line {line_count}: Position {end} exceeds chromosome {chrom} length {ref_chroms[chrom]}"

            valid_lines += 1

        if valid_lines == 0:
            return False, "No valid data lines found in bedMethyl file"

        # Check chromosome overlap with reference
        if ref_chroms:
            overlap = chroms_found & set(ref_chroms.keys())
            if not overlap:
                return False, f"No matching chromosomes between bedMethyl and reference"

    finally:
        fh.close()

    return True, f"Valid bedMethyl format ({valid_lines} lines validated, {len(chroms_found)} chromosomes)"


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    input_type = sys.argv[1]
    input_file = sys.argv[2]
    chrom_sizes = sys.argv[3] if len(sys.argv) > 3 else None

    if input_type == "modBAM":
        is_valid, message = validate_modbam(input_file, chrom_sizes)
    elif input_type == "bedMethyl":
        is_valid, message = validate_bedmethyl(input_file, chrom_sizes)
    else:
        print(f"Unknown input type: {input_type}. Expected 'modBAM' or 'bedMethyl'")
        sys.exit(1)

    if is_valid:
        print(f"PASS: {message}")
        sys.exit(0)
    else:
        print(f"FAIL: {message}")
        sys.exit(1)


if __name__ == "__main__":
    main()
