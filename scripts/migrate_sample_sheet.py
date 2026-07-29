#!/usr/bin/env python3
"""Convert an old-format sample sheet to the new format.

Old columns: data_type, line, tissue, sample_type, replicate, seq_id, fastq_path, paired, ref_genome
New columns: Sample_ID, Assay, Genome, Levels, Replicate_ID, Read_files, Read_layout, IP_target, Control

Usage:
    python scripts/migrate_sample_sheet.py old_samples.tsv > new_samples.tsv
    python scripts/migrate_sample_sheet.py old_samples.tsv -o new_samples.tsv
"""

import argparse
import sys
import pandas as pd


OLD_COLNAMES = [
    "data_type", "line", "tissue", "sample_type", "replicate",
    "seq_id", "fastq_path", "paired", "ref_genome",
]

# Map old data_type + sample_type → new Assay
# ChIP samples: sample_type determines broad/narrow based on legacy peaktype patterns
# TF_ samples: historically always narrow
_BROAD_MARKS = {
    "H3K9me2", "H3K9me3", "H3K27me1", "H3K27me3",
    "H3K4me1", "H3K36me3", "H3K79me2",
    "CenH3",
}

_NARROW_MARKS = {
    "H3K4me3", "H3K27ac", "H3K9ac", "H3K4me2",
}


def infer_assay(data_type, sample_type):
    """Infer the new Assay value from old data_type and sample_type."""
    if data_type == "RNAseq":
        return "RNAseq"
    if data_type == "RAMPAGE":
        return "RAMPAGE"
    if data_type in ("sRNA", "smallRNA", "shRNA"):
        return "sRNA"
    if data_type == "ATAC":
        return "ATAC"
    if data_type == "mC" or data_type.startswith("mC"):
        if sample_type in ("dmC", "bedMethyl"):
            return "dmC"
        elif sample_type == "WGBS":
            return "WGBS"
        elif sample_type == "EMseq":
            return "EMseq"
        elif sample_type == "Pico":
            return "WGBS"  # Pico → WGBS
        else:
            return "WGBS"  # default for mC

    # ChIP / TF
    if data_type.startswith("TF_"):
        # TF samples are always narrow
        return "ChIP_narrow"
    if data_type.startswith("ChIP"):
        # Determine broad/narrow from sample_type
        if sample_type == "Input":
            return "ChIP_broad"  # placeholder; will be updated based on paired IP
        if sample_type in ("IP", "IPb"):
            return "ChIP_narrow" if sample_type == "IP" else "ChIP_broad"
        if sample_type in _BROAD_MARKS:
            return "ChIP_broad"
        if sample_type in _NARROW_MARKS:
            return "ChIP_narrow"
        # Default: broad (most histone marks are broad)
        return "ChIP_broad"

    raise ValueError(f"Cannot infer Assay from data_type='{data_type}', sample_type='{sample_type}'")


def infer_ip_target(data_type, sample_type):
    """Infer IP_target from old fields."""
    if data_type.startswith("TF_"):
        # TF_<name>: IP_target = sample_type (IP, IPb, Input)
        if sample_type == "Input":
            return "Input"
        return data_type.split("_", 1)[1]  # e.g. "TB1"
    if data_type.startswith("ChIP"):
        return sample_type  # H3K9me2, Input, etc.
    return ""  # non-ChIP assays


def build_read_files(seq_id, fastq_path):
    """Convert old seq_id + fastq_path into new Read_files format."""
    if str(fastq_path) == "SRA":
        # SRA accession(s) — comma → plus for merging
        return str(seq_id).replace(",", "+")
    else:
        # Local path: encode as "directory/seq_id" so the bridge function
        # can reconstruct seq_id and fastq_path for the download rules
        import os
        fpath = str(fastq_path)
        sid = str(seq_id)
        if fpath.endswith(".bam") or fpath.endswith(".bed.gz"):
            # BAM or bedMethyl: fastq_path IS the file path
            return fpath
        else:
            # FASTQ directory: combine path and seq_id
            return os.path.join(fpath, sid)


def generate_sample_id(row):
    """Auto-generate a concise Sample_ID from old fields."""
    parts = [row["line"], row["tissue"]]

    if row["data_type"].startswith("TF_"):
        tf_name = row["data_type"].split("_", 1)[1]
        if row["sample_type"] == "Input":
            parts.append("Input")
        else:
            parts.append(tf_name)
    elif row["data_type"].startswith("ChIP"):
        parts.append(row["sample_type"])
    elif row["data_type"] == "ATAC":
        parts.append("ATAC")
    elif row["data_type"] == "RNAseq":
        parts.append("RNA")
    elif row["data_type"] == "RAMPAGE":
        parts.append("RAMPAGE")
    elif row["data_type"] in ("sRNA", "smallRNA", "shRNA"):
        parts.append("sRNA")
    elif row["data_type"] == "mC" or row["data_type"].startswith("mC"):
        parts.append(row["sample_type"])
    else:
        parts.append(row["data_type"])

    parts.append(str(row["replicate"]))

    return "_".join(parts)


def migrate(old_df):
    """Convert an old-format DataFrame to new-format."""
    rows = []
    # First pass: build rows without Control
    sample_id_map = {}  # (data_type, line, tissue, sample_type, replicate, ref_genome) → Sample_ID

    for _, row in old_df.iterrows():
        sid = generate_sample_id(row)
        key = (row["data_type"], row["line"], row["tissue"],
               row["sample_type"], row["replicate"], row["ref_genome"])
        sample_id_map[key] = sid

    for _, row in old_df.iterrows():
        sid = generate_sample_id(row)
        assay = infer_assay(row["data_type"], row["sample_type"])
        ip_target = infer_ip_target(row["data_type"], row["sample_type"])
        read_files = build_read_files(row["seq_id"], row["fastq_path"])

        new_row = {
            "Sample_ID": sid,
            "Assay": assay,
            "Genome": row["ref_genome"],
            "Levels": f"genotype:{row['line']},tissue:{row['tissue']}",
            "Replicate_ID": str(row["replicate"]),
            "Read_files": read_files,
            "Read_layout": row["paired"],
            "IP_target": ip_target,
            "Control": "",
        }
        rows.append(new_row)

    new_df = pd.DataFrame(rows)

    # Second pass: populate Control column
    # For ChIP/TF IP samples: find matching Input
    for idx, row in old_df.iterrows():
        if row["data_type"].startswith("ChIP") or row["data_type"].startswith("TF_"):
            if row["sample_type"] != "Input":
                # Find the matching Input
                input_key = (row["data_type"], row["line"], row["tissue"],
                             "Input", row["replicate"], row["ref_genome"])
                if input_key in sample_id_map:
                    new_df.at[idx, "Control"] = sample_id_map[input_key]
                else:
                    # Try any replicate of Input
                    for rep_key, rep_sid in sample_id_map.items():
                        if (rep_key[0] == row["data_type"] and
                            rep_key[1] == row["line"] and
                            rep_key[2] == row["tissue"] and
                            rep_key[3] == "Input" and
                            rep_key[5] == row["ref_genome"]):
                            new_df.at[idx, "Control"] = rep_sid
                            break

    # Fix Input Assay to match paired IP's Assay
    for idx, row in new_df.iterrows():
        if row["IP_target"] == "Input":
            # Find any IP that references this as control
            referencing = new_df[new_df["Control"] == row["Sample_ID"]]
            if not referencing.empty:
                new_df.at[idx, "Assay"] = referencing.iloc[0]["Assay"]

    return new_df


def main():
    parser = argparse.ArgumentParser(description="Migrate old-format sample sheet to new format")
    parser.add_argument("input", help="Old-format sample sheet TSV")
    parser.add_argument("-o", "--output", help="Output file (default: stdout)")
    args = parser.parse_args()

    # Read old format
    with open(args.input) as f:
        first_line = f.readline().strip().split("\t")

    if first_line == OLD_COLNAMES:
        old_df = pd.read_csv(args.input, sep="\t", header=0, dtype=str)
    else:
        old_df = pd.read_csv(args.input, sep="\t", header=None, names=OLD_COLNAMES, dtype=str)

    new_df = migrate(old_df)

    if args.output:
        new_df.to_csv(args.output, sep="\t", index=False)
        print(f"Migrated {len(new_df)} samples to {args.output}", file=sys.stderr)
    else:
        new_df.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
