#!/usr/bin/env python3
"""Convert SEACR / epic2 peak-caller output to UCSC broadPeak / narrowPeak.

All input formats use 0-based half-open coordinates (BED convention),
so coordinates pass through unchanged.

broadPeak  (9 cols):  chrom start end name score strand signalValue pValue qValue
narrowPeak (10 cols): broadPeak + peak (summit offset relative to start)

Score columns follow UCSC track conventions: ``score`` is an integer in
[0, 1000]; ``pValue`` and ``qValue`` are -log10 transforms (-1 means "not
provided"); ``signalValue`` is a caller-specific magnitude (SEACR total
signal, epic2 log2 fold change).
"""
import argparse
import math
import re
import sys


def _safe_neg_log10(p):
    if p is None or p <= 0:
        return -1.0
    return -math.log10(p)


def _scale_score(value, max_seen):
    """Map a non-negative numeric value to a UCSC track score in [0, 1000]."""
    if max_seen <= 0:
        return 0
    return int(round(min(1000.0, 1000.0 * value / max_seen)))


def convert_seacr(in_path, out_path, peaktype):
    """SEACR .stringent.bed / .relaxed.bed → broadPeak or narrowPeak.

    Input columns (6, tab-separated):
        chrom, start, end, total_signal, max_signal, max_signal_region

    where ``max_signal_region`` is "chrom:start-end" naming the contiguous
    sub-interval of maximum coverage. We convert that to a peak-summit
    offset (``narrowPeak`` only).
    """
    rows = []
    with open(in_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 6:
                continue
            chrom, start_s, end_s, total_s, max_s, region = cols[:6]
            try:
                start = int(start_s)
                end = int(end_s)
                total = float(total_s)
                mx = float(max_s)
            except ValueError:
                continue
            m = re.match(r"^[^:]+:(\d+)-(\d+)$", region)
            if m:
                ss, se = int(m.group(1)), int(m.group(2))
                summit_offset = (ss + se) // 2 - start
                if summit_offset < 0 or summit_offset >= (end - start):
                    summit_offset = -1
            else:
                summit_offset = -1
            rows.append((chrom, start, end, total, mx, summit_offset))

    max_total = max((r[3] for r in rows), default=0.0)
    with open(out_path, "w") as out:
        for i, (chrom, start, end, total, _mx, summit) in enumerate(rows, start=1):
            score = _scale_score(total, max_total)
            if peaktype == "narrow":
                out.write(
                    f"{chrom}\t{start}\t{end}\tpeak_{i}\t{score}\t.\t"
                    f"{total:.6g}\t-1\t-1\t{summit}\n"
                )
            else:
                out.write(
                    f"{chrom}\t{start}\t{end}\tpeak_{i}\t{score}\t.\t"
                    f"{total:.6g}\t-1\t-1\n"
                )


def convert_epic2(in_path, out_path):
    """epic2 default output → broadPeak.

    Input columns (epic2 ≥ 0.0.40):
        Chromosome Start End PValue Score Strand ChIPCount InputCount FDR log2FoldChange
    """
    rows = []
    with open(in_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            # Skip header row if present
            if cols[0].lower() in ("chromosome", "chrom"):
                continue
            try:
                chrom = cols[0]
                start = int(cols[1])
                end = int(cols[2])
                pvalue = float(cols[3])
                score_raw = float(cols[4])
                fdr = float(cols[8])
            except (ValueError, IndexError):
                continue
            lfc = 0.0
            if len(cols) > 9:
                try:
                    lfc = float(cols[9])
                except ValueError:
                    lfc = 0.0
            score = max(0, min(1000, int(round(score_raw))))
            rows.append((chrom, start, end, score, lfc, pvalue, fdr))

    with open(out_path, "w") as out:
        for i, (chrom, start, end, score, lfc, pvalue, fdr) in enumerate(rows, start=1):
            out.write(
                f"{chrom}\t{start}\t{end}\tpeak_{i}\t{score}\t.\t"
                f"{lfc:.6g}\t{_safe_neg_log10(pvalue):.6g}\t"
                f"{_safe_neg_log10(fdr):.6g}\n"
            )


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--caller", choices=["seacr", "epic2"], required=True)
    p.add_argument("--peaktype", choices=["broad", "narrow"], required=True)
    p.add_argument("input")
    p.add_argument("output")
    args = p.parse_args()

    if args.caller == "seacr":
        convert_seacr(args.input, args.output, args.peaktype)
    elif args.caller == "epic2":
        if args.peaktype != "broad":
            sys.exit("epic2 only outputs broadPeak (use seacr or macs2 for narrow)")
        convert_epic2(args.input, args.output)


if __name__ == "__main__":
    main()
