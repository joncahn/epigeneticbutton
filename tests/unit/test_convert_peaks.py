"""Unit tests for SEACR/epic2 → broadPeak/narrowPeak conversion.

Coordinate sanity is the main thing under test: SEACR's 6-col BED, epic2's
default 10-col TSV, and the resulting broadPeak/narrowPeak output all use
0-based half-open coordinates, so a row's start/end should pass through
unchanged. The narrowPeak summit offset must be derived from SEACR's
``max_signal_region`` (a chrom:start-end string) and must fall inside
[0, end-start).
"""
import os
import math
import tempfile

import pytest

from workflow.scripts.convert_peaks import convert_seacr, convert_epic2


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write(content):
    fd, path = tempfile.mkstemp(suffix=".tsv")
    os.close(fd)
    with open(path, "w") as fh:
        fh.write(content)
    return path


def _read_rows(path):
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


# ---------------------------------------------------------------------------
# SEACR
# ---------------------------------------------------------------------------

class TestConvertSEACR:

    def _fixture(self):
        # chrom start end total_signal max_signal max_signal_region
        return _write(
            "chr1\t100\t300\t12.5\t3.0\tchr1:180-220\n"
            "chr1\t500\t900\t40.0\t5.5\tchr1:700-750\n"
        )

    def test_narrow_preserves_coords_and_columns(self):
        in_path = self._fixture()
        out_path = in_path + ".narrowPeak"
        convert_seacr(in_path, out_path, "narrow")
        rows = _read_rows(out_path)
        assert len(rows) == 2
        # narrowPeak = 10 columns
        assert all(len(r) == 10 for r in rows), rows
        # Coords pass through (0-based half-open)
        assert (rows[0][0], rows[0][1], rows[0][2]) == ("chr1", "100", "300")
        assert (rows[1][0], rows[1][1], rows[1][2]) == ("chr1", "500", "900")
        # Score is a non-negative integer in [0, 1000]
        for r in rows:
            s = int(r[4])
            assert 0 <= s <= 1000
        # Strand placeholder
        assert rows[0][5] == "."
        # signalValue is the SEACR total (col 4 in input)
        assert float(rows[0][6]) == pytest.approx(12.5)
        assert float(rows[1][6]) == pytest.approx(40.0)
        # pValue / qValue placeholders are -1 (SEACR doesn't compute)
        assert rows[0][7] == "-1" and rows[0][8] == "-1"
        # Summit offset = midpoint(180, 220) - start(100) = 200 - 100 = 100
        assert int(rows[0][9]) == 100
        # Second row: midpoint(700,750)=725; 725-500 = 225
        assert int(rows[1][9]) == 225

    def test_broad_drops_summit_column(self):
        in_path = self._fixture()
        out_path = in_path + ".broadPeak"
        convert_seacr(in_path, out_path, "broad")
        rows = _read_rows(out_path)
        # broadPeak = 9 columns (no peak/summit)
        assert all(len(r) == 9 for r in rows), rows
        assert (rows[0][0], rows[0][1], rows[0][2]) == ("chr1", "100", "300")

    def test_handles_missing_or_malformed_region(self):
        # If max_signal_region isn't parseable, summit defaults to -1
        path = _write("chr2\t10\t90\t5.0\t1.5\tNOT_A_REGION\n")
        out = path + ".narrowPeak"
        convert_seacr(path, out, "narrow")
        rows = _read_rows(out)
        assert int(rows[0][9]) == -1


# ---------------------------------------------------------------------------
# epic2
# ---------------------------------------------------------------------------

class TestConvertEpic2:

    def _fixture(self):
        return _write(
            "Chromosome\tStart\tEnd\tPValue\tScore\tStrand\tChIPCount\tInputCount\tFDR\tlog2FoldChange\n"
            "chr1\t100\t500\t1e-08\t450\t.\t250\t40\t1e-06\t2.5\n"
            "chr2\t1000\t2000\t0.001\t150\t.\t60\t20\t0.01\t1.2\n"
        )

    def test_basic_conversion_and_minus_log10(self):
        in_path = self._fixture()
        out_path = in_path + ".broadPeak"
        convert_epic2(in_path, out_path)
        rows = _read_rows(out_path)
        assert len(rows) == 2
        # broadPeak = 9 columns
        assert all(len(r) == 9 for r in rows), rows
        # Coords pass through (0-based half-open)
        assert (rows[0][0], rows[0][1], rows[0][2]) == ("chr1", "100", "500")
        assert (rows[1][0], rows[1][1], rows[1][2]) == ("chr2", "1000", "2000")
        # signalValue is log2FoldChange (col 10 in input)
        assert float(rows[0][6]) == pytest.approx(2.5)
        # pValue is -log10(1e-08) = 8
        assert float(rows[0][7]) == pytest.approx(8.0, abs=0.01)
        # qValue is -log10(1e-06) = 6
        assert float(rows[0][8]) == pytest.approx(6.0, abs=0.01)
        # Score clipped to [0, 1000]
        assert 0 <= int(rows[0][4]) <= 1000

    def test_skips_header_row(self):
        # Header row (Chromosome\tStart\t...) should not appear in the output
        in_path = self._fixture()
        out_path = in_path + ".broadPeak"
        convert_epic2(in_path, out_path)
        rows = _read_rows(out_path)
        assert all(r[0] != "Chromosome" for r in rows)

    def test_zero_or_missing_pvalue_uses_minus_one(self):
        path = _write("chr1\t10\t20\t0\t5\t.\t100\t10\t0\t1.0\n")
        out = path + ".broadPeak"
        convert_epic2(path, out)
        rows = _read_rows(out)
        # log10(0) is undefined → sentinel -1
        assert float(rows[0][7]) == -1.0
        assert float(rows[0][8]) == -1.0
