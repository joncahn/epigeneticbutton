#!/usr/bin/env bash
# ena_download.sh — Download FASTQ files from ENA for a given SRA/ERA/DRA accession.
#
# Usage: bash ena_download.sh <accession> <output_dir> <SE|PE>
#
# Arguments:
#   accession   SRR/ERR/DRR accession (e.g. SRR20678305)
#   output_dir  Directory to write downloaded .fastq.gz files into
#   layout      SE or PE — determines which files to download
#
# Downloads:
#   SE: {output_dir}/{accession}.fastq.gz
#   PE: {output_dir}/{accession}_1.fastq.gz + {output_dir}/{accession}_2.fastq.gz
#
# Exit codes:
#   0  All files downloaded successfully
#   1  Usage error or download failure

set -euo pipefail

if [[ $# -ne 3 ]]; then
    printf "Usage: %s <accession> <output_dir> <SE|PE>\n" "$0" >&2
    exit 1
fi

accession="$1"
output_dir="$2"
layout="$3"

# --- Construct ENA FTP URL ---
# URL pattern: https://ftp.sra.ebi.ac.uk/vol1/fastq/{first6}/{intermediate}/{accession}/
# The intermediate directory depends on the number of digits after the prefix letters.

# Extract the numeric suffix (strip leading letters)
prefix="${accession%%[0-9]*}"
digits="${accession#"$prefix"}"
ndigits="${#digits}"

first6="${accession:0:6}"

case "$ndigits" in
    6) intermediate="" ;;
    7) intermediate="00${digits: -1}" ;;
    8) intermediate="0${digits: -2}" ;;
    9) intermediate="${digits: -3}" ;;
    *)
        printf "ERROR: Unexpected accession format '%s' (%d digits)\n" "$accession" "$ndigits" >&2
        exit 1
        ;;
esac

if [[ -n "$intermediate" ]]; then
    base_url="https://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${intermediate}/${accession}"
else
    base_url="https://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${accession}"
fi

# --- Download files ---
curl_opts=(-fsSL --retry 2 --connect-timeout 30)

download_file() {
    local url="$1"
    local dest="$2"
    printf "ENA: downloading %s\n" "$url"
    if ! curl "${curl_opts[@]}" -o "$dest" "$url"; then
        printf "ENA: failed to download %s\n" "$url" >&2
        rm -f "$dest"
        return 1
    fi
    if ! gzip -t "$dest" 2>/dev/null; then
        printf "ENA: gzip integrity check failed for %s\n" "$dest" >&2
        rm -f "$dest"
        return 1
    fi
}

mkdir -p "$output_dir"

if [[ "$layout" == "PE" ]]; then
    download_file "${base_url}/${accession}_1.fastq.gz" "${output_dir}/${accession}_1.fastq.gz" || exit 1
    download_file "${base_url}/${accession}_2.fastq.gz" "${output_dir}/${accession}_2.fastq.gz" || exit 1
elif [[ "$layout" == "SE" ]]; then
    download_file "${base_url}/${accession}.fastq.gz" "${output_dir}/${accession}.fastq.gz" || exit 1
else
    printf "ERROR: layout must be SE or PE, got '%s'\n" "$layout" >&2
    exit 1
fi

printf "ENA: successfully downloaded %s (%s)\n" "$accession" "$layout"
