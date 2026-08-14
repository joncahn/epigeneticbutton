"""Centralized sample-sheet parsing and metadata utilities.

Replaces the scattered parsing logic from the Snakefile (lines 62-193)
with a single importable module. All sample-sheet column names, assay
vocabularies, environment mappings, and helper functions live here.
"""

import io
import re
from collections import OrderedDict
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NEW_COLNAMES = [
    "Sample_ID", "Assay", "Genome", "Levels", "Replicate_ID",
    "Read_files", "Read_layout", "IP_target", "Control",
]

VALID_ASSAYS = [
    "ChIP_broad", "ChIP_narrow",
    "CUT_RUN_broad", "CUT_RUN_narrow",
    "CUT_TAG_broad", "CUT_TAG_narrow",
    "ATAC",
    "RNAseq", "RAMPAGE",
    "sRNA",
    "WGBS", "WGBS_nd", "PBAT", "EMseq", "dmC",
]

# Assays that pull down a target via antibody and call peaks against an
# (Input/WCE/IgG) control. They share the ChIP env, peak-type machinery,
# and IP_target/Control sample-sheet semantics.
IP_PEAK_ASSAYS = {
    "ChIP_broad", "ChIP_narrow",
    "CUT_RUN_broad", "CUT_RUN_narrow",
    "CUT_TAG_broad", "CUT_TAG_narrow",
}

ASSAY_TO_ENV = {
    "ChIP_broad": "ChIP",
    "ChIP_narrow": "ChIP",
    "CUT_RUN_broad": "ChIP",
    "CUT_RUN_narrow": "ChIP",
    "CUT_TAG_broad": "ChIP",
    "CUT_TAG_narrow": "ChIP",
    "ATAC": "ATAC",
    "RNAseq": "RNA",
    "RAMPAGE": "RNA",
    "sRNA": "sRNA",
    "WGBS": "mC",
    "WGBS_nd": "mC",
    "PBAT": "mC",
    "EMseq": "mC",
    "dmC": "mC",
}

ASSAY_TO_PEAKTYPE = {
    "ChIP_broad": "broad",
    "ChIP_narrow": "narrow",
    "CUT_RUN_broad": "broad",
    "CUT_RUN_narrow": "narrow",
    "CUT_TAG_broad": "broad",
    "CUT_TAG_narrow": "narrow",
    "ATAC": "narrow",
}

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _row_apply(df, func):
    """Like ``df.apply(func, axis=1)`` but safe on an empty DataFrame.

    ``DataFrame.apply(..., axis=1)`` returns an empty *DataFrame* (not a
    Series) when there are no rows, which cannot be assigned back to a single
    column. An empty sheet is a legitimate input (e.g. ``snakemake --unlock``
    parses the workflow without a real sample sheet), so guard against it.
    """
    if len(df) == 0:
        return pd.Series(dtype=object)
    return df.apply(func, axis=1)


# ---------------------------------------------------------------------------
# Levels helpers
# ---------------------------------------------------------------------------

def parse_levels(levels_str):
    """Parse a Levels string into an OrderedDict of factor:level pairs.

    >>> parse_levels("genotype:WT,tissue:root")
    OrderedDict([('genotype', 'WT'), ('tissue', 'root')])
    """
    result = OrderedDict()
    for pair in levels_str.split(","):
        pair = pair.strip()
        if ":" not in pair:
            raise ValueError(
                f"Invalid Levels entry '{pair}' — expected 'factor:level' format"
            )
        factor, level = pair.split(":", 1)
        result[factor.strip()] = level.strip()
    return result


def levels_to_label(levels_str):
    """Convert Levels string to a filesystem-safe label using level values.

    >>> levels_to_label("genotype:WT,tissue:root")
    'WT_root'
    """
    parsed = parse_levels(levels_str)
    return "_".join(parsed.values())


def levels_to_factors(levels_str):
    """Return the list of factor names from a Levels string.

    >>> levels_to_factors("genotype:WT,tissue:root")
    ['genotype', 'tissue']
    """
    return list(parse_levels(levels_str).keys())


# ---------------------------------------------------------------------------
# Read-file helpers
# ---------------------------------------------------------------------------

def parse_read_files(read_files_str, read_layout):
    """Parse the Read_files field into structured information.

    Returns (components, is_sra) where:
    - components: list of individual read-file entries (after splitting on '+')
      Each component is a string (SRA ID or file path, possibly comma-separated for PE)
    - is_sra: True if the entries are SRA accessions

    For SRA with '+' merge: "SRR111+SRR222" → (["SRR111", "SRR222"], True)
    For local PE: "/a/r1.fq,/a/r2.fq" → (["/a/r1.fq,/a/r2.fq"], False)
    For local SE: "/a/reads.fq" → (["/a/reads.fq"], False)
    """
    parts = [p.strip() for p in read_files_str.split("+")]
    # Check if any part looks like an SRA/ENA/DDBJ accession (SRR/ERR/DRR)
    is_sra = bool(re.match(r"^[SED]RR\d+$", parts[0]))
    return parts, is_sra


_FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _is_url(path):
    """Return True if the path looks like an HTTP(S) URL."""
    return path.startswith("http://") or path.startswith("https://")


def _url_basename(url):
    """Extract the filename from a URL, stripping query parameters."""
    import os
    path_part = url.split("?")[0].split("#")[0]
    return os.path.basename(path_part)


def get_seq_id_and_path(read_files_str, read_layout):
    """Bridge function: derive old-style seq_id and fastq_path from Read_files.

    Returns (seq_id, fastq_path) for backward compatibility with existing
    download rules that use these fields.

    Read_files formats handled:
    - SRA accession(s): "SRR111" or "SRR111+SRR222" → seq_id="SRR111,SRR222", path="SRA"
    - BAM file: "/path/to/file.bam" → seq_id="file", path="/path/to/file.bam"
    - BAM URL: "https://host/file.bam" → seq_id="file", path="https://host/file.bam"
    - bedMethyl: "/path/to/file.bed.gz" → seq_id="file", path="/path/to/file.bed.gz"
    - bedMethyl URL: "https://host/file.bed.gz" → seq_id="file", path="https://..."
    - Explicit FASTQ paths: "r1.fq.gz,r2.fq.gz" → seq_id="EXPLICIT", path="r1.fq.gz,r2.fq.gz"
    - FASTQ URL(s): "https://host/r.fq.gz" → seq_id="URL", path="https://host/r.fq.gz"
    """
    parts, is_sra = parse_read_files(read_files_str, read_layout)
    if is_sra:
        # SRA: seq_id is the comma-joined accessions, path is "SRA"
        return ",".join(parts), "SRA"
    else:
        import os
        first_file = parts[0].split(",")[0]  # first mate for PE
        # Strip query params for extension checks on URLs
        check_name = _url_basename(first_file) if _is_url(first_file) else first_file
        if check_name.endswith(".bam"):
            # BAM input (local path or URL): path IS the file/URL
            seq_id = os.path.splitext(os.path.basename(check_name))[0]
            return seq_id, first_file
        elif check_name.endswith(".bed.gz") or check_name.endswith(".bedmethyl.gz"):
            # bedMethyl input (local path or URL): path IS the file/URL
            base = os.path.basename(check_name)
            for suffix in (".bedmethyl.gz", ".bed.gz"):
                if base.endswith(suffix):
                    seq_id = base[:-len(suffix)]
                    break
            return seq_id, first_file
        elif any(check_name.endswith(ext) for ext in _FASTQ_EXTENSIONS):
            # Carry ALL '+'-components so the download rules can merge them
            # (each component is one SE file or a comma-separated R1,R2 pair).
            if _is_url(first_file):
                # FASTQ URL(s): use "URL" sentinel so download rules can dispatch
                return "URL", "+".join(parts)
            else:
                # Explicit FASTQ path(s): pass through the full Read_files string
                return "EXPLICIT", "+".join(parts)
        else:
            raise ValueError(
                f"Unrecognized Read_files format: '{read_files_str}'. "
                f"Expected SRA accession, .bam, .bed.gz, or explicit FASTQ "
                f"path ending in .fastq.gz/.fq.gz/.fastq/.fq"
            )


# ---------------------------------------------------------------------------
# Sample sheet I/O
# ---------------------------------------------------------------------------

def read_sample_sheet(filepath):
    """Read and validate a new-format sample sheet TSV.

    Full-line comments (lines whose first non-whitespace character is '#')
    are skipped — useful for parking samples without deleting their rows.
    A '#' *within* a field (e.g. a URL fragment or a path) is preserved.

    Returns a DataFrame sorted by [Genome, Assay, Levels, Sample_ID].
    """
    # Drop full-line comments before parsing. Note: we deliberately do NOT
    # use pandas' comment="#", which truncates every line at the first '#'
    # and would corrupt any Read_files URL/path containing one.
    with open(filepath) as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", header=0, dtype=str)
    df.columns = df.columns.str.strip()

    # Fill optional columns with empty strings
    for col in NEW_COLNAMES:
        if col not in df.columns:
            if col in ("IP_target", "Control"):
                df[col] = ""
            else:
                raise ValueError(f"Required column '{col}' missing from sample sheet")

    # Normalize NaN to empty string for optional columns
    df["IP_target"] = df["IP_target"].fillna("")
    df["Control"] = df["Control"].fillna("")

    df = explode_genomes(df)

    df = df.sort_values(
        by=["Genome", "Assay", "Levels", "Sample_ID"]
    ).reset_index(drop=True)

    return df


def parse_genomes(genome_str):
    """Split a Genome cell into its list of reference genomes.

    ``Genome`` accepts a comma-separated list so one sample can be mapped to
    several references (e.g. ``B73,W22``) — the reads are identical, only the
    alignment target differs. A single genome is just a one-element list.
    """
    return [g.strip() for g in str(genome_str or "").split(",") if g.strip()]


def explode_genomes(df):
    """One input row per sample becomes one internal row per (sample, genome).

    Downstream code — env mapping, analysis keys, every rule — assumes exactly
    one reference genome per row. Exploding here keeps that assumption intact
    while letting the sheet express multi-genome mapping on a single line, so
    the reads are declared (and fetched, and trimmed) only once.

    ``Sample_ID`` is deliberately NOT made unique per genome: it identifies the
    library, not the alignment. Post-alignment paths are disambiguated by the
    separate ``mapped_name`` token (see add_compat_columns).
    """
    if not len(df):
        return df
    if not df["Genome"].astype(str).str.contains(",").any():
        return df  # fast path: nothing to explode
    df = df.copy()
    df["Genome"] = df["Genome"].apply(parse_genomes)
    # Duplicates must be caught here: after the explode they are indistinguishable
    # from two legitimately separate rows, and would silently produce a doubled
    # DAG for the same (sample, genome).
    for _, row in df.iterrows():
        names = row["Genome"]
        if len(set(names)) != len(names):
            raise ValueError(
                f"Sample '{row['Sample_ID']}': Genome lists the same reference "
                f"more than once ({', '.join(names)})"
            )
    return df.explode("Genome", ignore_index=True)


# ---------------------------------------------------------------------------
# Analysis key / name building
# ---------------------------------------------------------------------------

def build_analysis_key(row):
    """Build an analysis-level grouping key from a DataFrame row.

    Returns a tuple: (Assay, levels_label, IP_target, Genome)
    """
    levels_label = levels_to_label(row["Levels"])
    ip_target = row.get("IP_target", "") or ""
    return (row["Assay"], levels_label, ip_target, row["Genome"])


def build_analysis_name(row):
    """Build an analysis-level name string from a DataFrame row.

    Components are joined with '__', omitting empty parts (e.g. blank IP_target).
    Examples:
      ChIP: "ChIP_broad__WT_root__H3K9me2__ColCEN"
      RNA:  "RNAseq__WT_root__ColCEN"
    """
    key = build_analysis_key(row)
    return "__".join(part for part in key if part)


def build_control_merge_key(row):
    """Build a merge-key for control samples that ignores Assay.

    Returns a tuple: (levels_label, IP_target, Genome). Used to group
    control-sample rows (Input/WCE/IgG) into one merged-replicate group
    when they represent the same biological control material, regardless
    of which Assay the rows happen to be labeled with. A single Input
    that serves as control for both broad and narrow ChIP IPs merges
    correctly as long as Levels, IP_target, and Genome match — even if
    rep1 is labeled ChIP_broad and rep2 ChIP_narrow.
    """
    levels_label = levels_to_label(row["Levels"])
    ip_target = row.get("IP_target", "") or ""
    return (levels_label, ip_target, row["Genome"])


def build_analysis_to_replicates(df):
    """Build a dict mapping analysis_key → list of Replicate_IDs.

    Only includes non-control samples.
    """
    non_control = df[~df["Sample_ID"].isin(identify_control_samples(df))]
    result = {}
    for _, row in non_control.iterrows():
        key = build_analysis_key(row)
        result.setdefault(key, []).append(row["Replicate_ID"])
    return result


def get_replicate_sample_ids(analysis_name, df):
    """Get the genome-qualified replicate names for an analysis group.

    Returns ``mapped_name`` values ('{Sample_ID}__{Genome}'), not bare
    Sample_IDs: every caller uses these to build post-alignment paths, where the
    reference is part of the file's identity. The *lookup key* is still a bare
    analysis name or control Sample_ID.


    The analysis_name can be:
    - An analysis name like "ChIP_broad__WT_cell__H3K9me2__Spombe"
    - A control Sample_ID like "WT_cell_Input_rep1" (for merging controls)

    Returns the list of Sample_IDs that belong to that analysis group.
    """
    controls = identify_control_samples(df)
    non_control = df[~df["Sample_ID"].isin(controls)]

    # Try matching non-control analysis names first
    result = []
    for _, row in non_control.iterrows():
        if build_analysis_name(row) == analysis_name:
            result.append(row.get("mapped_name", row["Sample_ID"]))
    if result:
        return result

    # If no match, check if analysis_name is a control Sample_ID
    # and find all controls in the same analysis group. Control merging
    # uses build_control_merge_key (Assay-agnostic) so a single biological
    # Input/IgG/WCE merges across all reps whose Levels/IP_target/Genome
    # match, even if individual rep rows are labeled with different Assay
    # values (e.g. one rep ChIP_broad, another ChIP_narrow).
    if analysis_name in set(controls):
        control_df = df[df["Sample_ID"].isin(controls)]
        match = control_df[control_df["Sample_ID"] == analysis_name]
        if not match.empty:
            ctrl_key = build_control_merge_key(match.iloc[0])
            for _, row in control_df.iterrows():
                if build_control_merge_key(row) == ctrl_key:
                    result.append(row.get("mapped_name", row["Sample_ID"]))

    return result


# ---------------------------------------------------------------------------
# Control handling
# ---------------------------------------------------------------------------

def identify_control_samples(df):
    """Return the set of Sample_IDs that are referenced as controls."""
    controls = set()
    for val in df["Control"]:
        if pd.notna(val) and str(val).strip():
            controls.add(str(val).strip())
    return controls


def match_sample_rows(df, name):
    """Rows for a sample named either bare or genome-qualified.

    Pre-alignment code passes a bare ``Sample_ID``; post-alignment code passes
    the ``mapped_name`` ('{Sample_ID}__{Genome}'). Accepting both here means the
    lookup helpers work on either side of the alignment boundary without every
    caller having to know which token it holds.
    """
    match = df.loc[df["Sample_ID"] == name]
    if match.empty and "mapped_name" in df.columns:
        match = df.loc[df["mapped_name"] == name]
    return match


def get_control_sample_id(sample_id, df):
    """Get the control for a given sample.

    The returned name matches the form of the input: asked with a
    genome-qualified name, the control comes back qualified with the SAME
    genome, because the caller is building a post-alignment path (the control's
    BAM) and a control aligned to a different reference would be wrong.
    """
    match = match_sample_rows(df, sample_id)
    if match.empty:
        return None
    row = match.iloc[0]
    control = str(row["Control"]).strip()
    if not control or control == "nan":
        return None
    # Qualify the control iff the caller used a qualified name.
    if "mapped_name" in df.columns and sample_id == row.get("mapped_name"):
        return f"{control}__{row['ref_genome']}"
    return control


# ---------------------------------------------------------------------------
# Analysis samples (non-control)
# ---------------------------------------------------------------------------

def get_analysis_samples(df):
    """Return a DataFrame of analysis-level samples (controls filtered out).

    Deduplicated by analysis key (one row per unique analysis group).
    When an analysis group has mixed PE/SE replicates, the analysis-level
    Read_layout is forced to SE and a warning is printed (PE reads are
    treated as SE for merged peak calling to avoid MACS2 errors).
    """
    import sys

    controls = identify_control_samples(df)
    non_control = df[~df["Sample_ID"].isin(controls)].copy()

    # Deduplicate by analysis key
    non_control["_analysis_key"] = _row_apply(
        non_control, lambda row: build_analysis_key(row)
    )

    # Detect mixed PE/SE within analysis groups before deduplication
    mixed_groups = []
    for key, group in non_control.groupby("_analysis_key"):
        layouts = group["Read_layout"].unique()
        if len(layouts) > 1:
            name = build_analysis_name(group.iloc[0])
            sids = group["Sample_ID"].tolist()
            mixed_groups.append((key, name, sids))

    deduped = non_control.drop_duplicates(subset=["_analysis_key"]).copy()

    # Force SE for mixed groups (PE reads treated as SE for merged analysis)
    if mixed_groups:
        for key, name, sids in mixed_groups:
            mask = deduped["_analysis_key"] == key
            deduped.loc[mask, "Read_layout"] = "SE"
            sid_list = ", ".join(sids)
            print(
                f"WARNING: Analysis group '{name}' has mixed PE/SE "
                f"replicates ({sid_list}). Merged peak calling will treat "
                f"all reads as SE. Per-replicate peak calls use each "
                f"sample's own layout.",
                file=sys.stderr,
            )

    deduped = deduped.drop(columns=["_analysis_key"])
    return deduped


# ---------------------------------------------------------------------------
# Field lookup
# ---------------------------------------------------------------------------

def get_sample_field(sample_id, df, field):
    """Look up a single field value by bare or genome-qualified sample name."""
    match = match_sample_rows(df, sample_id)
    if match.empty:
        return None
    return match[field].iloc[0]


# ---------------------------------------------------------------------------
# Environment and peak type
# ---------------------------------------------------------------------------

def get_env(assay):
    """Map an Assay value to its pipeline environment folder name."""
    env = ASSAY_TO_ENV.get(assay)
    if env is None:
        raise ValueError(f"Unknown assay '{assay}' — not in ASSAY_TO_ENV")
    return env


def get_peaktype(assay, config_override=None):
    """Get peak type (broad/narrow) for a given Assay.

    If config_override is provided (a dict mapping assay→peaktype),
    it takes precedence.
    """
    if config_override and assay in config_override:
        return config_override[assay]
    pt = ASSAY_TO_PEAKTYPE.get(assay)
    if pt is None:
        raise ValueError(
            f"No peak type defined for assay '{assay}'. "
            "Peak types are only defined for ChIP_broad, ChIP_narrow, "
            "CUT_RUN_broad, CUT_RUN_narrow, CUT_TAG_broad, "
            "CUT_TAG_narrow, and ATAC."
        )
    return pt


# ---------------------------------------------------------------------------
# Backward-compatibility bridge columns
# ---------------------------------------------------------------------------

def add_compat_columns(df):
    """Add backward-compatible columns derived from new-format fields.

    These columns bridge the gap during rule migration, allowing rules
    that still reference old column names to work with new-format data.

    Adds: data_type, ref_genome, replicate, paired, line, tissue,
          sample_type, extra_info, levels_label, env, sample_name
    """
    df = df.copy()

    # Direct mappings
    df["data_type"] = df["Assay"]
    df["ref_genome"] = df["Genome"]
    df["replicate"] = df["Replicate_ID"]
    df["paired"] = df["Read_layout"]

    # Levels → line (1st factor value) and tissue (2nd factor value)
    def _extract_level(levels_str, index):
        parsed = parse_levels(levels_str)
        vals = list(parsed.values())
        if index < len(vals):
            return vals[index]
        return ""

    df["line"] = df["Levels"].apply(lambda x: _extract_level(x, 0))
    df["tissue"] = df["Levels"].apply(lambda x: _extract_level(x, 1))

    # sample_type: IP_target for IP-with-peaks assays, else Assay
    def _derive_sample_type(row):
        if row["Assay"] in IP_PEAK_ASSAYS:
            return row["IP_target"] if row["IP_target"] else row["Assay"]
        return row["Assay"]

    df["sample_type"] = _row_apply(df, _derive_sample_type)

    # extra_info: IP_target or "N/A"
    df["extra_info"] = df["IP_target"].apply(
        lambda x: x if x else "N/A"
    )

    # levels_label
    df["levels_label"] = df["Levels"].apply(levels_to_label)

    # env
    df["env"] = df["Assay"].map(ASSAY_TO_ENV)

    # sample_name = Sample_ID (the key bridge)
    df["sample_name"] = df["Sample_ID"]

    # Post-alignment identity. Everything up to and including read trimming is
    # genome-independent and keeps the bare sample_name, so one download+trim
    # serves every reference. From alignment onward the reference is part of
    # what the file *is*, so those paths use mapped_name and two genomes can no
    # longer collide on one filename (issue #39). '__' is the reserved
    # delimiter — Sample_ID and Genome are both validated to exclude it, so
    # rsplit('__', 1) recovers the pair unambiguously.
    df["mapped_name"] = df["Sample_ID"].astype(str) + "__" + df["ref_genome"].astype(str)

    # Derive seq_id and fastq_path from Read_files
    def _derive_seq_id(row):
        seq_id, _ = get_seq_id_and_path(row["Read_files"], row["Read_layout"])
        return seq_id

    def _derive_fastq_path(row):
        _, fq_path = get_seq_id_and_path(row["Read_files"], row["Read_layout"])
        return fq_path

    df["seq_id"] = _row_apply(df, _derive_seq_id)
    df["fastq_path"] = _row_apply(df, _derive_fastq_path)

    return df
