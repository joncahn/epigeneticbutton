"""Centralized sample-sheet parsing and metadata utilities.

Replaces the scattered parsing logic from the Snakefile (lines 62-193)
with a single importable module. All sample-sheet column names, assay
vocabularies, environment mappings, and helper functions live here.
"""

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
            if _is_url(first_file):
                # FASTQ URL(s): use "URL" sentinel so download rules can dispatch
                return "URL", parts[0]
            else:
                # Explicit FASTQ path(s): pass through the full Read_files string
                return "EXPLICIT", parts[0]
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

    Lines starting with '#' (and any '# ...' tail of a non-comment line)
    are skipped as comments — useful for parking samples without deleting
    their rows.

    Returns a DataFrame sorted by [Genome, Assay, Levels, Sample_ID].
    """
    df = pd.read_csv(filepath, sep="\t", header=0, dtype=str, comment="#")
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

    df = df.sort_values(
        by=["Genome", "Assay", "Levels", "Sample_ID"]
    ).reset_index(drop=True)

    return df


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
    """Get list of replicate Sample_IDs for a given analysis group.

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
            result.append(row["Sample_ID"])
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
                    result.append(row["Sample_ID"])

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


def get_control_sample_id(sample_id, df):
    """Get the Control sample_id for a given sample.

    Returns the Sample_ID of the control, or None if no control.
    """
    match = df.loc[df["Sample_ID"] == sample_id]
    if match.empty:
        return None
    control = str(match["Control"].iloc[0]).strip()
    if not control:
        return None
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
    non_control["_analysis_key"] = non_control.apply(
        lambda row: build_analysis_key(row), axis=1
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
    """Look up a single field value for a sample by Sample_ID."""
    match = df.loc[df["Sample_ID"] == sample_id]
    if match.empty:
        # Also try sample_name column (backward compatibility bridge)
        if "sample_name" in df.columns:
            match = df.loc[df["sample_name"] == sample_id]
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

    df["sample_type"] = df.apply(_derive_sample_type, axis=1)

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

    # Derive seq_id and fastq_path from Read_files
    def _derive_seq_id(row):
        seq_id, _ = get_seq_id_and_path(row["Read_files"], row["Read_layout"])
        return seq_id

    def _derive_fastq_path(row):
        _, fq_path = get_seq_id_and_path(row["Read_files"], row["Read_layout"])
        return fq_path

    df["seq_id"] = df.apply(_derive_seq_id, axis=1)
    df["fastq_path"] = df.apply(_derive_fastq_path, axis=1)

    return df
