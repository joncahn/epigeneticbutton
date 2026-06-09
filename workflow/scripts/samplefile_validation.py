"""Validation for new-format sample sheets.

Implements the rules defined in dev/docs/sample-sheet-spec.md.
"""

import os
import re
from scripts.sample_sheet import VALID_ASSAYS, ASSAY_TO_ENV, IP_PEAK_ASSAYS

# Characters that are unsafe for filesystem use in Sample_ID
_UNSAFE_CHARS = re.compile(r'[/\\\s\'\";&|<>$`!{}()\[\]?*~#]')
_DOUBLE_UNDERSCORE = re.compile(r'__')
_SRA_REGEX = re.compile(r'^[SDE]RR\d+$')


def _merge_component_kind(comp):
    """Classify one '+'-merge component for merge-compatibility checks.

    Returns 'sra', 'fastq', 'bam', 'bedmethyl', or 'other'. SRA accessions
    and FASTQ files can be '+'-merged (concatenated after download); BAM
    (needs samtools merge) and bedMethyl (counts must be summed) cannot.
    """
    first = comp.split(",")[0].strip()
    base = first.split("?")[0].split("#")[0].lower()
    if _SRA_REGEX.match(first):
        return "sra"
    if base.endswith(".bam"):
        return "bam"
    if base.endswith(".bed.gz") or base.endswith(".bedmethyl.gz"):
        return "bedmethyl"
    if any(base.endswith(e) for e in (".fastq.gz", ".fq.gz", ".fastq", ".fq")):
        return "fastq"
    return "other"

# Assays that may declare a Control sample. RAMPAGE is normalized against
# an RNA-seq Control rather than an IP-style Input/IgG, but the Control
# field semantics (Sample_ID reference, no chaining) are the same.
_CONTROL_ASSAYS = IP_PEAK_ASSAYS | {"RAMPAGE"}


def _is_url(path):
    """Return True if the path looks like an HTTP(S) URL."""
    return path.startswith("http://") or path.startswith("https://")


def _is_local_path(token):
    """True if a Read_files component entry is a local filesystem path.

    Excludes empty strings, SRA-style accessions, and HTTP(S) URLs.
    """
    if not token:
        return False
    if _SRA_REGEX.match(token):
        return False
    if _is_url(token):
        return False
    return True


def check_table(tab, check_paths=True):
    """Validate a new-format sample sheet DataFrame.

    When ``check_paths`` is True (default), local Read_files paths are checked
    for existence on disk. Pass ``check_paths=False`` for dry-run-style
    validation where the inputs may not yet be staged.

    Raises ValueError with all collected error messages if validation fails.
    """
    errors = []
    warnings = []

    # --- Sample_ID: required, unique, filesystem-safe ---
    sample_ids = set()
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        sid = str(row.get("Sample_ID", "")).strip()
        if not sid or sid == "nan":
            errors.append(f"[X] Row #{i}: Sample_ID is required")
            continue
        if _UNSAFE_CHARS.search(sid):
            errors.append(f"[X] Row #{i} '{sid}': Sample_ID contains unsafe characters")
        if _DOUBLE_UNDERSCORE.search(sid):
            errors.append(f"[X] Row #{i} '{sid}': Sample_ID must not contain '__'")
        if sid in sample_ids:
            errors.append(f"[X] Row #{i} '{sid}': duplicate Sample_ID")
        sample_ids.add(sid)

    # --- Assay: controlled vocabulary ---
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        assay = str(row.get("Assay", "")).strip()
        if assay not in VALID_ASSAYS:
            errors.append(
                f"[X] Row #{i} '{row.get('Sample_ID', '')}': "
                f"Assay '{assay}' not in {VALID_ASSAYS}"
            )

    # --- Genome: required ---
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        genome = str(row.get("Genome", "")).strip()
        if not genome or genome == "nan":
            errors.append(f"[X] Row #{i} '{row.get('Sample_ID', '')}': Genome is required")

    # --- Levels: required, consistent factors ---
    all_factors = []
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        levels = str(row.get("Levels", "")).strip()
        if not levels or levels == "nan":
            errors.append(f"[X] Row #{i} '{row.get('Sample_ID', '')}': Levels is required")
            continue
        pairs = levels.split(",")
        factor_names = []
        for pair in pairs:
            pair = pair.strip()
            if ":" not in pair:
                errors.append(
                    f"[X] Row #{i} '{row.get('Sample_ID', '')}': "
                    f"Levels entry '{pair}' must be 'factor:level' format"
                )
            else:
                factor, level = pair.split(":", 1)
                if not factor.strip() or not level.strip():
                    errors.append(
                        f"[X] Row #{i} '{row.get('Sample_ID', '')}': "
                        f"empty factor name or level in Levels"
                    )
                factor_names.append(factor.strip())
        # Store the sheet row number alongside the factors: rows with blank
        # Levels are skipped above, so a positional index would misreport
        # which row mismatched.
        all_factors.append((i, factor_names))

    # Check consistent factor count and names
    if all_factors:
        ref_row, ref_factors = all_factors[0]
        for row_num, factors in all_factors[1:]:
            if len(factors) != len(ref_factors):
                errors.append(
                    f"[X] Row #{row_num}: Levels has {len(factors)} factors, "
                    f"expected {len(ref_factors)} (matching row {ref_row})"
                )
            elif factors != ref_factors:
                errors.append(
                    f"[X] Row #{row_num}: Levels factor names {factors} "
                    f"don't match row {ref_row} {ref_factors}"
                )

    # --- Replicate_ID: required ---
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        rep = str(row.get("Replicate_ID", "")).strip()
        if not rep or rep == "nan":
            errors.append(f"[X] Row #{i} '{row.get('Sample_ID', '')}': Replicate_ID is required")

    # --- Read_layout: SE or PE ---
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        layout = str(row.get("Read_layout", "")).strip()
        if layout not in ("SE", "PE"):
            errors.append(
                f"[X] Row #{i} '{row.get('Sample_ID', '')}': "
                f"Read_layout must be 'SE' or 'PE', got '{layout}'"
            )

    # --- Read_files: required, SRA or paths ---
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        read_files = str(row.get("Read_files", "")).strip()
        layout = str(row.get("Read_layout", "")).strip()
        sid = str(row.get("Sample_ID", "")).strip()
        if not read_files or read_files == "nan":
            errors.append(f"[X] Row #{i} '{sid}': Read_files is required")
            continue

        components = [c.strip() for c in read_files.split("+")]
        if len(components) > 1:
            # '+'-merge concatenates components after download. Supported for
            # SRA accessions and FASTQ files; BAM (needs samtools merge) and
            # bedMethyl (counts must be summed, not concatenated) are not, and
            # components of different types may not be mixed.
            kinds = {_merge_component_kind(c) for c in components}
            unmergeable = kinds & {"bam", "bedmethyl", "other"}
            if len(kinds) > 1:
                errors.append(
                    f"[X] Row #{i} '{sid}': '+'-merge components must all be the "
                    f"same type (got {sorted(kinds)})"
                )
            elif unmergeable:
                bad = sorted(unmergeable)[0]
                errors.append(
                    f"[X] Row #{i} '{sid}': '+'-merge is not supported for {bad} "
                    f"inputs (only SRA accessions and FASTQ files); merge upstream instead"
                )
        for comp in components:
            # Each component is either an SRA ID or file path(s)
            files_in_comp = [f.strip() for f in comp.split(",")]
            first = files_in_comp[0]
            if _SRA_REGEX.match(first):
                # SRA: all parts should be SRA IDs
                for f in files_in_comp:
                    if not _SRA_REGEX.match(f):
                        errors.append(
                            f"[X] Row #{i} '{sid}': mixed SRA/path in Read_files component"
                        )
            else:
                # Local paths: validate PE has comma-separated pair
                if layout == "PE" and len(files_in_comp) == 1 and not first.endswith(".bam"):
                    errors.append(
                        f"[X] Row #{i} '{sid}': Read_layout is PE but Read_files "
                        f"has only one path (expected comma-separated pair)"
                    )
                if layout == "SE" and len(files_in_comp) > 1:
                    errors.append(
                        f"[X] Row #{i} '{sid}': Read_layout is SE but Read_files "
                        f"has multiple comma-separated paths"
                    )

    # --- Read_files: cross-row duplicate check ---
    seen_inputs = {}  # path/accession -> Sample_ID
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        read_files = str(row.get("Read_files", "")).strip()
        sid = str(row.get("Sample_ID", "")).strip()
        if not read_files or read_files == "nan":
            continue
        components = [c.strip() for c in read_files.split("+")]
        for comp in components:
            files_in_comp = [f.strip() for f in comp.split(",")]
            for f in files_in_comp:
                if not f:
                    continue
                if f in seen_inputs:
                    errors.append(
                        f"[X] Row #{i} '{sid}': Read_files entry '{f}' is also "
                        f"used by '{seen_inputs[f]}'"
                    )
                else:
                    seen_inputs[f] = sid

    # --- Read_files: local path existence check ---
    # For SRA accessions and HTTP(S) URLs we can't (or won't) probe ahead of
    # time; for local paths, fail fast if they don't exist on disk.
    if check_paths:
        for i, (_, row) in enumerate(tab.iterrows(), start=1):
            read_files = str(row.get("Read_files", "")).strip()
            sid = str(row.get("Sample_ID", "")).strip()
            if not read_files or read_files == "nan":
                continue
            for comp in (c.strip() for c in read_files.split("+")):
                for f in (x.strip() for x in comp.split(",")):
                    if _is_local_path(f) and not os.path.exists(f):
                        errors.append(
                            f"[X] Row #{i} '{sid}': Read_files path '{f}' does not exist"
                        )

    # --- IP_target: required for ChIP, blank for others ---
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        assay = str(row.get("Assay", "")).strip()
        ip_target = str(row.get("IP_target", "")).strip()
        sid = str(row.get("Sample_ID", "")).strip()
        if ip_target == "nan":
            ip_target = ""
        if assay in IP_PEAK_ASSAYS:
            if not ip_target:
                errors.append(
                    f"[X] Row #{i} '{sid}': IP_target is required for {assay}"
                )
        else:
            if ip_target:
                errors.append(
                    f"[X] Row #{i} '{sid}': IP_target must be blank for {assay}"
                )

    # --- Control: reference validation ---
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        assay = str(row.get("Assay", "")).strip()
        control = str(row.get("Control", "")).strip()
        sid = str(row.get("Sample_ID", "")).strip()
        if control == "nan":
            control = ""
        if control:
            if assay not in _CONTROL_ASSAYS:
                errors.append(
                    f"[X] Row #{i} '{sid}': Control is not allowed for {assay}"
                )
            if control not in sample_ids:
                errors.append(
                    f"[X] Row #{i} '{sid}': Control '{control}' does not match "
                    f"any Sample_ID in the sheet"
                )
            # No chaining
            ctrl_row = tab[tab["Sample_ID"] == control]
            if not ctrl_row.empty:
                ctrl_ctrl = str(ctrl_row["Control"].iloc[0]).strip()
                if ctrl_ctrl and ctrl_ctrl != "nan":
                    errors.append(
                        f"[X] Row #{i} '{sid}': Control '{control}' itself has "
                        f"a Control (chaining not allowed)"
                    )

    # --- Print warnings ---
    for w in warnings:
        print(w)

    if errors:
        full_message = "\n".join(errors)
        raise ValueError(
            f"[X] Validation failed — please fix the errors below in your "
            f"samplefile and rerun.\n{full_message}\n\n"
        )


# Envs that previously required genomesize — now auto-computed from FASTA
# (kept for reference; genomesize and star_index are no longer validated as required)


def check_genome_config(tab, config, check_paths=True):
    """Validate that config['genomes'] has required fields for all genomes in the sample sheet.

    When ``check_paths`` is True (default), local genome-config file paths
    (``fasta_file``, ``gff_file``, ``gtf_file``, ``te_file``,
    ``structural_rna_fafile``, ``gaf_file``, ``gene_info_file``) are checked
    for existence on disk. URLs and the ``<auto>`` sentinel are skipped.

    Raises ValueError with all collected error messages if validation fails.
    """
    errors = []
    warnings = []
    genomes_cfg = config.get("genomes", {})

    # Collect which envs are used per genome
    genome_envs = {}
    for _, row in tab.iterrows():
        genome = str(row.get("Genome", "")).strip()
        assay = str(row.get("Assay", "")).strip()
        env = ASSAY_TO_ENV.get(assay)
        if genome and genome != "nan":
            genome_envs.setdefault(genome, set())
            if env:
                genome_envs[genome].add(env)

    for genome, envs in genome_envs.items():
        if genome not in genomes_cfg:
            errors.append(
                f"[X] Genome '{genome}': no entry in config['genomes']. "
                f"Add it to the genomes: section of your options file."
            )
            continue

        gcfg = genomes_cfg[genome]

        # Always required
        for field in ("fasta_file", "gff_file"):
            if field not in gcfg:
                errors.append(f"[X] Genome '{genome}': missing required field '{field}'")

        # Local path existence checks (URLs and "<auto>" sentinels are skipped).
        # gtf_file is optional and defaults to "<auto>" (derived from gff_file).
        # te_file is optional. fasta_file and gff_file are required above.
        if check_paths:
            for field in ("fasta_file", "gff_file", "gtf_file", "te_file"):
                val = gcfg.get(field)
                if val and val != "<auto>" and not _is_url(val) and not os.path.exists(val):
                    errors.append(
                        f"[X] Genome '{genome}': {field} '{val}' does not exist"
                    )

        # genomesize and star_index are auto-computed from the reference FASTA;
        # user-provided values in the options file are optional overrides

        # structural_rna_fafile is optional: auto-derived via Infernal when absent or "<auto>"
        if check_paths and "sRNA" in envs and config.get("structural_rna_depletion", True):
            srna_fa = gcfg.get("structural_rna_fafile", "<auto>")
            if (srna_fa and srna_fa != "<auto>"
                    and not _is_url(srna_fa)
                    and not os.path.exists(srna_fa)):
                errors.append(
                    f"[X] Genome '{genome}': structural_rna_fafile '{srna_fa}' "
                    f"does not exist"
                )

        # GO fields when GO: true
        if config.get("GO", False) and "RNA" in envs:
            for field in ("genus", "species", "gaf_file", "gene_info_file"):
                if field not in gcfg:
                    errors.append(
                        f"[X] Genome '{genome}': missing '{field}' (required when GO: true)"
                    )
            # Existence check for the GO annotation files
            if check_paths:
                for field in ("gaf_file", "gene_info_file"):
                    val = gcfg.get(field)
                    if val and not _is_url(val) and not os.path.exists(val):
                        errors.append(
                            f"[X] Genome '{genome}': {field} '{val}' does not exist"
                        )

    # motif_ref_genome must reference a valid genome when motifs are enabled
    if config.get("motifs", False):
        motif_genome = config.get("motif_ref_genome", "")
        if motif_genome and motif_genome not in genomes_cfg:
            errors.append(
                f"[X] motif_ref_genome '{motif_genome}' not found in config['genomes']"
            )

    # --- Print warnings ---
    for w in warnings:
        print(w)

    if errors:
        full_message = "\n".join(errors)
        raise ValueError(
            f"[X] Genome config validation failed — please fix the errors below "
            f"in your options file and rerun.\n{full_message}\n\n"
        )


# ---------------------------------------------------------------------------
# Extra output target file format validation (Snakefile-driven analyses)
# ---------------------------------------------------------------------------

def _open_maybe_gzip(path):
    """Open a text file, transparently decompressing if it's gzipped."""
    import gzip
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _iter_data_lines(path):
    """Yield (line_number, columns) for non-comment, non-blank rows in a TSV.

    The first non-blank line is yielded even if it looks like a header — the
    caller decides whether to treat it as one based on column-1 content.
    """
    with _open_maybe_gzip(path) as fh:
        for i, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            yield i, line.split("\t")


def _looks_like_header(cols):
    """Heuristic: a BED first row is a header if cols[1]/cols[2] aren't ints."""
    if len(cols) < 3:
        return False
    try:
        int(cols[1]); int(cols[2])
        return False
    except (ValueError, IndexError):
        return True


def _validate_browser_target_file(path, errors):
    """Validate the browser_target_file's extended-BED schema.

    Each data row must be:
        chrom \t start \t end \t label \t binsize [\t htstart \t htwidth]

    where label is a non-empty string not beginning with '-' (so deeptools
    doesn't mis-parse it as a flag), binsize is an integer >= 1, and the
    optional htstart/htwidth columns are comma-separated coordinates and
    widths. The first row may be a header (column 2/3 non-integer).
    """
    saw_data = False
    for lineno, cols in _iter_data_lines(path):
        if not saw_data and _looks_like_header(cols):
            saw_data = True
            continue
        saw_data = True
        if len(cols) < 5:
            errors.append(
                f"[X] browser_target_file '{path}' line {lineno}: expected "
                f"at least 5 tab-separated columns "
                f"(chrom, start, end, label, binsize), got {len(cols)}"
            )
            continue
        chrom, start_s, end_s, label, binsize_s = cols[:5]
        try:
            start, end = int(start_s), int(end_s)
            if start < 0 or end <= start:
                errors.append(
                    f"[X] browser_target_file '{path}' line {lineno}: "
                    f"invalid coordinates start={start_s} end={end_s}"
                )
        except ValueError:
            errors.append(
                f"[X] browser_target_file '{path}' line {lineno}: "
                f"start/end not integers (start={start_s!r}, end={end_s!r})"
            )
        if not label or label.startswith("-"):
            errors.append(
                f"[X] browser_target_file '{path}' line {lineno}: label "
                f"{label!r} must be non-empty and must not start with '-' "
                f"(would be misinterpreted as a CLI flag downstream)"
            )
        try:
            bs = int(binsize_s)
            if bs < 1:
                errors.append(
                    f"[X] browser_target_file '{path}' line {lineno}: "
                    f"binsize must be >= 1, got {bs}"
                )
        except ValueError:
            errors.append(
                f"[X] browser_target_file '{path}' line {lineno}: binsize "
                f"{binsize_s!r} is not an integer"
            )
        # Optional htstart/htwidth: pair-or-neither, comma-separated ints
        ht_start = cols[5] if len(cols) > 5 else ""
        ht_width = cols[6] if len(cols) > 6 else ""
        if (ht_start and not ht_width) or (ht_width and not ht_start):
            errors.append(
                f"[X] browser_target_file '{path}' line {lineno}: htstart "
                f"and htwidth must both be present or both absent"
            )
        for label_, val in (("htstart", ht_start), ("htwidth", ht_width)):
            if not val:
                continue
            for tok in val.split(","):
                try:
                    v = int(tok)
                    if v < 0:
                        raise ValueError
                except ValueError:
                    errors.append(
                        f"[X] browser_target_file '{path}' line {lineno}: "
                        f"{label_} entry {tok!r} is not a non-negative integer"
                    )
    if not saw_data:
        errors.append(
            f"[X] browser_target_file '{path}' has no data rows"
        )


def _validate_bed_target_file(path, key, errors):
    """Validate a BED-format target file (chrom, start, end[, name[, score, strand]])."""
    saw_data = False
    for lineno, cols in _iter_data_lines(path):
        if not saw_data and _looks_like_header(cols):
            saw_data = True
            continue
        saw_data = True
        if len(cols) < 3:
            errors.append(
                f"[X] {key} '{path}' line {lineno}: expected ≥3 "
                f"tab-separated columns, got {len(cols)}"
            )
            continue
        try:
            start, end = int(cols[1]), int(cols[2])
            if start < 0 or end <= start:
                errors.append(
                    f"[X] {key} '{path}' line {lineno}: invalid "
                    f"coordinates start={cols[1]} end={cols[2]}"
                )
        except ValueError:
            errors.append(
                f"[X] {key} '{path}' line {lineno}: start/end not integers"
            )
    if not saw_data:
        errors.append(f"[X] {key} '{path}' has no data rows")


def _validate_geneid_target_file(path, key, errors):
    """Validate a TSV gene-id list (col 1 = gene id, optional col 2 = label)."""
    saw_data = False
    for lineno, cols in _iter_data_lines(path):
        saw_data = True
        if not cols[0].strip():
            errors.append(
                f"[X] {key} '{path}' line {lineno}: first column "
                f"(gene ID) is empty"
            )
    if not saw_data:
        errors.append(f"[X] {key} '{path}' has no data rows")


def check_extra_output_files(config, check_paths=True):
    """Validate optional 'extra output' target files referenced in the options.

    The pipeline supports several optional inputs that drive add-on outputs
    (motif scans, gene-expression plots, GO enrichment, sRNA cluster
    analysis, heatmaps, browser plots). Each is gated by an analysis flag
    in the options file; we only validate when the relevant analysis is
    enabled and the configured path exists. Missing-file errors are
    surfaced too, so users learn at startup rather than mid-run that a
    referenced file isn't where they said it was.

    Pass ``check_paths=False`` to skip on-disk existence (mirrors the
    other validators' behavior).
    """
    errors = []

    full_analysis = config.get("full_analysis", False)
    motifs_on = config.get("motifs", False)
    go_on = config.get("GO", False)

    # (key, gate, validator)
    checks = [
        ("browser_target_file", full_analysis, _validate_browser_target_file),
        ("heatmap_target_file", full_analysis, lambda p, e: _validate_bed_target_file(p, "heatmap_target_file", e)),
        ("motif_target_file", motifs_on, lambda p, e: _validate_bed_target_file(p, "motif_target_file", e)),
        ("rnaseq_target_file", True, lambda p, e: _validate_geneid_target_file(p, "rnaseq_target_file", e)),
        ("srna_target_file", True, None),  # ShortStack accepts gff/bed/tab; existence only
    ]

    for key, gated_on, validator in checks:
        if not gated_on:
            continue
        val = config.get(key)
        if not val:
            continue
        if not check_paths:
            continue
        if not os.path.exists(val):
            # Skip default placeholder paths that the user clearly hasn't
            # touched; warn loudly otherwise. The default values shipped
            # in epicc-options.yaml all live under "data/" and end with
            # the well-known stub names, so we only flag when the path
            # has been customized away from those defaults.
            default_stubs = {
                "browser_target_file": "data/target_loci.bed",
                "heatmap_target_file": "data/target_genes.bed",
                "motif_target_file": "data/target_genes.bed",
                "rnaseq_target_file": "data/target_genes.txt",
                "srna_target_file": "config/ath.gff3",
            }
            if val == default_stubs.get(key):
                continue
            errors.append(f"[X] {key}: '{val}' does not exist")
            continue
        if validator is not None:
            validator(val, errors)

    # Background file for GO enrichment is optional; "default" is a sentinel
    if go_on:
        bg = config.get("rnaseq_background_file", "")
        if bg and bg != "default":
            if check_paths and not os.path.exists(bg):
                errors.append(
                    f"[X] rnaseq_background_file: '{bg}' does not exist"
                )
            elif check_paths:
                _validate_geneid_target_file(bg, "rnaseq_background_file", errors)

    if errors:
        full_message = "\n".join(errors)
        raise ValueError(
            f"[X] Extra output target-file validation failed — please fix "
            f"the errors below in your options file and rerun.\n"
            f"{full_message}\n\n"
        )
