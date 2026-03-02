"""Validation for new-format sample sheets.

Implements the rules defined in dev/docs/sample-sheet-spec.md.
"""

import re
from scripts.sample_sheet import VALID_ASSAYS

# Characters that are unsafe for filesystem use in Sample_ID
_UNSAFE_CHARS = re.compile(r'[/\\\s\'\";&|<>$`!{}()\[\]?*~#]')
_DOUBLE_UNDERSCORE = re.compile(r'__')
_SRA_REGEX = re.compile(r'^[SDE]RR\d+$')

_CHIP_ASSAYS = {"ChIP_broad", "ChIP_narrow"}
_CONTROL_ASSAYS = {"ChIP_broad", "ChIP_narrow", "RAMPAGE"}


def check_table(tab):
    """Validate a new-format sample sheet DataFrame.

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
        all_factors.append(factor_names)

    # Check consistent factor count and names
    if all_factors:
        ref_factors = all_factors[0]
        for i, factors in enumerate(all_factors[1:], start=2):
            if len(factors) != len(ref_factors):
                errors.append(
                    f"[X] Row #{i}: Levels has {len(factors)} factors, "
                    f"expected {len(ref_factors)} (matching row 1)"
                )
            elif factors != ref_factors:
                errors.append(
                    f"[X] Row #{i}: Levels factor names {factors} "
                    f"don't match row 1 {ref_factors}"
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
                    warnings.append(
                        f"[!] Row #{i} '{sid}': Read_layout is PE but Read_files "
                        f"has only one path (expected comma-separated pair)"
                    )
                if layout == "SE" and len(files_in_comp) > 1:
                    warnings.append(
                        f"[!] Row #{i} '{sid}': Read_layout is SE but Read_files "
                        f"has multiple comma-separated paths"
                    )

    # --- IP_target: required for ChIP, blank for others ---
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        assay = str(row.get("Assay", "")).strip()
        ip_target = str(row.get("IP_target", "")).strip()
        sid = str(row.get("Sample_ID", "")).strip()
        if ip_target == "nan":
            ip_target = ""
        if assay in _CHIP_ASSAYS:
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
