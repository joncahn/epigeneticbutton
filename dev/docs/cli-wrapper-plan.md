# EPICC CLI Wrapper — Implementation Plan

## Overview

Replace direct `snakemake` invocations with a user-friendly `epicc` wrapper script. The wrapper exposes key pipeline options as CLI flags, auto-detects execution mode (local vs SLURM), supports configurable output directories, and passes through arbitrary snakemake arguments.

## Language Choice: Python

- Consistent with the rest of the codebase (Snakefile, sample_sheet.py, profiler)
- `argparse` provides structured subcommands with auto-generated help
- Native YAML reading (PyYAML already in the conda env)
- No new dependencies — runs in the existing `epicc` conda environment

## CLI Interface

```
epicc run [options] [-- SNAKEMAKE_ARGS...]
epicc dry-run [options] [-- SNAKEMAKE_ARGS...]
epicc validate [options]
epicc profile [options]
epicc clean [options]
epicc init [DIR]
```

### Global options (before subcommand)

| Flag | Default | Description |
|------|---------|-------------|
| `--options FILE` | `config/epicc-options.yaml` | Path to options file |
| `--samples FILE` | Value from options file | Path to sample sheet TSV |
| `--output-dir DIR` | `results` | Output directory prefix |
| `--genome-dir DIR` | `genomes` | Genome directory prefix |
| `--cores N` | Half of `nproc` | Cores for local execution |
| `--profile PROFILE` | Auto-detect | Snakemake profile directory |
| `--keep-intermediates TIER` | From options file | none/standard/custom/all |
| `--verbose` / `--quiet` | Normal | Control output level |

### Subcommands

**`epicc run`** — Execute the pipeline.
- `--target TARGET` — Named target (default: `all`)
- `--until RULE` — Stop after this rule
- `--forcerun RULE [RULE ...]` — Force re-execution of specific rules
- Everything after `--` passed verbatim to snakemake

**`epicc dry-run`** — Same as `run` but adds `--dry-run`.

**`epicc validate`** — Load options + sample sheet, check config, print summary. Optionally lint with `snakemake --dry-run --lint`.

**`epicc profile`** — Wrapper around `dev/profile_snakemake_log.py`.
- `--latest | LOGFILE` — Which log to profile
- `--html OUTPUT.html` — Generate HTML report

**`epicc clean`** — Selective cleanup.
- `--intermediates` — Remove temp files based on retention tier
- `--conda-envs` — Run `snakemake --sdm conda --conda-cleanup-envs`
- `--logs` — Clean old `.snakemake/log/` files
- `--all` — All of the above

**`epicc init [DIR]`** — Scaffold a new analysis directory with template options file and sample sheet.

## Config Override Mechanism

CLI flags map to snakemake `--config key=value` overrides:

| CLI flag | Config key |
|----------|-----------|
| `--output-dir DIR` | `output_dir=DIR` |
| `--genome-dir DIR` | `genome_dir=DIR` |
| `--samples FILE` | `sample_file=FILE` |
| `--keep-intermediates TIER` | `keep_intermediates=TIER` |

The `--options` flag maps to snakemake's `--configfile`.

## Profile Auto-Detection

1. If `--profile` explicitly given → use it
2. If `--cores` explicitly given (no `--profile`) → local mode
3. Otherwise: `shutil.which("sbatch")` → `profiles/slurm` if found, else local with half of `nproc`

## Output Directory (Prerequisite)

This is the largest prerequisite and should be tackled first.

### Current state

`results/` appears ~616 times and `genomes/` ~100 times across 8 rule files + Snakefile, all as hardcoded string prefixes.

### Approach

Define config-derived variables in the Snakefile:

```python
RESULTS_DIR = config.get("output_dir", "results")
GENOMES_DIR = config.get("genome_dir", "genomes")
```

Mechanically replace all occurrences in rule files. Two patterns:

1. **f-string paths** (input functions, expand calls): `f"results/{env}/..."` → `f"{RESULTS_DIR}/{env}/..."`
2. **Snakemake rule strings** (input/output directives): `"results/{env}/mapped/..."` → needs f-string wrapping with doubled wildcards: `f"{RESULTS_DIR}/{{env}}/mapped/..."`

### Alternatives considered

- **Symlinks** (`results/` → target): Fragile, breaks on move, doesn't support concurrent runs. Rejected.
- **`workdir:` directive**: Changes CWD for all rules, breaks relative paths in config. Rejected.

## Snakemake Command Construction

```python
cmd = ["snakemake"]
cmd += ["--configfile", options_file]

if profile:
    cmd += ["--profile", profile]
else:
    cmd += ["--use-conda", "--conda-frontend", "conda", "--cores", str(cores)]

config_overrides = []
if output_dir != "results":
    config_overrides.append(f"output_dir={output_dir}")
# ... other overrides
if config_overrides:
    cmd += ["--config"] + config_overrides

if target != "all":
    cmd.append(target)

cmd += extra_snakemake_args
```

Uses `subprocess.run()` (not `os.exec`) so the wrapper can capture exit codes and print a completion summary.

## Implementation Phases

### Phase 1: Configurable output directory

1. Add `output_dir` / `genome_dir` config keys to Snakefile with defaults
2. Replace all `results/` occurrences across rule files (~616)
3. Replace all `genomes/` occurrences across rule files (~100)
4. Update test configs and dry-run tests to exercise non-default output dir
5. Verify with pombe dry-run and full run

### Phase 2: Core CLI wrapper

1. Create `epicc` script with argparse skeleton and subcommands
2. Implement `run` and `dry-run` — command construction, profile auto-detection, passthrough
3. Implement `validate` — config + sample sheet loading, summary output
4. Wire all global flags to config overrides
5. Test locally and on cluster

### Phase 3: Ancillary subcommands

1. `epicc profile` wrapping `dev/profile_snakemake_log.py`
2. `epicc clean` for intermediate/log/conda cleanup
3. `epicc init` for scaffolding new analysis directories

### Phase 4: Polish

1. Update README.md, CLAUDE.md with new invocation patterns
2. Optionally refactor `validate_pombe.sh` to use `epicc`
3. Add `dev/docs/design-decisions.md` entry
4. Optional: `pyproject.toml` for pip-installable entry point

## Files to Create

- `epicc` (repo root) — Main CLI script (~300-400 lines), executable with `#!/usr/bin/env python3`

## Files to Modify

### Phase 1 (output_dir) — all rule files + Snakefile

| File | `results/` refs | `genomes/` refs |
|------|----------------|-----------------|
| `workflow/rules/combined_analysis.smk` | ~177 | ~9 |
| `workflow/rules/ChIPseq.smk` | ~121 | ~16 |
| `workflow/rules/RNAseq.smk` | ~101 | ~14 |
| `workflow/rules/mC.smk` | ~98 | ~15 |
| `workflow/rules/smallRNA.smk` | ~67 | ~14 |
| `workflow/rules/ATACseq.smk` | ~25 | ~1 |
| `workflow/rules/sample_download.smk` | ~23 | — |
| `workflow/rules/environment_setup.smk` | ~4 | ~30 |
| `workflow/Snakefile` | targets | config setup |
| `config/epicc-options.yaml` | new keys | new keys |

### Phase 2 (wrapper)

- `README.md` — Update running instructions
- `CLAUDE.md` — Update running instructions
- `dev/docs/TODO.md` — Mark items done

## Design Decisions

**`--` passthrough, not `--smk`**: Standard POSIX convention for end-of-options. More intuitive than a custom flag. Supported by argparse via `REMAINDER`.

**Subcommands over flat flags**: Clearer interface, mirrors git/docker/conda conventions. Each subcommand has focused help text.

**Auto-detect SLURM**: Avoids requiring users to remember different invocations. Explicit `--profile` and `--cores` serve as overrides.

**Full path replacement over symlinks**: More upfront work, but robust — supports concurrent runs, no filesystem edge cases, transparent to users.

**`subprocess.run()` over `os.exec`**: Allows the wrapper to print a post-run summary (elapsed time, log path, suggestion to run `epicc profile`).
