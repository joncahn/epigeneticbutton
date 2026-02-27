Run S. pombe pipeline validation. Uses `scripts/validate_pombe.sh` to orchestrate three stages:

- `--dry` — Snakemake dry-run integration tests (fast, ~42s)
- `--full` — Execute the full pipeline (auto-detects SLURM vs local)
- `--check` — Validate completed pipeline outputs against expected files

## Workflow

1. **After code changes**: always run `--dry` first to catch config/rule errors
2. **If runtime behavior changed** (shell commands, new rules, resource changes): run `--full`
3. **After a full run completes**: run `--check` to validate all outputs

## When to skip

- Documentation-only changes (README, CLAUDE.md)
- Changes to non-pombe configs or sample sheets
- Python-only refactors already covered by unit tests

## Usage

```bash
bash scripts/validate_pombe.sh $ARGUMENTS
```

If no arguments are provided, pass `--dry` as the default. Common patterns:

- Quick validation: `--dry`
- Full cycle: `--all` (runs dry → full → check, stops on failure)
- Post-run only: `--check`
