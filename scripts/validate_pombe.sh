#!/usr/bin/env bash
#
# validate_pombe.sh — Orchestrate S. pombe pipeline validation stages.
#
# Usage:
#   bash scripts/validate_pombe.sh --dry     # Run dry-run integration tests
#   bash scripts/validate_pombe.sh --full    # Execute the full Snakemake pipeline
#   bash scripts/validate_pombe.sh --check   # Validate completed pipeline outputs
#   bash scripts/validate_pombe.sh --all     # All three in sequence
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

CONFIG_FILE="tests/integration/data/test_config_pombe.yaml"
SAMPLE_FILE="tests/integration/data/test_samples_pombe.tsv"

# Colors for output (if terminal supports it)
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    YELLOW='\033[0;33m'
    NC='\033[0m'
else
    GREEN='' RED='' YELLOW='' NC=''
fi

usage() {
    echo "Usage: $0 [--dry] [--full] [--check] [--all]"
    echo ""
    echo "  --dry    Run dry-run integration tests (pytest)"
    echo "  --full   Execute the full Snakemake pipeline"
    echo "  --check  Validate completed pipeline outputs (pytest)"
    echo "  --all    Run all three stages in sequence"
    echo ""
    echo "Stages abort on failure. --all stops before --full if --dry fails."
    exit 1
}

run_dry() {
    echo -e "${YELLOW}=== Stage: Dry-Run Tests ===${NC}"
    conda run -n smk9 pytest tests/integration/test_pombe_dryrun.py -v
    echo -e "${GREEN}=== Dry-run tests PASSED ===${NC}"
}

run_full() {
    echo -e "${YELLOW}=== Stage: Full Pipeline Run ===${NC}"

    if ! [ -f "$CONFIG_FILE" ]; then
        echo -e "${RED}ERROR: Config file not found: $CONFIG_FILE${NC}" >&2
        exit 1
    fi
    if ! [ -f "$SAMPLE_FILE" ]; then
        echo -e "${RED}ERROR: Sample file not found: $SAMPLE_FILE${NC}" >&2
        exit 1
    fi

    if command -v sbatch &>/dev/null; then
        echo -e "  SLURM detected — using profiles/slurm"
        conda run -n smk9 snakemake \
            --profile profiles/slurm \
            --configfile "$CONFIG_FILE"
    else
        CORES=$(( $(nproc) / 2 ))
        [ "$CORES" -lt 1 ] && CORES=1
        echo -e "  No SLURM — running locally with $CORES cores"
        conda run -n smk9 snakemake \
            --use-conda --conda-frontend conda \
            --cores "$CORES" \
            --configfile "$CONFIG_FILE"
    fi

    echo -e "${GREEN}=== Full pipeline run COMPLETED ===${NC}"
}

run_check() {
    echo -e "${YELLOW}=== Stage: Post-Run Validation ===${NC}"
    conda run -n smk9 pytest tests/integration/test_pombe_postrun.py -v -m slow
    echo -e "${GREEN}=== Post-run validation PASSED ===${NC}"
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
if [ $# -eq 0 ]; then
    usage
fi

DO_DRY=false
DO_FULL=false
DO_CHECK=false

while [ $# -gt 0 ]; do
    case "$1" in
        --dry)   DO_DRY=true ;;
        --full)  DO_FULL=true ;;
        --check) DO_CHECK=true ;;
        --all)   DO_DRY=true; DO_FULL=true; DO_CHECK=true ;;
        -h|--help) usage ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}" >&2
            usage
            ;;
    esac
    shift
done

# ---------------------------------------------------------------------------
# Execute stages in order
# ---------------------------------------------------------------------------
if $DO_DRY; then
    run_dry
fi

if $DO_FULL; then
    run_full
fi

if $DO_CHECK; then
    run_check
fi

echo -e "${GREEN}=== All requested stages completed ===${NC}"
