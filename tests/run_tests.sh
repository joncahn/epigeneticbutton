#!/bin/bash
# Run unit tests for EpigeneticButton dmC (direct methylation) feature
#
# Usage:
#   ./tests/run_tests.sh              # Run all tests
#   ./tests/run_tests.sh -v           # Verbose output
#   ./tests/run_tests.sh --cov        # With coverage report
#   ./tests/run_tests.sh --help       # Show help

set -euo pipefail

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "$REPO_ROOT"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}EpigeneticButton Test Suite${NC}"
echo "============================="
echo ""

# Check if pytest is installed
if ! command -v pytest &> /dev/null; then
    echo -e "${RED}Error: pytest is not installed${NC}"
    echo "Please install test dependencies:"
    echo "  pip install -r tests/requirements-test.txt"
    exit 1
fi

# Parse arguments
RUN_COVERAGE=false
PYTEST_ARGS=()

for arg in "$@"; do
    case $arg in
        --cov|--coverage)
            RUN_COVERAGE=true
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --cov, --coverage    Generate coverage report"
            echo "  -v, --verbose        Verbose output"
            echo "  -k EXPRESSION        Run tests matching expression"
            echo "  --help, -h           Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                              # Run all tests"
            echo "  $0 -v                           # Verbose output"
            echo "  $0 --cov                        # With coverage"
            echo "  $0 -k test_validate_bedmethyl   # Run specific tests"
            exit 0
            ;;
        *)
            PYTEST_ARGS+=("$arg")
            ;;
    esac
done

# Build pytest command
PYTEST_CMD="pytest tests/unit/"

if [ "$RUN_COVERAGE" = true ]; then
    PYTEST_CMD="$PYTEST_CMD --cov=workflow/scripts --cov-report=html --cov-report=term"
    echo -e "${YELLOW}Running tests with coverage...${NC}"
else
    echo -e "${YELLOW}Running tests...${NC}"
fi

# Add any additional arguments
if [ ${#PYTEST_ARGS[@]} -gt 0 ]; then
    PYTEST_CMD="$PYTEST_CMD ${PYTEST_ARGS[*]}"
fi

echo "Command: $PYTEST_CMD"
echo ""

# Run tests
if $PYTEST_CMD; then
    echo ""
    echo -e "${GREEN}✓ All tests passed!${NC}"

    if [ "$RUN_COVERAGE" = true ]; then
        echo ""
        echo -e "${YELLOW}Coverage report generated in: htmlcov/index.html${NC}"
    fi

    exit 0
else
    echo ""
    echo -e "${RED}✗ Some tests failed${NC}"
    exit 1
fi
