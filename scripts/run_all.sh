#!/usr/bin/env bash
# Reproduce every result in results/ from scratch.
#   bash scripts/run_all.sh
#
# Exit codes are NOT swallowed: an environment check that FAILs (something
# installed is broken) stops the run.  SKIP (not installed yet) does not.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PY="${PYTHON:-python}"

echo "### Module 0 : environment validation"
$PY environment/verify_environment.py --write environment/versions_installed.txt

echo "### tests : unit + small-sample pipeline"
if $PY -c "import pytest" 2>/dev/null; then
    $PY -m pytest -q
else
    echo "  pytest not installed -- skipping (it is in environment/environment.yml)"
fi

echo "### Module 1 : generate the toy event sample (ROOT file)"
$PY scripts/m01_make_toy_events.py --config config/analysis.yaml

echo "### Module 1 : analysis (uproot arm)"
$PY analysis/m01_analysis_uproot.py --config config/analysis.yaml

echo "### Module 1 : closure test"
$PY analysis/m01_closure_test.py

echo "### Module 1 : ROOT arms (skipped automatically if ROOT is absent)"
if $PY -c "import ROOT" 2>/dev/null; then
    $PY analysis/m01_analysis_pyroot.py
    root -l -b -q 'analysis/m01_analysis_root.C("data/generated/m01_toy_events.root")'
    echo "### Module 1 : cross-check the three arms"
    $PY analysis/m01_compare_arms.py
else
    echo "  ROOT not installed -- skipping (see environment/setup_wsl2.sh)"
fi

echo "### done.  figures in results/figures, tables in results/tables"
