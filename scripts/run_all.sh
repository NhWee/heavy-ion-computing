#!/usr/bin/env bash
# Reproduce every result in results/ from scratch.
#   bash scripts/run_all.sh
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PY="${PYTHON:-python}"

echo "### Module 0 : environment validation"
$PY environment/verify_environment.py --write environment/versions_installed.txt || true

echo "### Module 1 : generate the toy event sample (ROOT file)"
$PY scripts/m01_make_toy_events.py --config config/analysis.yaml

echo "### Module 1 : analysis (uproot arm)"
$PY analysis/m01_analysis_uproot.py --config config/analysis.yaml

echo "### Module 1 : closure test"
$PY analysis/m01_closure_test.py

echo "### Module 1 : ROOT arm (skipped automatically if ROOT is absent)"
if $PY -c "import ROOT" 2>/dev/null; then
    $PY analysis/m01_analysis_pyroot.py
    root -l -b -q 'analysis/m01_analysis_root.C("data/generated/m01_toy_events.root")'
else
    echo "  ROOT not installed -- skipping (see environment/setup_wsl2.sh)"
fi

echo "### done.  figures in results/figures, tables in results/tables"
