#!/usr/bin/env bash
# Reproduce every result in results/ from scratch.
#
#   bash scripts/run_all.sh                  # any machine (CI, cloud container, WSL2)
#   RECORD_ENV=1 bash scripts/run_all.sh     # WSL2 `hic` env: also rewrite the
#                                            # environment record (see below)
#
# Exit codes are NOT swallowed: an environment check that FAILs (something
# installed is broken) stops the run.  SKIP (not installed yet) does not.
# Everything printed is also kept in results/logs/run_all_<UTC time>.log
# (ignored by git; attach it when reporting a problem).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PY="${PYTHON:-python}"

mkdir -p results/logs
LOG="results/logs/run_all_$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee -a "$LOG") 2>&1
echo "### log -> $LOG   (git $(git describe --always --dirty 2>/dev/null || echo unknown))"

echo "### Module 0 : environment validation"
if [ "${RECORD_ENV:-0}" = "1" ]; then
    $PY environment/verify_environment.py --write environment/versions_installed.txt
else
    # environment/versions_installed.txt is the record of the WSL2 `hic`
    # environment.  Running elsewhere (CI, a cloud container) must not
    # overwrite it -- docs/next_steps.md C3.
    $PY environment/verify_environment.py
fi

if [ "${RUN_TESTS:-1}" = "1" ]; then
    echo "### tests : unit + small-sample pipeline"
    if $PY -c "import pytest" 2>/dev/null; then
        $PY -m pytest -q
    else
        echo "  pytest not installed -- skipping (it is in environment/environment.yml)"
    fi
fi

echo "### Module 1 : generate the toy event sample (ROOT file)"
$PY scripts/m01_make_toy_events.py --config config/analysis.yaml

echo "### Module 1 : analysis (uproot arm)"
$PY analysis/m01_analysis_uproot.py --config config/analysis.yaml

echo "### Module 1 : closure test"
$PY analysis/m01_closure_test.py --config config/analysis.yaml

echo "### Module 1 : fit-method validation (pseudo-experiments, ~1 min)"
$PY analysis/m01_fit_validation.py --config config/analysis.yaml

echo "### Module 1 : Tsallis fit systematics (fit-range / binning variations)"
$PY analysis/m01_fit_systematics.py --config config/analysis.yaml

echo "### Module 1 : ROOT arms (skipped automatically if ROOT is absent)"
if $PY -c "import ROOT" 2>/dev/null; then
    $PY analysis/m01_analysis_pyroot.py
    # The C++ macro cannot read YAML: pass it the cut values from the config,
    # so the only numbers duplicated in the macro are its bin definitions
    # (and m01_compare_arms.py compares bin edges and cut values).
    MACRO_ARGS="$($PY - <<'EOF'
import yaml
c = yaml.safe_load(open("config/analysis.yaml"))
t, m = c["analysis"]["track_cuts"], c["analysis"]["muon_cuts"]
lo, hi = c["analysis"]["fits"]["tsallis_pt_range"]
print(f'"{c["paths"]["generated_file"]}", "results/tables/m01_root_histograms.root", '
      f'"results/tables/m01_root_summary.json", {t["eta_abs_max"]}, {t["pt_min"]}, '
      f'{t["pt_max"]}, {m["pt_min"]}, {m["eta_abs_max"]}, {lo}, {hi}')
EOF
)"
    root -l -b -q "analysis/m01_analysis_root.C(${MACRO_ARGS})"
    echo "### Module 1 : cross-check the three arms"
    $PY analysis/m01_compare_arms.py
else
    echo "  ROOT not installed -- skipping (see environment/setup_wsl2.sh)"
fi

echo "### done.  figures in results/figures, tables in results/tables, log in $LOG"
