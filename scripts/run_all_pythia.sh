#!/usr/bin/env bash
# Reproduce every PYTHIA-based result (Modules 2-4) in the WSL2 `hic` environment.
#   bash scripts/run_all_pythia.sh            # ~30 min on 4 cores
#
# The cloud container produced the committed Module 2 sensitivity / scan /
# pT-comparison and the Module 3-4 results with PyPI pythia8mc 8.318; this
# script regenerates them with the conda PYTHIA (8.312) so that every number
# in results/ comes from one environment.  Compare before committing:
#   git diff --stat results/
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PY="${PYTHON:-python}"
mkdir -p results/logs

echo "### Module 2: nominal minimum-bias sample (20k) and comparisons"
$PY scripts/m02_run_pythia.py > results/logs/m02_nominal.log 2>&1; tail -2 results/logs/m02_nominal.log
$PY analysis/m02_compare_data.py
$PY analysis/m02_compare_pt.py

echo "### Module 2: one-setting variations + sensitivity table"
bash scripts/m02_run_variations.sh
G=data/generated
$PY analysis/m02_sensitivity.py nominal=$G/m02_var_nominal.root cr_off=$G/m02_var_cr_off.root \
    tune4C=$G/m02_var_tune4C.root pt0ref240=$G/m02_var_pt0ref240.root decay_all=$G/m02_var_decay_all.root

echo "### Module 2: pT0Ref x CR-range scan (20 points x 10k)"
rm -f results/tables/m02/scan/pt0_*.json       # regenerate, do not reuse cloud points
$PY analysis/m02_tune_scan.py generate --n-events 10000
$PY analysis/m02_tune_scan.py fit

echo "### Module 3: HepMC event record"
$PY scripts/m03_write_hepmc.py --n-events 100
$PY analysis/m03_event_record.py

echo "### Module 4: jets"
$PY scripts/m02_run_pythia.py --card generators/pythia_hardqcd_13TeV.cmnd \
    --out $G/m04_hardqcd_13TeV.root > results/logs/m04_gen.log 2>&1; tail -2 results/logs/m04_gen.log
$PY analysis/m04_jets.py
echo "### done -- review: git diff --stat results/"
