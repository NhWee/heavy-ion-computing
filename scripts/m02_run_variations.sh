#!/usr/bin/env bash
# Generate the PYTHIA samples for the dN/deta sensitivity study (Module 2).
#   bash scripts/m02_run_variations.sh [N_EVENTS]
# Then:
#   python analysis/m02_sensitivity.py nominal=data/generated/m02_var_nominal.root \
#       cr_off=data/generated/m02_var_cr_off.root tune4C=data/generated/m02_var_tune4C.root \
#       pt0ref240=data/generated/m02_var_pt0ref240.root decay_all=data/generated/m02_var_decay_all.root
# Each variation changes ONE setting on top of generators/pythia_mb_13TeV.cmnd;
# the settings are recorded in each file's provenance (extra_settings).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PY="${PYTHON:-python}"
N="${1:-20000}"
OUT=data/generated
run() {  # name, then zero or more --set arguments
    local name="$1"; shift
    echo "### variation: $name $*"
    "$PY" scripts/m02_run_pythia.py --n-events "$N" --out "$OUT/m02_var_${name}.root" "$@" \
        > "results/logs/m02_var_${name}.log" 2>&1
    tail -3 "results/logs/m02_var_${name}.log"
}
run nominal
run cr_off    --set "ColourReconnection:reconnect = off"
run tune4C    --set "Tune:pp = 5"
run pt0ref240 --set "MultipartonInteractions:pT0Ref = 2.40"
run decay_all --set "ParticleDecays:limitTau0 = off"
