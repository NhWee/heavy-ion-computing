#!/usr/bin/env bash
# Reproduce every PYTHIA-based result (Modules 2-4) in the WSL2 `hic` environment.
#
#   bash scripts/run_all_pythia.sh                 # all stages (~30 min on 4 cores)
#   bash scripts/run_all_pythia.sh compare scan m03 m04   # only these stages
#
# Stages, in order:
#   nominal     generate the 20k minimum-bias sample
#   compare     dN/deta and pT spectrum vs ALICE (needs the nominal sample)
#   variations  generate the one-setting variations (5 x 20k)
#   sensitivity sensitivity table from the variations
#   scan        pT0Ref x CR-range scan (20 x 10k) + fit
#   stringpt    StringPT:sigma scan at pT0Ref 2.40 (6 x 20k) + fit
#   m03         HepMC event record (Module 3)
#   m04         hard-QCD sample + jets (Module 4)
#   m05         proton and nuclear PDFs (Module 5; no generation)
#
# Before anything runs, the ALICE tables (HEPData payloads, NOT in git) are
# checked: present, and with the sha256 recorded in their provenance JSON.
# A missing table stops the script here, before any result is touched
# (it used to fail half-way, after deleting the old scan results).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PY="${PYTHON:-python}"
G=data/generated
mkdir -p results/logs

ALL="nominal compare variations sensitivity scan stringpt m03 m04 m05"
STAGES="${*:-$ALL}"
want() { [[ " $STAGES " == *" $1 "* ]]; }

# ---------------------------------------------------------------- preflight
need_data=0
for s in compare sensitivity scan stringpt; do want "$s" && need_data=1; done
if [[ $need_data == 1 ]]; then
    "$PY" - <<'EOF'
import hashlib, json, pathlib, sys
bad = []
for j in ("data/raw/alice_pp13TeV_dndeta.json", "data/raw/alice_pp13TeV_pt.json"):
    prov = json.loads(pathlib.Path(j).read_text())
    f = pathlib.Path(j).parent / prov["file"]
    if not f.exists():
        bad.append(f"MISSING  {f}   (HEPData {prov['hepdata']['record']} {prov['hepdata']['table']})")
    elif prov.get("sha256") and hashlib.sha256(f.read_bytes()).hexdigest() != prov["sha256"]:
        bad.append(f"CHANGED  {f}   (sha256 differs from {j})")
if bad:
    print("preflight: ALICE tables not ready -- nothing has been run or changed.")
    print("\n".join("  " + b for b in bad))
    print("Copy them from where the browser saved them, e.g.\n"
          "  cp /mnt/c/Users/<you>/Downloads/HEPData-ins1395253-v1-Table_1.yaml data/raw/alice_pp13TeV_dndeta_table1.yaml\n"
          "  cp /mnt/c/Users/<you>/Downloads/HEPData-ins1395253-v1-Table_2.yaml data/raw/alice_pp13TeV_pt_table2.yaml")
    sys.exit(1)
print("preflight: ALICE tables present, sha256 match")
EOF
fi
for s in compare; do
    if want "$s" && ! want nominal && [[ ! -f $G/m02_pythia_mb_13TeV.root ]]; then
        echo "stage '$s' needs $G/m02_pythia_mb_13TeV.root -- add the 'nominal' stage"; exit 1
    fi
done
if want sensitivity && ! want variations && [[ ! -f $G/m02_var_nominal.root ]]; then
    echo "stage 'sensitivity' needs the variation samples -- add the 'variations' stage"; exit 1
fi

# ---------------------------------------------------------------- stages
if want nominal; then
    echo "### Module 2: nominal minimum-bias sample (20k)"
    $PY scripts/m02_run_pythia.py > results/logs/m02_nominal.log 2>&1; tail -2 results/logs/m02_nominal.log
fi
if want compare; then
    echo "### Module 2: dN/deta and pT spectrum vs ALICE"
    $PY analysis/m02_compare_data.py
    $PY analysis/m02_compare_pt.py
fi
if want variations; then
    echo "### Module 2: one-setting variations"
    bash scripts/m02_run_variations.sh
fi
if want sensitivity; then
    echo "### Module 2: sensitivity table"
    $PY analysis/m02_sensitivity.py nominal=$G/m02_var_nominal.root cr_off=$G/m02_var_cr_off.root \
        tune4C=$G/m02_var_tune4C.root pt0ref240=$G/m02_var_pt0ref240.root decay_all=$G/m02_var_decay_all.root
fi
if want scan; then
    echo "### Module 2: pT0Ref x CR-range scan (20 points x 10k)"
    rm -f results/tables/m02/scan/pt0_*.json   # regenerate; the preflight guarantees the data exist
    $PY analysis/m02_tune_scan.py generate --n-events 10000
    $PY analysis/m02_tune_scan.py fit
fi
if want stringpt; then
    echo "### Module 2: StringPT:sigma scan (6 points x 20k)"
    rm -f results/tables/m02/stringpt/stringpt_*.json
    $PY analysis/m02_stringpt_scan.py generate --n-events 20000
    $PY analysis/m02_stringpt_scan.py fit
fi
if want m03; then
    echo "### Module 3: HepMC event record"
    $PY scripts/m03_write_hepmc.py --n-events 100
    $PY analysis/m03_event_record.py
fi
if want m04; then
    echo "### Module 4: jets"
    $PY scripts/m02_run_pythia.py --card generators/pythia_hardqcd_13TeV.cmnd \
        --out $G/m04_hardqcd_13TeV.root > results/logs/m04_gen.log 2>&1; tail -2 results/logs/m04_gen.log
    $PY analysis/m04_jets.py
fi
if want m05; then
    echo "### Module 5: proton and nuclear PDFs"
    $PY analysis/m05_pdfs.py
fi
echo "### done ($STAGES) -- review: git diff --stat results/"
