#!/usr/bin/env bash
# Reproduce every PYTHIA-based result (Modules 2-6) in the WSL2 `hic` environment.
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
#   m05         proton and nuclear PDFs (Module 5; no generation; uses .venv-pdf, see the stage)
#   m06         Monte Carlo Glauber, Pb-Pb and p-Pb (Module 6; numpy only, ~10 min on 4 cores)
#   origin      where PYTHIA pi/K/p come from + 8 settings (direct vs resonance decays,
#               vector-meson fraction, Bose-Einstein, rescattering); after "identified";
#               NOT in the default list; ~3 min on 8+ cores
#   identified  pi/K/p spectra vs ALICE (NOT in the default list: needs the HEPData
#               ins1784041 tables, docs/howto/fetch_hepdata.md appendix; needs 'nominal'
#               sample; generates one 20k tuned sample)
#
# Before anything runs, the ALICE tables (HEPData payloads, NOT in git) are
# checked: present, and with the sha256 recorded in their provenance JSON.
# A missing table stops the script here, before any result is touched
# (it used to fail half-way, after deleting the old scan results).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PY="${PYTHON:-python}"
# Fail before touching anything if this is not the hic environment (a fresh
# terminal starts outside it: no `python`, or a system python without numpy).
if ! command -v "$PY" >/dev/null 2>&1 || ! "$PY" -c "import numpy, uproot, awkward" >/dev/null 2>&1; then
    echo "'$PY' with numpy/uproot/awkward not found -- activate the environment first:"
    echo "    micromamba activate hic"
    exit 1
fi
G=data/generated
mkdir -p results/logs
# Steps whose output goes to a log: on failure print the log's tail and stop.
# (Plain `cmd > log 2>&1` under `set -e` used to end the script with no message.)
logged() {
    local log=$1 rc=0; shift
    "$@" > "$log" 2>&1 || rc=$?
    if (( rc != 0 )); then
        echo "FAILED (exit $rc): $*"; echo "--- last lines of $log:"; tail -n 25 "$log"; exit 1
    fi
}

ALL="nominal compare variations sensitivity scan stringpt m03 m04 m05 m06"
STAGES="${*:-$ALL}"
want() { [[ " $STAGES " == *" $1 "* ]]; }

# ---------------------------------------------------------------- preflight
# PYTHIA's Python code and its XML data must be the same version.  They stop
# matching if the PyPI wheel pythia8mc is pip-installed INTO the conda env: it
# writes its own share/Pythia8/xmldoc over conda's (troubleshooting #24).
need_gen=0
for s in nominal variations scan stringpt m03 m04 identified origin; do want "$s" && need_gen=1; done
if [[ $need_gen == 1 ]]; then
    out=$("$PY" -c "
try:
    import pythia8
except ImportError:
    import pythia8mc as pythia8
p = pythia8.Pythia('', False)
print('PYTHIA', pythia8.__name__, p.settings.parm('Pythia:versionNumber'))" 2>&1 || true)
    if [[ $out == *Abort* || $out != *PYTHIA* ]]; then
        echo "preflight: PYTHIA cannot start -- nothing has been run or changed."
        echo "$out" | grep -v "^ *$" | tail -n 3
        echo "If it says 'unmatched version numbers', repair the env (troubleshooting #24):"
        echo "    pip uninstall -y pythia8mc"
        echo "    micromamba install -y -n hic -c conda-forge --force-reinstall pythia8=8.312"
        exit 1
    fi
    echo "preflight: ${out##*$'\n'}"
fi
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
if want identified || want origin; then
    # Not unpacked yet: look for the archive in the Windows Downloads folder(s)
    # (WSL2 mounts C: at /mnt/c) and unpack it -- no user name to type.
    if [[ ! -f data/raw/ins1784041/Table_1.yaml ]]; then
        arc=$(ls -t /mnt/c/Users/*/Downloads/HEPData-ins1784041*yaml*.tar.gz 2>/dev/null | head -n 1 || true)
        if [[ -n $arc ]]; then
            echo "preflight: unpacking $arc -> data/raw/ins1784041/"
            mkdir -p data/raw/ins1784041
            tar xzf "$arc" --strip-components=1 -C data/raw/ins1784041
        fi
    fi
    for n in 1 2 3 4 5 6 7 8 9 10 11 12; do
        if [[ ! -f data/raw/ins1784041/Table_$n.yaml ]]; then
            echo "preflight: data/raw/ins1784041/Table_$n.yaml missing -- unpack the HEPData ins1784041"
            echo "           'Download All -> YAML' archive there (docs/howto/fetch_hepdata.md, appendix)"; exit 1
        fi
    done
fi
for s in compare identified; do
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
    logged results/logs/m02_nominal.log $PY scripts/m02_run_pythia.py; tail -n 2 results/logs/m02_nominal.log
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
    logged results/logs/m04_gen.log $PY scripts/m02_run_pythia.py --card generators/pythia_hardqcd_13TeV.cmnd \
        --out $G/m04_hardqcd_13TeV.root; tail -n 2 results/logs/m04_gen.log
    $PY analysis/m04_jets.py
fi
if want m05; then
    echo "### Module 5: proton and nuclear PDFs"
    # The PDF classes come from pythia8mc, which lives in its own venv so that
    # it cannot overwrite the conda PYTHIA data (#22, #24).  Create it once:
    #   python -m venv --system-site-packages .venv-pdf && .venv-pdf/bin/pip install pythia8mc
    PDFPY=$PY
    [[ -x .venv-pdf/bin/python ]] && PDFPY=.venv-pdf/bin/python
    if ! "$PDFPY" -c "import sys; sys.path.insert(0, 'src'); from hicpy.pdfs import pythia_module; pythia_module()" 2>/dev/null; then
        echo "no PYTHIA module with PDF classes for '$PDFPY' -- create the PDF venv once:"
        echo "    python -m venv --system-site-packages .venv-pdf && .venv-pdf/bin/pip install pythia8mc"
        exit 1
    fi
    $PDFPY analysis/m05_pdfs.py
fi
if want m06; then
    echo "### Module 6: Monte Carlo Glauber"
    logged results/logs/m06_glauber.log $PY analysis/m06_glauber.py; grep -E "^\[" results/logs/m06_glauber.log
fi
if want identified; then
    echo "### Module 2: identified pi/K/p vs ALICE"
    $PY scripts/m02_derive_identified_inel0.py
    $PY analysis/m02_compare_identified.py --label Monash
    logged results/logs/id_gen_best.log $PY scripts/m02_run_pythia.py --set "MultipartonInteractions:pT0Ref = 2.40" --set "StringPT:sigma = 0.31" \
        --out $G/m02_id_pt0240_sig031.root; tail -n 1 results/logs/id_gen_best.log
    $PY analysis/m02_compare_identified.py --mc $G/m02_id_pt0240_sig031.root \
        --label "pT0Ref 2.40, sigma 0.31" --tag _pt0240_sig031
fi
if want origin; then
    echo "### Module 2: origin of pi/K/p in PYTHIA, 8 settings in parallel (logs: results/logs/origin_*.log)"
    $PY scripts/m02_derive_identified_inel0.py > /dev/null
    run_origin() {   # tag label [--set ...]: one background job, log tail on failure
        local tag=$1 label=$2; shift 2
        logged "results/logs/origin${tag:-_monash}.log" $PY analysis/m02_pion_origin.py --tag "$tag" --label "$label" "$@"
    }
    pids=()
    run_origin ""                         "Monash" & pids+=($!)
    run_origin _pt0240_sig031             "pT0Ref 2.40, sigma 0.31" \
        --set "MultipartonInteractions:pT0Ref = 2.40" --set "StringPT:sigma = 0.31" & pids+=($!)
    run_origin _pt0240_sig031_rescatter   "pT0Ref 2.40, sigma 0.31, rescattering" \
        --set "MultipartonInteractions:pT0Ref = 2.40" --set "StringPT:sigma = 0.31" \
        --set "HadronLevel:Rescatter = on" & pids+=($!)
    run_origin _udv025 "UDvector 0.25" --set "StringFlav:mesonUDvector = 0.25" & pids+=($!)
    run_origin _udv10  "UDvector 1.0"  --set "StringFlav:mesonUDvector = 1.0"  & pids+=($!)
    run_origin _udv15  "UDvector 1.5"  --set "StringFlav:mesonUDvector = 1.5"  & pids+=($!)
    run_origin _bose      "Bose-Einstein on" --set "HadronLevel:BoseEinstein = on" & pids+=($!)
    run_origin _rescatter "rescattering on"  --set "HadronLevel:Rescatter = on"    & pids+=($!)
    for pid in "${pids[@]}"; do wait "$pid" || exit 1; done
    grep -h "^\[closure\]" results/logs/origin_*.log || true
    $PY analysis/m02_pion_origin_summary.py
fi
echo "### done ($STAGES) -- review: git diff --stat results/"
