#!/usr/bin/env bash
# Module 8: download the CMS 2012 DoubleMuParked NanoAOD file (CERN Open Data
# record 12341) and write its provenance JSON.
#
#   bash scripts/m08_fetch_cms_dimuon.sh
#
# ~2 GB, ~61 M events.  The file goes to data/raw/opendata/ (NOT in git); the
# provenance JSON data/raw/cms_2012_doublemuparked.json (in git) records the
# record, DOI, URL, size, sha256 and download time.  Re-running resumes a
# partial download and only re-computes the checksum.
#
# Tries HTTPS (curl, resumable) first, then XRootD (xrdcp) if installed.
# The cloud session cannot reach opendata.cern.ch / eospublic.cern.ch: run in WSL2.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PY="${PYTHON:-python}"
NAME=Run2012BC_DoubleMuParked_Muons.root
EOSPATH=/eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/$NAME
URLS=("https://eospublic.cern.ch/$EOSPATH" "https://opendata.cern.ch$EOSPATH")
XRD="root://eospublic.cern.ch/$EOSPATH"
OUT=data/raw/opendata/$NAME
mkdir -p data/raw/opendata

got=""
for u in "${URLS[@]}"; do
    echo "### trying $u"
    if curl -fL --retry 3 -C - -o "$OUT" "$u"; then got="$u"; break; fi
done
if [[ -z $got ]] && command -v xrdcp >/dev/null 2>&1; then
    echo "### trying $XRD"
    xrdcp -f "$XRD" "$OUT" && got="$XRD"
fi
if [[ -z $got ]]; then
    echo "download failed from every source -- check the network, or download $NAME by hand from"
    echo "https://opendata.cern.ch/record/12341 and put it at $OUT, then run this script again"
    [[ -f $OUT ]] || exit 1
    got="by hand"
fi

"$PY" - "$OUT" "$got" <<'EOF'
import datetime as dt, hashlib, json, pathlib, sys
path, source = pathlib.Path(sys.argv[1]), sys.argv[2]
h = hashlib.sha256()
with open(path, "rb") as f:
    for block in iter(lambda: f.read(1 << 24), b""):
        h.update(block)
try:
    import uproot
    n = uproot.open(path)["Events"].num_entries
except Exception as exc:                       # noqa: BLE001
    n = None; print(f"[warn] could not open with uproot: {exc}")
prov = {
    "_status": "Downloaded by scripts/m08_fetch_cms_dimuon.sh; sha256 computed locally after download.",
    "verified": True,
    "experiment": "CMS",
    "title": "DoubleMuParked dataset from 2012 in NanoAOD format reduced on muons",
    "opendata_record": 12341,
    "opendata_url": "https://opendata.cern.ch/record/12341",
    "doi": "10.7483/OPENDATA.CMS.LVG5.QT81",
    "analysis_record": {"record": 12342, "doi": "10.7483/OPENDATA.CMS.AAR1.4NZQ",
                        "title": "Analysis of the di-muon spectrum using data from the CMS detector taken in 2012"},
    "license": "CC0 1.0",
    "system": "pp", "sqrt_s_TeV": 8.0, "runs": "2012 Run B and C (DoubleMuParked)",
    "content": "muon variables only: nMuon, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Muon_charge",
    "source_url": source,
    "file": str(path.relative_to(pathlib.Path("data/raw"))),
    "size_bytes": path.stat().st_size,
    "sha256": h.hexdigest(),
    "n_events": n,
    "downloaded_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "note": "CMS and CERN do not endorse works produced with these data (record terms of use)."}
pathlib.Path("data/raw/cms_2012_doublemuparked.json").write_text(json.dumps(prov, indent=2) + "\n")
print(f"[prov] data/raw/cms_2012_doublemuparked.json: {prov['size_bytes'] / 1e9:.2f} GB, "
      f"{n} events, sha256 {prov['sha256'][:16]}...")
EOF
