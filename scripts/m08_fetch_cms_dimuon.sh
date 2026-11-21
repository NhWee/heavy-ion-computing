#!/usr/bin/env bash
# Module 8: download the CMS 2012 DoubleMuParked NanoAOD file (CERN Open Data
# record 12341) and write its provenance JSON.
#
#   bash scripts/m08_fetch_cms_dimuon.sh
#   nohup bash scripts/m08_fetch_cms_dimuon.sh > results/logs/m08_fetch.log 2>&1 &   # in the background
#
# 2.08 GB, ~61 M events.  The file goes to data/raw/opendata/ (NOT in git); the
# provenance JSON data/raw/cms_2012_doublemuparked.json (in git) is written ONLY
# when the file is complete (size equal to the server's Content-Length) and
# opens as a ROOT file with an "Events" tree.  A partial file is kept and the
# next run resumes it.  Idempotent: a complete file is not downloaded again.
#
# Source: https://opendata.cern.ch/eos/opendata/... (HTTPS, resumable).  The
# eospublic.cern.ch HTTPS endpoint uses the CERN Grid CA, which is not in the
# usual system trust store ("self-signed certificate in certificate chain");
# TLS verification is NOT switched off for it.  xrdcp (XRootD) is tried last
# if it is installed (no resume: it restarts).
# The cloud session cannot reach opendata.cern.ch at all: run this in WSL2.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PY="${PYTHON:-python}"
NAME=Run2012BC_DoubleMuParked_Muons.root
EOSPATH=/eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/$NAME
URL="https://opendata.cern.ch$EOSPATH"
XRD="root://eospublic.cern.ch/$EOSPATH"
OUT=data/raw/opendata/$NAME
PROV=data/raw/cms_2012_doublemuparked.json
mkdir -p data/raw/opendata

size_of() { stat -c %s "$1" 2>/dev/null || echo 0; }
remote=$(curl -sfIL --max-time 60 "$URL" | tr -d '\r' | awk 'tolower($1)=="content-length:"{s=$2} END{print s}')
if [[ -z $remote ]] || (( remote < 1000000 )); then
    echo "could not get the file size from $URL (HTTP error or blocked network: '${remote:-none}')"; exit 1
fi
echo "### $NAME: server size $remote bytes ($(awk "BEGIN{printf \"%.2f\", $remote/1e9}") GB), have $(size_of "$OUT")"

source_used="$URL"
for attempt in $(seq 1 50); do
    have=$(size_of "$OUT")
    [[ $have == "$remote" ]] && break
    if (( have > remote )); then echo "local file larger than the server's -- removing it"; rm -f "$OUT"; fi
    echo "### attempt $attempt: resuming at $(awk "BEGIN{printf \"%.1f\", $have/1e6}") MB"
    curl -fL --retry 5 --retry-delay 10 --speed-limit 1000 --speed-time 120 -C - -o "$OUT" "$URL"; rc=$?
    if (( rc == 22 )); then echo "HTTP error from the server -- not retrying"; break; fi
    (( rc != 0 )) && sleep 15
done

if [[ $(size_of "$OUT") != "$remote" ]] && command -v xrdcp >/dev/null 2>&1; then
    echo "### HTTPS did not finish; trying $XRD (restarts from zero)"
    xrdcp -f "$XRD" "$OUT.xrd" && mv "$OUT.xrd" "$OUT" && source_used="$XRD"
fi
if [[ $(size_of "$OUT") != "$remote" ]]; then
    echo "INCOMPLETE: $(size_of "$OUT") of $remote bytes.  Run the script again to resume, or download"
    echo "$NAME from https://opendata.cern.ch/record/12341 by browser and put it at $OUT."
    echo "No provenance written."
    exit 1
fi

"$PY" - "$OUT" "$source_used" "$remote" "$PROV" <<'EOF' || exit 1
import datetime as dt, hashlib, json, pathlib, sys
path, source, remote, prov_path = pathlib.Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), pathlib.Path(sys.argv[4])
import uproot
try:
    n = uproot.open(path)["Events"].num_entries
except Exception as exc:                       # noqa: BLE001
    sys.exit(f"{path} has the right size but does not open as a ROOT file with 'Events': {exc}")
h = hashlib.sha256()
with open(path, "rb") as f:
    for block in iter(lambda: f.read(1 << 24), b""):
        h.update(block)
prov = {
    "_status": "Downloaded by scripts/m08_fetch_cms_dimuon.sh; complete (size = server Content-Length), "
               "opens with uproot; sha256 computed locally.",
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
    "size_bytes": path.stat().st_size, "server_size_bytes": remote,
    "sha256": h.hexdigest(),
    "n_events": n,
    "downloaded_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "note": "CMS and CERN do not endorse works produced with these data (record terms of use)."}
prov_path.write_text(json.dumps(prov, indent=2) + "\n")
print(f"[prov] {prov_path}: {prov['size_bytes'] / 1e9:.2f} GB, {n:,} events, sha256 {prov['sha256'][:16]}...")
EOF
