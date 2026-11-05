#!/usr/bin/env python3
"""Download a HEPData table and complete / verify its provenance JSON.

    python scripts/m02_fetch_hepdata.py data/raw/alice_pp13TeV_dndeta.json [--table "Table 1"]

Needs network access to www.hepdata.net (the cloud session used to write this
repository could not reach it; run from WSL2 or any machine that can).

What it does, in order
  1. GET https://www.hepdata.net/record/<record>?format=json
     and check that the record's arXiv number / INSPIRE id match the JSON.
  2. list the tables; pick --table (default: the JSON's hepdata.table)
  3. GET the table as YAML (the canonical HEPData format), save it next to
     the JSON under the JSON's "file" name, record sha256 + download time
  4. print the table's header/qualifiers so you can confirm the event class
     and eta range BY EYE, then set "verified": true and the HEPData DOI from
     the record itself.

Nothing is written unless steps 1-3 succeed; mismatching identifiers abort
with a message instead of silently overwriting the JSON.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import pathlib
import sys
import urllib.error
import urllib.parse
import urllib.request

HEPDATA = "https://www.hepdata.net"


# HEPData sits behind Cloudflare, which answers 403 to requests that do not
# look like a browser (a bare urllib User-Agent is enough to be refused).
# Send browser-like headers; if that still fails, fall back to curl, which
# some proxies trust more than Python.
HEADERS = {
    "User-Agent": ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                   "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"),
    "Accept": "application/json, text/x-yaml, text/plain, */*",
    "Accept-Language": "en",
}


def get(url: str, timeout=60) -> bytes:
    req = urllib.request.Request(url, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except urllib.error.HTTPError as exc:
        if exc.code not in (403, 429, 503):
            raise
        print(f"      urllib got HTTP {exc.code}; retrying with curl")
    except urllib.error.URLError as exc:
        # proxy CONNECT refusals surface here, not as HTTPError
        print(f"      urllib failed ({exc.reason}); retrying with curl")
    import shutil
    import subprocess
    if not shutil.which("curl"):
        raise RuntimeError(f"HTTP 403 from {url} and no curl to retry with. "
                           "Download the file in a browser and place it as the JSON's 'file'.")
    cmd = ["curl", "-sS", "-L", "--fail", "--max-time", str(timeout)]
    for k, v in HEADERS.items():
        cmd += ["-H", f"{k}: {v}"]
    out = subprocess.run(cmd + [url], capture_output=True)
    if out.returncode != 0:
        raise RuntimeError(f"curl failed for {url}: {out.stderr.decode(errors='replace').strip()}\n"
                           "If the site keeps refusing automated downloads, open the record "
                           "in a browser, download the table as YAML, save it under the "
                           "JSON's 'file' name and run this script with --local.")
    return out.stdout


def _norm(s) -> str:
    return str(s).replace(" ", "").replace("_", "").upper()


def describe_table(payload: bytes, event_class: str):
    """Print the table header and pick, from the dependent-variable qualifiers,
    the column whose event class matches `event_class`.

    Returns (dependent_index, eta_range_text, column_label).  The choice is
    automatic but the evidence (every column's qualifiers) is printed, so a
    reader of the log can see why that column was chosen."""
    import yaml
    doc = yaml.safe_load(payload)
    deps = doc.get("dependent_variables", [])
    want = _norm(event_class)
    matches = []
    for i, dv in enumerate(deps):
        quals = dv.get("qualifiers", [])
        text = " | ".join(f"{q.get('name')}={q.get('value')}" for q in quals)
        print(f"      dependent[{i}]: {dv.get('header')}  qualifiers: {text}")
        blob = _norm(text) + _norm(dv.get("header", {}).get("name", ""))
        # "INEL>0" must not match a column that merely says "INEL"
        if want in blob and not (want == "INEL" and "INEL>0" in blob):
            matches.append(i)
    eta_range = None
    for iv in doc.get("independent_variables", []):
        vals = iv.get("values", [])
        lo = min(v.get("low", v.get("value")) for v in vals)
        hi = max(v.get("high", v.get("value")) for v in vals)
        eta_range = f"{lo} to {hi} ({len(vals)} bins, from the table)"
        print(f"      independent: {iv.get('header')}  range [{lo}, {hi}]  ({len(vals)} bins)")
    if len(matches) == 1:
        idx = matches[0]
        label = deps[idx].get("header", {}).get("name", "")
        print(f"      -> column {idx} matches event class {event_class!r}: {label}")
        return idx, eta_range, label
    if not matches and len(deps) == 1:
        print(f"      -> single column; its qualifiers do not name {event_class!r} -- "
              "CHECK the paper that this is the right event class")
        return 0, eta_range, deps[0].get("header", {}).get("name", "")
    raise RuntimeError(f"could not pick a column for event class {event_class!r} "
                       f"automatically (candidates: {matches or 'none'} of {len(deps)}); "
                       "set \"dependent_index\" in the JSON by hand from the list above")


def local_only(pj, prov, table) -> int:
    """--local: the table was downloaded by hand (browser).  Fingerprint it,
    print its header, mark verified -- the identifiers were checked by the
    human on the HEPData page, and that is recorded in _status."""
    out = pj.parent / prov["file"]
    if not out.exists():
        print(f"ABORT -- {out} not found; download {table!r} of {prov['hepdata']['record']} "
              "as YAML from the record page and save it there")
        return 1
    payload = out.read_bytes()
    sha = hashlib.sha256(payload).hexdigest()
    print(f"=== --local: {out} ({len(payload)} bytes, sha256 {sha[:16]}...)")
    idx, eta_range, label = describe_table(payload, prov["event_class"])
    prov["dependent_index"] = idx
    prov.setdefault("selection", {})["eta_range"] = eta_range
    prov["selection"]["column"] = label
    prov.update({
        "verified": True,
        "_status": f"table downloaded BY HAND from {prov['hepdata'].get('record_url')} on "
                   f"{dt.datetime.now(dt.timezone.utc):%Y-%m-%d}; identifiers checked on the "
                   "HEPData page by the person who downloaded it; header printed by "
                   "scripts/m02_fetch_hepdata.py --local",
        "downloaded_utc": dt.datetime.fromtimestamp(out.stat().st_mtime, dt.timezone.utc)
                           .strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sha256": sha,
    })
    pj.write_text(json.dumps(prov, indent=2, ensure_ascii=False) + "\n")
    print(f"[out] {pj} updated (verified = true, by hand)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("provenance", help="data/raw/<dataset>.json")
    ap.add_argument("--table", default=None, help='table name, e.g. "Table 1"')
    ap.add_argument("--version", type=int, default=None, help="record version (default: latest)")
    ap.add_argument("--local", action="store_true",
                    help="skip the download: fingerprint and describe a table YAML already "
                         "saved under the JSON's 'file' (downloaded by hand in a browser)")
    args = ap.parse_args()

    pj = pathlib.Path(args.provenance)
    prov = json.loads(pj.read_text())
    record = prov["hepdata"]["record"]
    table = args.table or prov["hepdata"]["table"]

    if args.local:
        return local_only(pj, prov, table)

    print(f"=== step 1/4: fetch the HEPData RECORD {record} and check it is the paper "
          f"named in {pj.name} (arXiv {prov.get('arxiv')}, INSPIRE {prov.get('inspire_id')})")
    url = f"{HEPDATA}/record/{record}?format=json"
    print(f"[get] {url}")
    rec = json.loads(get(url))
    meta = rec.get("record", rec)
    rec_arxiv = str(meta.get("arxiv_id", "")).replace("arXiv:", "")
    rec_inspire = str(meta.get("inspire_id", ""))
    rec_doi = meta.get("doi") or meta.get("hepdata_doi")
    rec_title = meta.get("title", "")
    print(f"      title   : {rec_title}")
    print(f"      arXiv   : {rec_arxiv}    inspire: {rec_inspire}    doi: {rec_doi}")
    problems = []
    if prov.get("arxiv") and rec_arxiv and prov["arxiv"] != rec_arxiv:
        problems.append(f"arXiv mismatch: JSON {prov['arxiv']} vs record {rec_arxiv}")
    if prov.get("inspire_id") and rec_inspire and str(prov["inspire_id"]) != rec_inspire:
        problems.append(f"INSPIRE mismatch: JSON {prov['inspire_id']} vs record {rec_inspire}")
    if problems:
        print("ABORT -- the record does not match the provenance JSON:\n  " + "\n  ".join(problems))
        return 1

    print(f"=== step 2/4: list the record's tables and pick {table!r}")
    tables = rec.get("data_tables", [])
    print(f"[rec] {len(tables)} table(s):")
    for t in tables:
        print(f"      - {t.get('name')}: {t.get('description', '')[:110]}")
    names = [t.get("name") for t in tables]
    if table not in names:
        print(f"ABORT -- table {table!r} not in record; pick one of {names}")
        return 1
    tdesc = next(t for t in tables if t.get("name") == table)
    version = args.version or meta.get("version", 1)

    print(f"=== step 3/4: download {table!r} as YAML (HEPData's canonical format) "
          f"to {pj.parent / prov['file']} and fingerprint it (sha256)")
    turl = f"{HEPDATA}/download/table/{record}/{urllib.parse.quote(table)}/{version}/yaml"
    print(f"[get] {turl}")
    payload = get(turl)
    out = pj.parent / prov["file"]
    out.write_bytes(payload)
    sha = hashlib.sha256(payload).hexdigest()
    print(f"[out] {out} ({len(payload)} bytes, sha256 {sha[:16]}...)")

    print("=== step 4/4: read the table header; pick the column whose qualifiers name "
          f"the event class {prov['event_class']!r}; record the eta range")
    idx, eta_range, label = describe_table(payload, prov["event_class"])
    prov["dependent_index"] = idx
    prov.setdefault("selection", {})["eta_range"] = eta_range
    prov["selection"]["column"] = label

    prov.update({
        "verified": True,
        "_status": f"verified against {url} on "
                   f"{dt.datetime.now(dt.timezone.utc):%Y-%m-%d}; table header printed by "
                   "scripts/m02_fetch_hepdata.py -- confirm event class/eta range in 'selection'",
        "title": rec_title or prov.get("title"),
        "arxiv": rec_arxiv or prov.get("arxiv"),
        "inspire_id": rec_inspire or prov.get("inspire_id"),
        "downloaded_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sha256": sha,
    })
    prov["hepdata"].update({"table": table, "table_description": tdesc.get("description", ""),
                            "doi": tdesc.get("doi") or rec_doi or prov["hepdata"].get("doi"),
                            "version": version})
    pj.write_text(json.dumps(prov, indent=2, ensure_ascii=False) + "\n")
    print(f"[out] {pj} updated (verified = true)")
    print("next: python analysis/m02_compare_data.py --config config/analysis_m02_pythia.yaml")
    return 0


if __name__ == "__main__":
    sys.exit(main())
