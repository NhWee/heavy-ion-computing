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
import urllib.parse
import urllib.request

HEPDATA = "https://www.hepdata.net"


def get(url: str, timeout=60) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "heavy-ion-computing/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("provenance", help="data/raw/<dataset>.json")
    ap.add_argument("--table", default=None, help='table name, e.g. "Table 1"')
    ap.add_argument("--version", type=int, default=None, help="record version (default: latest)")
    args = ap.parse_args()

    pj = pathlib.Path(args.provenance)
    prov = json.loads(pj.read_text())
    record = prov["hepdata"]["record"]
    table = args.table or prov["hepdata"]["table"]

    # 1. record metadata ---------------------------------------------------
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

    # 2. tables --------------------------------------------------------------
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

    # 3. download YAML -------------------------------------------------------
    turl = f"{HEPDATA}/download/table/{record}/{urllib.parse.quote(table)}/{version}/yaml"
    print(f"[get] {turl}")
    payload = get(turl)
    out = pj.parent / prov["file"]
    out.write_bytes(payload)
    sha = hashlib.sha256(payload).hexdigest()
    print(f"[out] {out} ({len(payload)} bytes, sha256 {sha[:16]}...)")

    # 4. show the header so the event class can be confirmed by eye -------
    import yaml
    doc = yaml.safe_load(payload)
    for dv in doc.get("dependent_variables", []):
        print(f"      dependent: {dv.get('header')}  qualifiers: {dv.get('qualifiers')}")
    for iv in doc.get("independent_variables", []):
        vals = iv.get("values", [])
        lo = min(v.get("low", v.get("value")) for v in vals)
        hi = max(v.get("high", v.get("value")) for v in vals)
        print(f"      independent: {iv.get('header')}  range [{lo}, {hi}]  ({len(vals)} bins)")
    n_dep = len(doc.get("dependent_variables", []))
    if n_dep > 1:
        print(f"NOTE: {n_dep} dependent variables in this table (e.g. INEL and INEL>0). "
              "analysis/m02_compare_data.py reads the FIRST; add "
              '"dependent_index" to the JSON if the wanted class is another one.')

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
    return 0


if __name__ == "__main__":
    sys.exit(main())
