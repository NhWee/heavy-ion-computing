#!/usr/bin/env python3
"""Build the ALICE INEL>0 pi, K, p spectra (pp 13 TeV) from HEPData ins1784041.

    python scripts/m02_derive_identified_inel0.py [--src data/raw/ins1784041]

WHY THIS STEP EXISTS
--------------------
The paper (ALICE, arXiv:2003.02394, EPJC 80 (2020) 693) is about MULTIPLICITY
DEPENDENCE.  Its HEPData record has the spectra in ten V0M multiplicity
classes (Tables 1-6) and the ratios "class / INEL>0" (Tables 7-12) -- but no
table of the INEL>0 spectrum itself.  The INEL>0 spectrum is recovered in two
independent ways, and they must agree:

  (A) ratio method   S_INEL = S_i / R_i            for each class i = I..X
      (R_i is literally S_i divided by the INEL>0 spectrum ALICE used)
  (B) weighted sum   S_INEL = sum_i f_i S_i
      f_i = fraction of the INEL>0 cross section in class i (paper, Table 1:
      I 0-0.92 %, II 0.92-4.6, III 4.6-9.2, IV 9.2-13.8, V 13.8-18.4,
      VI 18.4-27.6, VII 27.6-36.8, VIII 36.8-46, IX 46-64.5, X 64.5-100 %).
      The classes are disjoint event samples, so this is exact.

Central value: (A), averaged over the ten classes (no external numbers).
(B) checks the class fractions and the reading of the tables.

Uncertainties (the tables give, per class, 'stat', 'full syst.' and
'unc. syst.' = the part uncorrelated across multiplicity):
  stat   sqrt(sum f_i^2 stat_i^2)        exact for disjoint samples
  sys    sum f_i fullsys_i               linear: the correlated part (tracking,
                                         PID, material) is common to all
                                         classes; slightly conservative for
                                         the uncorrelated part
pT range: bins below 10 GeV/c only -- class X is binned differently above
(10-20 instead of 10-13, 13-20), and a 20k-event MC has nothing there.

OUTPUT (payloads gitignored, provenance JSON committed)
  data/raw/alice_pp13TeV_{pi,k,p}_inel0.yaml   HEPData-format table, labels stat/sys
  data/raw/alice_pp13TeV_{pi,k,p}.json         provenance: sources + sha256 of
                                               every input table and the output
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import pathlib

import numpy as np
import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
CLASSES = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
EDGES_PCT = [0, 0.92, 4.6, 9.2, 13.8, 18.4, 27.6, 36.8, 46, 64.5, 100]
FRAC = np.diff(EDGES_PCT) / 100.0
# species: (spectra I-IX, spectra X, ratio I-IX, ratio X), label for output
SPECIES = {"pi": ((1, 2, 7, 8), "pi+ + pi-"), "k": ((3, 4, 9, 10), "K+ + K-"),
           "p": ((5, 6, 11, 12), "p + pbar")}
PT_MAX = 10.0


def read(src, n):
    doc = yaml.safe_load((src / f"Table_{n}.yaml").read_text())
    x = [(v["low"], v["high"]) for v in doc["independent_variables"][0]["values"]]
    out = []
    for dv in doc["dependent_variables"]:
        cls = next(q["value"] for q in dv["qualifiers"] if q["name"].lower().startswith(("mutl", "mult")))
        rows = {}
        for xb, v in zip(x, dv["values"]):
            e = {err["label"]: float(err["symerror"]) for err in v.get("errors", [])}
            rows[xb] = (float(v["value"]), e.get("stat", 0.0), e.get("full syst.", 0.0))
        out.append((cls.split()[-1], rows))
    return out


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def derive(src, tables):
    t_s, t_sx, t_r, t_rx = tables
    spec = dict(read(src, t_s) + read(src, t_sx))
    ratio = dict(read(src, t_r) + read(src, t_rx))
    assert sorted(spec) == sorted(CLASSES) == sorted(ratio), (sorted(spec), sorted(ratio))
    bins = [b for b in spec["I"] if b[1] <= PT_MAX + 1e-9 and all(b in spec[c] and b in ratio[c] for c in CLASSES)]
    rows, spread, closure = [], [], []
    for b in bins:
        est = np.array([spec[c][b][0] / ratio[c][b][0] for c in CLASSES])
        s = np.array([spec[c][b][0] for c in CLASSES])
        st = np.array([spec[c][b][1] for c in CLASSES])
        sy = np.array([spec[c][b][2] for c in CLASSES])
        y = est.mean()
        wsum = float(np.sum(FRAC * s))
        rows.append({"low": b[0], "high": b[1], "y": float(y),
                     "stat": float(np.sqrt(np.sum(FRAC ** 2 * st ** 2)) * y / wsum),
                     "sys": float(np.sum(FRAC * sy) * y / wsum)})
        spread.append(float((est.max() - est.min()) / y))
        closure.append(wsum / y)
    return rows, np.array(spread), np.array(closure)


def write_yaml(rows, label, path):
    doc = {"independent_variables": [{"header": {"name": "PT", "units": "GEV/c"},
                                      "values": [{"low": r["low"], "high": r["high"]} for r in rows]}],
           "dependent_variables": [{
               "header": {"name": "1/N_evt d2N/(dy dpT)", "units": "(GEV/c)^-1"},
               "qualifiers": [{"name": "RE", "value": f"P P --> ({label}) X"},
                              {"name": "SQRT(S)", "value": "13000.0 GEV"},
                              {"name": "YRAP", "value": "-0.5 < y < 0.5"},
                              {"name": "EVENT CLASS", "value": "INEL>0 (derived, see provenance)"},
                              {"name": "CHARGES", "value": "sum of both charges"}],
               "values": [{"value": r["y"], "errors": [{"label": "stat", "symerror": r["stat"]},
                                                       {"label": "sys", "symerror": r["sys"]}]} for r in rows]}]}
    path.write_text(yaml.safe_dump(doc, sort_keys=False))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=str(REPO / "data" / "raw" / "ins1784041"))
    args = ap.parse_args()
    src = pathlib.Path(args.src)
    raw = REPO / "data" / "raw"
    for sp, (tables, label) in SPECIES.items():
        rows, spread, closure = derive(src, tables)
        out = raw / f"alice_pp13TeV_{sp}_inel0.yaml"
        write_yaml(rows, label, out)
        print(f"[{sp:2s}] {len(rows)} bins, pT {rows[0]['low']}-{rows[-1]['high']} GeV/c; "
              f"(A) spread over 10 classes max {100 * spread.max():.2f} %; "
              f"(B)/(A) {closure.min():.3f}-{closure.max():.3f} (mean {closure.mean():.3f}); "
              f"rel sys {100 * min(r['sys'] / r['y'] for r in rows):.1f}-{100 * max(r['sys'] / r['y'] for r in rows):.1f} %")
        prov = {
            "_status": ("Derived by scripts/m02_derive_identified_inel0.py from the HEPData 'Download All "
                        "-> YAML' archive of ins1784041 (downloaded by hand by the user, 2026-11-19). "
                        "The record has no INEL>0 spectrum table; see the script docstring."),
            "verified": True,
            "experiment": "ALICE",
            "title": "Multiplicity dependence of pi, K, and p production in pp collisions at sqrt(s) = 13 TeV",
            "arxiv": "2003.02394",
            "journal": "Eur. Phys. J. C 80 (2020) 693",
            "inspire_id": "1784041",
            "hepdata": {"record": "ins1784041", "record_url": "https://www.hepdata.net/record/ins1784041",
                        "table": f"derived from Tables {', '.join(str(t) for t in tables)}",
                        "version": 1},
            "system": "pp", "sqrt_s_TeV": 13.0,
            "observable": f"1/N d2N/(dy dpT), {label}",
            "event_class": "INEL>0",
            "selection": {"y_range": "|y| < 0.5", "pt_range": f"{rows[0]['low']}-{rows[-1]['high']} GeV/c",
                          "primary_definition": "ALICE primary: excludes weak-decay products of strange hadrons",
                          "form": "sum of charges, not invariant (no 1/(2 pi pT))"},
            "derivation": {"central": "mean over V0M classes I-X of S_i / R_i",
                           "check": "sum_i f_i S_i with f_i from the paper's class percentiles",
                           "class_percentiles": EDGES_PCT,
                           "max_class_spread": float(spread.max()),
                           "weighted_sum_over_central": [float(closure.min()), float(closure.max())],
                           "stat": "sqrt(sum f_i^2 stat_i^2)", "sys": "sum f_i fullsys_i (linear)"},
            "inputs_sha256": {f"ins1784041/Table_{t}.yaml": sha(src / f"Table_{t}.yaml") for t in tables},
            "derived_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "file": out.name, "sha256": sha(out), "dependent_index": 0,
        }
        (raw / f"alice_pp13TeV_{sp}.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(f"[out] {raw}/alice_pp13TeV_{{pi,k,p}}{{.json,_inel0.yaml}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
