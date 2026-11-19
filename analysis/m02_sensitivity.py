#!/usr/bin/env python3
"""Module 2, step 4: how sensitive is the PYTHIA / ALICE dN/deta ratio to the model?

    bash scripts/m02_run_variations.sh            # generates the samples (~5 min)
    python analysis/m02_sensitivity.py name=path.root [name=path.root ...]
           [--data data/raw/alice_pp13TeV_dndeta.json] [--out results/tables/m02/m02_sensitivity.json]

For every sample: <dN/deta> in the data's central bins for INEL>0 events, the
fraction of events that are INEL>0, and MC/data.  The first sample is the
reference ("nominal"); every other sample is also quoted RELATIVE to it, which
is robust against a different PYTHIA version or seed between machines.

Why these variations (scripts/m02_run_variations.sh)
----------------------------------------------------
  cr_off      colour reconnection off -- CR merges strings and is the main
              knob that lowers the multiplicity in Monash
  tune4C      the previous default tune (2010), tuned to earlier LHC data
  pt0ref240   MPI regulator pT0Ref 2.28 -> 2.40 GeV (+5 %); fewer MPI, fewer
              particles.  Monash extrapolates pT0 from 7 to 13 TeV with a power
              law, so this is the parameter an energy-extrapolation error lives in
  decay_all   every unstable particle decays (no c*tau cut): shows how much the
              ALICE "primary" definition is worth in the first place

and one data-side check, done with --inel: the INEL column of the same table.
In |eta| < 1, dN/deta(INEL) = f(INEL>0) x dN/deta(INEL>0) exactly (an event
without a charged primary in |eta| < 1 adds nothing there), so the data
themselves fix f(INEL>0) = ratio of the two columns, and it can be compared
with the generator's fraction directly.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import awkward as ak
import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "analysis"))

from m02_compare_data import load_dataset, read_mc, select_events, central_density   # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="+", help="name=path.root; the first is the reference")
    ap.add_argument("--data", default=str(REPO / "data" / "raw" / "alice_pp13TeV_dndeta.json"))
    ap.add_argument("--out", default=str(REPO / "results" / "tables" / "m02" / "m02_sensitivity.json"))
    ap.add_argument("--pt-data", default=str(REPO / "data" / "raw" / "alice_pp13TeV_pt.json"),
                    help="provenance JSON of the pT spectrum; skipped if absent or unverified")
    args = ap.parse_args()

    prov, d = load_dataset(pathlib.Path(args.data))
    central = (d["xlo"] >= -0.5 - 1e-9) & (d["xhi"] <= 0.5 + 1e-9)
    eta_half = float(max(np.abs(d["xlo"][central]).max(), np.abs(d["xhi"][central]).max()))
    w = (d["xhi"] - d["xlo"])[central]
    data_inel0 = float(np.sum(d["y"][central] * w) / w.sum())
    data_err = float(np.sum(d["sys"][central] * w) / w.sum())

    # the INEL column of the same table fixes the measured INEL>0 fraction
    import yaml
    table = yaml.safe_load((pathlib.Path(args.data).parent / prov["file"]).read_text())
    cols = {i: dv for i, dv in enumerate(table["dependent_variables"])}
    inel_col = next(i for i, dv in cols.items()
                    if any(str(q.get("value")) == "INEL" for q in dv.get("qualifiers", [])))
    y_inel = np.array([v["value"] for v in cols[inel_col]["values"]], float)
    order = np.argsort([v.get("low", v.get("value")) for v in table["independent_variables"][0]["values"]])
    y_inel = y_inel[order]
    data_inel = float(np.sum(y_inel[central] * w) / w.sum())
    f_data = data_inel / data_inel0

    # optional: the pT spectrum (Table 2) -- <pT> and the ratio in pT regions
    pt = None
    try:
        from m02_compare_pt import invariant_yield
        pprov, pd = load_dataset(pathlib.Path(args.pt_data))
        pe = np.concatenate([pd["xlo"], [pd["xhi"][-1]]])
        pc = 0.5 * (pe[1:] + pe[:-1])

        def mean_pt(y):
            wgt = y * 2 * np.pi * pc * np.diff(pe)
            return float((wgt * pc).sum() / wgt.sum())
        pt = {"edges": pe, "ctr": pc, "data": pd["y"], "data_mean_pt": mean_pt(pd["y"]),
              "regions": [(0.15, 0.25), (0.4, 0.8), (1.0, 4.0)], "mean_pt": mean_pt,
              "yield": invariant_yield}
    except (ValueError, FileNotFoundError) as exc:
        print(f"[pt] pT block skipped: {exc}")

    rows = []
    for spec in args.samples:
        name, path = spec.split("=", 1)
        mc = read_mc(pathlib.Path(path))
        mask = select_events(mc, "INEL>0")
        c, ce = central_density(mc, mask, 0.0, eta_half=eta_half)
        f = float(mask.mean())
        rows.append({"name": name, "file": path, "settings": mc["provenance"].get("extra_settings", []),
                     "pythia_version": mc["provenance"].get("pythia_version"),
                     "dndeta_inel0": c, "dndeta_inel0_err": ce, "ratio_inel0": c / data_inel0,
                     "f_inel0": f, "dndeta_inel": c * f, "ratio_inel": c * f / data_inel})
        if pt is not None:
            y, _ = pt["yield"](mc, mask, pt["edges"])
            r = y / pt["data"]
            rows[-1]["mean_pt"] = pt["mean_pt"](y)
            rows[-1]["mean_pt_ratio"] = rows[-1]["mean_pt"] / pt["data_mean_pt"]
            rows[-1]["pt_region_ratio"] = {f"{lo}-{hi}": float(r[(pt["ctr"] > lo) & (pt["ctr"] < hi)].mean())
                                           for lo, hi in pt["regions"]}
    ref = rows[0]
    for r in rows:
        r["rel_to_ref"] = r["dndeta_inel0"] / ref["dndeta_inel0"]

    print(f"data (ALICE, |eta| < {eta_half:g}): INEL>0 {data_inel0:.3f} +- {data_err:.3f}, "
          f"INEL {data_inel:.3f}  ->  f(INEL>0) = {f_data:.3f}")
    print(f"{'sample':12s} {'dN/deta':>14s} {'MC/data':>8s} {'vs ref':>7s} "
          f"{'f(INEL>0)':>9s} {'INEL MC/data':>12s}  settings")
    for r in rows:
        print(f"{r['name']:12s} {r['dndeta_inel0']:7.3f} ± {r['dndeta_inel0_err']:.3f} "
              f"{r['ratio_inel0']:8.3f} {r['rel_to_ref']:7.3f} {r['f_inel0']:9.3f} "
              f"{r['ratio_inel']:12.3f}  {'; '.join(r['settings']) or '(card)'}")

    if pt is not None:
        print(f"\npT spectrum (|eta| < 0.8, INEL>0, 0.15-20 GeV/c): ALICE <pT> = {pt['data_mean_pt']:.3f} GeV/c")
        print(f"{'sample':12s} {'<pT>':>7s} {'MC/data':>8s}  " + "  ".join(f"{lo}-{hi}".rjust(9) for lo, hi in pt["regions"]))
        for r in rows:
            print(f"{r['name']:12s} {r['mean_pt']:7.3f} {r['mean_pt_ratio']:8.3f}  "
                  + "  ".join(f"{v:9.3f}" for v in r["pt_region_ratio"].values()))

    out = pathlib.Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"eta_range": [-eta_half, eta_half],
                               "data": {"dndeta_inel0": data_inel0, "dndeta_inel0_sys": data_err,
                                        "dndeta_inel": data_inel, "f_inel0": f_data,
                                        "source": prov["hepdata"]},
                               "data_pt": ({"mean_pt": pt["data_mean_pt"], "source": pprov["hepdata"]}
                                           if pt is not None else None),
                               "samples": rows}, indent=2))
    print(f"[out] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
