#!/usr/bin/env python3
"""Module 2, step 7: PYTHIA vs ALICE identified pi, K, p spectra (pp 13 TeV, INEL>0).

    python analysis/m02_compare_identified.py [--mc data/generated/m02_pythia_mb_13TeV.root]
           [--data pi=data/raw/alice_pp13TeV_pi.json K=... p=...] [--set-label "Monash"]

Why
---
The charged-particle pT spectrum (m02_compare_pt.py) is +12-14 % high at
0.4-0.8 GeV/c whatever MPI, colour reconnection or string pT width does.
A charged spectrum is the SUM of pions, kaons and protons, and that pT range
is exactly where their mass-dependent shapes overlap (Module 1's closure
lesson).  The identified spectra (ALICE, arXiv:2003.02394, HEPData ins1784041,
|y| < 0.5) tell WHICH species carries the excess: a pion excess is a bulk
(MPI/fragmentation) problem, a kaon/proton excess a flavour- or baryon-
production problem.

Details that differ from the charged-particle comparison
  - RAPIDITY y = 0.5 ln((E + pz)/(E - pz)), not pseudorapidity: identified
    spectra are measured in |y| < 0.5.  The MC uses each particle's own E.
  - the observable form is read from the HEPData header: "1/(2 pi pT)" ->
    invariant yield, otherwise d2N/(dy dpT); "(pi+ + pi-)/2" or a qualifier
    with "/2" -> average of the charges, otherwise the sum.  What was detected
    is printed and recorded, because a factor 2 here would look like physics.
  - integrated yields over the MEASURED pT range only (no extrapolation), and
    the ratios K/pi and p/pi, which cancel the overall MPI normalisation.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import awkward as ak
import numpy as np
import uproot
import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from hicpy.style import apply_style, save_figure, OKABE_ITO, plt   # noqa: E402
from m02_compare_data import load_dataset                          # noqa: E402

PDG = {"pi": 211, "K": 321, "p": 2212}
Y_MAX = 0.5


def detect_form(prov_json: pathlib.Path):
    """Read the table header/qualifiers: invariant or not, sum or average."""
    prov = json.loads(prov_json.read_text())
    doc = yaml.safe_load((prov_json.parent / prov["file"]).read_text())
    dv = doc["dependent_variables"][int(prov.get("dependent_index", 0))]
    text = (str(dv.get("header", {})) + " " + " ".join(f"{q.get('name')}={q.get('value')}"
                                                       for q in dv.get("qualifiers", []))).upper()
    invariant = "2*PI*PT" in text.replace(" ", "") or "1/PT" in text.replace(" ", "")
    average = "/2" in text.replace(" ", "") or "AVERAGE" in text
    return {"invariant": invariant, "average": average, "header_text": text[:300]}


def read_mc(path):
    with uproot.open(path) as f:
        ev = f["events"].arrays(["particle_px", "particle_py", "particle_pz", "particle_E",
                                 "particle_pdg", "particle_charge"])
        prov = json.loads(str(f["provenance"]))
    f64 = lambda a: ak.values_astype(a, np.float64)                 # noqa: E731
    px, py, pz, e = map(f64, (ev.particle_px, ev.particle_py, ev.particle_pz, ev.particle_E))
    pt = np.hypot(px, py); p = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
    eta = np.arctanh((pz / p) * (1 - 1e-12))
    y = 0.5 * np.log((e + pz) / (e - pz))
    charged = ev.particle_charge != 0
    inel0 = ak.to_numpy(ak.sum(charged & (abs(eta) < 1.0), axis=1)) > 0
    return {"pt": pt, "y": y, "apdg": abs(ev.particle_pdg), "inel0": inel0, "provenance": prov}


def mc_spectrum(mc, species, edges, form):
    """Spectrum per INEL>0 event in |y| < Y_MAX, in the data's form."""
    m = mc["inel0"]
    sel = (mc["apdg"] == PDG[species]) & (abs(mc["y"]) < Y_MAX)
    pt = ak.to_numpy(ak.flatten(mc["pt"][m][sel[m]]))
    w = np.ones_like(pt)
    if form["invariant"]:
        w = w / (2 * np.pi * pt)
    if form["average"]:
        w = w / 2.0
    sumw, _ = np.histogram(pt, edges, weights=w)
    sumw2, _ = np.histogram(pt, edges, weights=w ** 2)
    norm = m.sum() * np.diff(edges) * (2 * Y_MAX)
    return sumw / norm, np.sqrt(sumw2) / norm


def integrate(y, edges, invariant):
    c = 0.5 * (edges[1:] + edges[:-1]); jac = 2 * np.pi * c if invariant else 1.0
    return float(np.sum(y * jac * np.diff(edges)))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mc", default=str(REPO / "data" / "generated" / "m02_pythia_mb_13TeV.root"))
    ap.add_argument("--data", nargs="+", default=[f"{s}=data/raw/alice_pp13TeV_{s.lower()}.json"
                                                  for s in ("pi", "K", "p")])
    ap.add_argument("--label", default="PYTHIA 8")
    ap.add_argument("--figdir", default=str(REPO / "results" / "figures" / "m02"))
    ap.add_argument("--tabdir", default=str(REPO / "results" / "tables" / "m02"))
    ap.add_argument("--tag", default="", help="suffix for output names (variations)")
    args = ap.parse_args()
    figdir, tabdir = pathlib.Path(args.figdir), pathlib.Path(args.tabdir)
    figdir.mkdir(parents=True, exist_ok=True); tabdir.mkdir(parents=True, exist_ok=True)

    mc = read_mc(args.mc)
    print(f"[mc]   {pathlib.Path(args.mc).name}: {len(mc['inel0'])} events, {int(mc['inel0'].sum())} INEL>0")
    res = {}
    for spec in args.data:
        sp, path = spec.split("=", 1)
        pj = pathlib.Path(path) if pathlib.Path(path).is_absolute() else REPO / path
        prov, d = load_dataset(pj)
        form = detect_form(pj)
        edges = np.concatenate([d["xlo"], [d["xhi"][-1]]])
        y, e = mc_spectrum(mc, sp, edges, form)
        r = y / d["y"]; tot = np.hypot(d["stat"], d["sys"])
        r_err = np.hypot(e / d["y"], r * tot / d["y"])
        res[sp] = {"edges": edges.tolist(), "mc": y.tolist(), "mc_err": e.tolist(), "data": d["y"].tolist(),
                   "data_tot": tot.tolist(), "ratio": r.tolist(), "ratio_err": r_err.tolist(), "form": form,
                   "yield_mc": integrate(y, edges, form["invariant"]),
                   "yield_data": integrate(d["y"], edges, form["invariant"]),
                   "source": prov["hepdata"]}
        c = 0.5 * (edges[1:] + edges[:-1]); k = (c > 0.4) & (c < 0.8)
        res[sp]["ratio_0.4-0.8"] = float(r[k].mean()) if k.any() else None
        print(f"[{sp:2s}]  form: {'invariant' if form['invariant'] else 'dN/dydpT'}, "
              f"{'average of charges' if form['average'] else 'sum of charges'}; "
              f"pT {edges[0]:g}-{edges[-1]:g}, {len(d['y'])} bins; yield MC/data "
              f"{res[sp]['yield_mc'] / res[sp]['yield_data']:.3f}; MC/data at 0.4-0.8 GeV/c "
              f"{res[sp]['ratio_0.4-0.8'] if res[sp]['ratio_0.4-0.8'] is not None else float('nan'):.3f}")
    if "pi" in res:
        for sp in ("K", "p"):
            if sp in res:
                rmc = res[sp]["yield_mc"] / res["pi"]["yield_mc"]
                rd = res[sp]["yield_data"] / res["pi"]["yield_data"]
                res[f"{sp}/pi"] = {"mc": rmc, "data": rd, "double_ratio": rmc / rd}
                print(f"[{sp}/pi] MC {rmc:.4f} vs data {rd:.4f} -> {rmc / rd:.3f} "
                      f"(measured pT ranges only; ranges differ per species)")

    apply_style()
    fig, axes = plt.subplots(1, len([s for s in ("pi", "K", "p") if s in res]), figsize=(12.6, 4.0), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, sp in zip(axes, [s for s in ("pi", "K", "p") if s in res]):
        r = res[sp]; e = np.array(r["edges"]); c = 0.5 * (e[1:] + e[:-1])
        ax.fill_between(c, 1 - np.array(r["data_tot"]) / np.array(r["data"]),
                        1 + np.array(r["data_tot"]) / np.array(r["data"]), color="0.88", step="mid", lw=0)
        ax.errorbar(c, r["ratio"], yerr=r["ratio_err"], fmt="o", ms=4, color=OKABE_ITO[1])
        ax.axhline(1, color="0.4", lw=1); ax.axvspan(0.4, 0.8, color=OKABE_ITO[4], alpha=0.12, lw=0)
        ax.set_xscale("log"); ax.set_xlabel(r"$p_{\rm T}$ (GeV/$c$)")
        ax.set_title({"pi": r"$\pi^\pm$", "K": r"K$^\pm$", "p": r"p, $\bar{\rm p}$"}[sp], fontsize=13)
        ax.set_ylim(0.5, 1.6)
    axes[0].set_ylabel(f"{args.label} / ALICE")
    save_figure(fig, figdir / f"m02_identified{args.tag}"); plt.close(fig)
    out = tabdir / f"m02_identified{args.tag}_summary.json"
    out.write_text(json.dumps({"mc_file": pathlib.Path(args.mc).name,
                               "pythia_version": mc["provenance"].get("pythia_version"),
                               "extra_settings": mc["provenance"].get("extra_settings"),
                               "y_max": Y_MAX, "species": res}, indent=1))
    print(f"[out] {figdir / f'm02_identified{args.tag}.pdf'}, {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
