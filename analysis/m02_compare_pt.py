#!/usr/bin/env python3
"""Module 2, second observable: PYTHIA pT spectrum vs ALICE pp 13 TeV (INEL>0).

    python analysis/m02_compare_pt.py [--config config/analysis_m02_pythia.yaml]
           [--mc data/generated/m02_pythia_mb_13TeV.root]
           [--data data/raw/alice_pp13TeV_pt.json] [--figdir ...] [--tabdir ...]

What is compared
----------------
ALICE (arXiv:1509.08734, HEPData ins1395253 Table 2) publishes the invariant
charged-particle yield

    1/(2 pi pT) 1/N_ev d^2N/(d eta d pT)      in |eta| < 0.8, INEL>0 events,

for primary charged particles with 0.15 < pT < 20 GeV/c.  The MC side applies
exactly that: the INEL>0 event class computed from the particles (same function
as m02_compare_data.py), charge != 0, |eta| < 0.8, the 1/(2 pi pT) Jacobian per
particle (troubleshooting #7), normalised per selected event and per bin width
and per unit eta, in the data's own pT bins.

Why this observable after dN/deta
---------------------------------
dN/deta fixes the NUMBER of particles; the pT spectrum says how that number is
shared between soft (bulk) and hard (jet-like) production.  If the 6 % dN/deta
excess is a normalisation of the soft component, the MC/data ratio is high at
low pT and closer to 1 at high pT; if PYTHIA is uniformly high, the ratio is
flat.  Integrated over pT the ratio must reproduce the dN/deta result, which
is printed as a cross-check (the integral covers only 0.15 < pT < 20 GeV/c).

The measurement is read only through its provenance JSON (verified = true), as
in m02_compare_data.py.  Without it the script runs MC-only.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import awkward as ak
import numpy as np
import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from hicpy.style import apply_style, stamp, save_figure, OKABE_ITO, plt   # noqa: E402
from m01_analysis_uproot import rel, hist_record                          # noqa: E402
from m02_compare_data import load_dataset, read_mc, select_events        # noqa: E402

ETA_MAX = 0.8


def invariant_yield(mc, mask, edges, eta_max=ETA_MAX):
    """1/(2 pi pT) 1/N_ev d^2N/(deta dpT) per bin, per-particle Jacobian."""
    sel = mc["charged"] & (abs(mc["eta"]) < eta_max)
    pt = mc["pt"][mask][sel[mask]]
    w_ev = mc["weight"][mask]
    w_ev_part = ak.to_numpy(ak.flatten(ak.broadcast_arrays(w_ev, pt)[0]))
    pt = ak.to_numpy(ak.flatten(pt))
    w = w_ev_part / (2.0 * np.pi * pt * 2.0 * eta_max)
    sumw, _ = np.histogram(pt, bins=edges, weights=w)
    sumw2, _ = np.histogram(pt, bins=edges, weights=w ** 2)
    norm = w_ev.sum() * np.diff(edges)
    return sumw / norm, np.sqrt(sumw2) / norm


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(REPO / "config" / "analysis_m02_pythia.yaml"))
    ap.add_argument("--mc", default=None)
    ap.add_argument("--data", default=None)
    ap.add_argument("--no-data", action="store_true")
    ap.add_argument("--figdir", default=None)
    ap.add_argument("--tabdir", default=None)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    cmp_cfg = cfg.get("comparison", {})
    mc_path = pathlib.Path(args.mc) if args.mc else REPO / cfg["paths"]["generated_file"]
    figdir = pathlib.Path(args.figdir) if args.figdir else REPO / cfg["paths"]["figure_dir"]
    tabdir = pathlib.Path(args.tabdir) if args.tabdir else REPO / cfg["paths"]["table_dir"]
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    data, dprov = None, None
    if not args.no_data:
        djson = pathlib.Path(args.data) if args.data else (
            REPO / cmp_cfg["pt_data"] if "pt_data" in cmp_cfg else None)
        if djson is None or not djson.exists():
            print(f"[data] no provenance JSON ({djson}); running MC only")
        else:
            try:
                dprov, data = load_dataset(djson)
                print(f"[data] {dprov['experiment']} {dprov['observable']} "
                      f"{dprov['event_class']}, {len(data['y'])} bins, from {djson.name}")
            except (ValueError, FileNotFoundError) as exc:
                print(f"[data] {exc}\n[data] running MC only")
                data = None

    mc = read_mc(mc_path)
    mask = select_events(mc, "INEL>0")
    print(f"[mc]   {mc_path.name}: {len(mc['weight'])} events, {int(mask.sum())} INEL>0")

    edges = (np.concatenate([data["xlo"], [data["xhi"][-1]]]) if data is not None
             else np.geomspace(0.15, 20.0, 41))
    y, e = invariant_yield(mc, mask, edges)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    dpt = np.diff(edges)

    # cross-check: integrating the invariant yield back to dN/deta (pT in range)
    dndeta_mc = float(np.sum(y * 2 * np.pi * ctr * dpt))
    print(f"[mc]   integral -> dN/deta(|eta|<{ETA_MAX}, {edges[0]:g}<pT<{edges[-1]:g}) "
          f"= {dndeta_mc:.3f}  (bin-centre Jacobian; cross-check only)")

    summary = {"mc_file": rel(mc_path), "event_class": "INEL>0", "eta_max": ETA_MAX,
               "yield_mc": hist_record(edges, y, e), "dndeta_from_yield_mc": dndeta_mc,
               "mc_provenance": mc["provenance"], "data": None, "ratio": None}

    apply_style()
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(6.2, 6.2), sharex=True,
                                  gridspec_kw={"height_ratios": [2.4, 1], "hspace": 0.06})
    ax.errorbar(ctr, y, yerr=e, fmt="s", ms=4, mfc="none", color=OKABE_ITO[1],
                label="PYTHIA 8 (INEL>0)")
    if data is not None:
        tot = np.hypot(data["stat"], data["sys"])
        r = y / data["y"]
        r_err = np.hypot(e / data["y"], r * tot / data["y"])
        ax.errorbar(data["x"], data["y"], yerr=tot, fmt="o", ms=4, color=OKABE_ITO[0],
                    label=f"{dprov['experiment']} ({dprov['event_class']})")
        axr.errorbar(ctr, r, yerr=r_err, fmt="s", ms=4, mfc="none", color=OKABE_ITO[1])
        dndeta_data = float(np.sum(data["y"] * 2 * np.pi * ctr * dpt))
        soft = ctr < 1.0
        summary["data"] = {"provenance": dprov, "y": data["y"].tolist(),
                           "stat": data["stat"].tolist(), "sys": data["sys"].tolist()}
        summary["ratio"] = {"per_bin": r.tolist(), "per_bin_err": r_err.tolist(),
                            "integral": dndeta_mc / dndeta_data,
                            "mean_below_1GeV": float(np.mean(r[soft])),
                            "mean_above_1GeV": float(np.mean(r[~soft]))}
        print(f"[cmp]  MC/data: integral {dndeta_mc / dndeta_data:.3f}; "
              f"pT < 1 GeV/c {np.mean(r[soft]):.3f}; pT > 1 GeV/c {np.mean(r[~soft]):.3f}")
    axr.axhline(1.0, color="0.4", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylabel(r"$\frac{1}{2\pi p_{\rm T}N_{\rm ev}}\frac{d^2N}{d\eta\,dp_{\rm T}}$ (GeV$^{-2}c^2$)")
    axr.set_xlabel(r"$p_{\rm T}$ (GeV/$c$)"); axr.set_ylabel("MC / data")
    axr.set_ylim(0.5, 1.5)
    stamp(ax, ["pp, $\\sqrt{s}$ = 13 TeV", rf"$|\eta| < {ETA_MAX}$, INEL>0"], loc="lower left")
    ax.legend(loc="upper right")
    save_figure(fig, figdir / "m02_pt_vs_data")
    plt.close(fig)

    out = tabdir / "m02_compare_pt_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"[out]  {rel(figdir / 'm02_pt_vs_data.pdf')}, {rel(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
