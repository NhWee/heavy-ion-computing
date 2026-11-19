#!/usr/bin/env python3
"""Module 5: proton and nuclear PDFs -- what the LHC probes, and what a lead
nucleus changes.

    python analysis/m05_pdfs.py [--figdir results/figures/m05] [--tabdir results/tables/m05]

  1. sum rules of the proton PDF the Monash tune uses (NNPDF2.3 LO):
     int (u - ubar) dx = 2, int (d - dbar) dx = 1, int (s - sbar) dx = 0,
     int x sum_i f_i dx = 1 -- at Q^2 = 10 and 10^4 GeV^2.  They test the grid
     AND the code that reads it, and show that DGLAP evolution conserves them
     while moving momentum from quarks to gluons.
  2. x f(x, Q^2) for u_v, d_v, ubar and g at two scales: at small x the gluon
     dominates and grows with Q^2 -- that is why LHC soft physics is gluon
     physics.
  3. the (x, Q^2) reached at the LHC: x_1,2 = (Q / sqrt(s)) exp(+-y) for
     sqrt(s) = 13 TeV pp and 5.02 TeV Pb-Pb, with three markers -- a
     mid-rapidity hadron (pT ~ 1 GeV, Q ~ 2 pT), the Z, a 100 GeV jet.
  4. EPS09 LO nuclear modification for Pb: R_g, R_valence, R_sea with the
     90 % CL Hessian band, and the values at the markers' x.  This is the
     "cold nuclear matter" baseline below everything hot that follows.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from hicpy.pdfs import ProtonPDF, NuclearRatio, integrate_logx, PROTON_SET   # noqa: E402
from hicpy.style import apply_style, save_figure, OKABE_ITO, plt              # noqa: E402

FLAVOURS = (-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21)
MARKERS = {                           # name: (Q2 in GeV^2, label)
    "hadron": (4.0, r"hadron $p_{\rm T}\approx 1$ GeV"),
    "Z": (91.1876 ** 2, r"$Z$ boson"),
    "jet": (200.0 ** 2, r"jet $p_{\rm T}$ = 100 GeV"),
}
ENERGIES = {"pp 13 TeV": 13000.0, "Pb-Pb 5.02 TeV": 5020.0}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figdir", default=str(REPO / "results" / "figures" / "m05"))
    ap.add_argument("--tabdir", default=str(REPO / "results" / "tables" / "m05"))
    args = ap.parse_args()
    figdir, tabdir = pathlib.Path(args.figdir), pathlib.Path(args.tabdir)
    figdir.mkdir(parents=True, exist_ok=True); tabdir.mkdir(parents=True, exist_ok=True)

    p = ProtonPDF()
    out = {"proton_set": PROTON_SET, "nuclear_set": "EPS09 LO, Pb-208 (PYTHIA built-in)"}

    # ---------------------------------------------------------- 1. sum rules
    rules = {}
    for q2 in (10.0, 1e4):
        rules[f"{q2:g}"] = {
            "u_valence": integrate_logx(lambda x: p.valence(2, x, q2) / x),
            "d_valence": integrate_logx(lambda x: p.valence(1, x, q2) / x),
            "s_valence": integrate_logx(lambda x: p.valence(3, x, q2) / x),
            "momentum": integrate_logx(lambda x: sum(p.xf(i, x, q2) for i in FLAVOURS)),
            "gluon_momentum": integrate_logx(lambda x: p.xf(21, x, q2)),
        }
    out["sum_rules"] = rules
    print("[1] sum rules (expected 2, 1, 0, 1):")
    for q2, r in rules.items():
        print(f"    Q2 = {q2:>5s} GeV2: u_v {r['u_valence']:.4f}  d_v {r['d_valence']:.4f}  "
              f"s_v {r['s_valence']:+.4f}  momentum {r['momentum']:.4f}  (gluon {r['gluon_momentum']:.3f})")

    # ---------------------------------------------------------- 2. xf(x)
    apply_style()
    x = np.geomspace(1e-4, 0.95, 300)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    curves = [(lambda q2: p.valence(2, x, q2), r"$xu_v$", OKABE_ITO[0]),
              (lambda q2: p.valence(1, x, q2), r"$xd_v$", OKABE_ITO[2]),
              (lambda q2: p.xf(-2, x, q2), r"$x\bar u$", OKABE_ITO[3]),
              (lambda q2: p.xf(21, x, q2) / 10, r"$xg/10$", OKABE_ITO[1])]
    for f, lab, col in curves:
        ax.plot(x, f(10.0), "-", color=col, lw=2, label=lab)
        ax.plot(x, f(1e4), "--", color=col, lw=1.6)
    ax.plot([], [], "k-", label=r"$Q^2$ = 10 GeV$^2$"); ax.plot([], [], "k--", label=r"$Q^2$ = 10$^4$ GeV$^2$")
    ax.set_xscale("log"); ax.set_xlim(1e-4, 1); ax.set_ylim(0, 1.2)
    ax.set_xlabel("$x$"); ax.set_ylabel(r"$x f(x, Q^2)$")
    ax.set_title(f"proton, {PROTON_SET.split('_0000')[0]}", fontsize=11)
    ax.legend(fontsize=10, ncol=2, loc="upper right")
    save_figure(fig, figdir / "m05_proton_pdfs"); plt.close(fig)
    xs = {f"{xx:g}": {"xg_Q2_10": float(p.xf(21, xx, 10.)[0]), "xg_Q2_1e4": float(p.xf(21, xx, 1e4)[0]),
                      "xuv_Q2_10": float(p.valence(2, xx, 10.)[0])} for xx in (1e-4, 1e-3, 1e-2, 0.1, 0.3)}
    out["xf_samples"] = xs
    print(f"[2] xg(1e-4): {xs['0.0001']['xg_Q2_10']:.1f} (Q2 = 10) -> {xs['0.0001']['xg_Q2_1e4']:.1f} (Q2 = 1e4); "
          f"xu_v peaks near x ~ 0.2")

    # ---------------------------------------------------------- 3. kinematics
    reach = {}
    for e_name, rs in ENERGIES.items():
        reach[e_name] = {m: {"Q2": q2, "x_y0": float(np.sqrt(q2) / rs),
                             "x_range_absy_lt_0.9": [float(np.sqrt(q2) / rs * np.exp(-0.9)),
                                                     float(np.sqrt(q2) / rs * np.exp(0.9))]}
                         for m, (q2, _) in MARKERS.items()}
    out["kinematics"] = reach
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for k, (e_name, rs) in enumerate(ENERGIES.items()):
        for y, ls in ((0.0, "-"), (2.0, "--"), (4.0, ":")):
            q = np.geomspace(1.0, 2e3, 200)
            ax.plot(q / rs * np.exp(-y), q ** 2, ls, color=OKABE_ITO[k], lw=1.6,
                    label=f"{e_name}, y = {y:g}" if True else None)
        for m, (q2, lab) in MARKERS.items():
            ax.plot(np.sqrt(q2) / rs, q2, "o", color=OKABE_ITO[k], ms=7)
            if k == 1:
                ax.annotate(lab, (np.sqrt(q2) / rs, q2), (8, -4), textcoords="offset points", fontsize=9)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-7, 1); ax.set_ylim(1, 1e7)
    ax.set_xlabel("$x$"); ax.set_ylabel(r"$Q^2$ (GeV$^2$)")
    ax.legend(fontsize=8, ncol=2, loc="upper left")
    save_figure(fig, figdir / "m05_kinematics"); plt.close(fig)
    print("[3] mid-rapidity x: " + "; ".join(
        f"{e}: hadron {r['hadron']['x_y0']:.1e}, Z {r['Z']['x_y0']:.1e}, jet {r['jet']['x_y0']:.1e}"
        for e, r in reach.items()))

    # ---------------------------------------------------------- 4. EPS09 for Pb
    nr = NuclearRatio(p)
    xr = np.geomspace(1e-4, 0.8, 60)       # above ~0.8 Fermi motion makes R blow up
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9), sharey=True)
    eps = {}
    for ax, (name, kind, col) in zip(axes, (("gluon", 21, OKABE_ITO[1]), ("valence", "v", OKABE_ITO[0]),
                                            ("sea (ubar)", -2, OKABE_ITO[3]))):
        for q2, ls in ((10.0, "-"), (1e4, "--")):
            if kind == "v":
                def ratio(m, q2=q2):
                    e = nr._member(m)
                    num = np.array([e.xf(2, xi, q2) - e.xf(-2, xi, q2) + e.xf(1, xi, q2) - e.xf(-1, xi, q2)
                                    for xi in xr])
                    den = (nr.nucleon_xf(2, xr, q2) - nr.nucleon_xf(-2, xr, q2)
                           + nr.nucleon_xf(1, xr, q2) - nr.nucleon_xf(-1, xr, q2))
                    return num / den
                r0 = ratio(0)
                d = 0.5 * np.sqrt(sum((ratio(2 * k - 1) - ratio(2 * k)) ** 2 for k in range(1, 16)))
            else:
                r0, d = nr.ratio_band(kind, xr, q2)
            if ls == "-":
                ax.fill_between(xr, r0 - d, r0 + d, color=col, alpha=0.25, lw=0)
            ax.plot(xr, r0, ls, color=col, lw=2, label=rf"$Q^2$ = {q2:g} GeV$^2$")
            eps[f"{name}_Q2_{q2:g}"] = {"x": xr.tolist(), "R": r0.tolist(), "dR_90cl": d.tolist()}
        ax.axhline(1, color="0.4", lw=1)
        ax.set_xscale("log"); ax.set_xlabel("$x$"); ax.set_title(f"Pb-208, {name}", fontsize=12)
        ax.set_ylim(0.4, 1.6)
    axes[0].set_ylabel(r"$R_i^{\rm Pb}(x, Q^2)$ (EPS09 LO)")
    axes[0].legend(fontsize=9, loc="upper left")
    for xx, lab in ((5e-4, "shadowing"), (0.035, "anti-\nshadowing"), (0.45, "EMC")):
        axes[1].text(xx, 0.47, lab, ha="center", fontsize=9, color="0.35")
    save_figure(fig, figdir / "m05_eps09_pb"); plt.close(fig)
    out["eps09_curves"] = eps

    at = {}
    rs = ENERGIES["Pb-Pb 5.02 TeV"]
    for m, (q2, _) in MARKERS.items():
        xx = np.sqrt(q2) / rs
        rg, dg = nr.ratio_band(21, [xx], max(q2, 1.69))
        at[m] = {"x": float(xx), "Q2": q2, "R_g": float(rg[0]), "dR_g_90cl": float(dg[0])}
    out["eps09_at_PbPb_5.02_midrapidity"] = at
    print("[4] EPS09 R_g for Pb at Pb-Pb 5.02 TeV mid-rapidity: " + "; ".join(
        f"{m} (x = {v['x']:.1e}) {v['R_g']:.2f} ± {v['dR_g_90cl']:.2f}" for m, v in at.items()))

    (tabdir / "m05_pdfs_summary.json").write_text(json.dumps(out, indent=1))
    print(f"[out] {figdir}/m05_*.pdf|png, {tabdir / 'm05_pdfs_summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
