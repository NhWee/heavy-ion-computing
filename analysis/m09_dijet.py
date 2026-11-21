#!/usr/bin/env python3
"""Module 9b: dijet momentum imbalance A_J from the same energy-loss geometry.

    python analysis/m09_dijet.py pp [--n-events 200000] [--workers N]   # PYTHIA pp 5.02 TeV, anti-kT R = 0.4
    python analysis/m09_dijet.py model                                  # quench, select, A_J; needs m09 model + glauber

QUESTION
--------
CMS (arXiv:1102.1957, PbPb 2.76 TeV) found that in central collisions many
back-to-back dijets are strongly unbalanced: A_J = (pT1 - pT2)/(pT1 + pT2)
for leading jets pT1 > 120 GeV/c, subleading pT2 > 50 GeV/c, Delta phi >
2 pi/3, and the fraction of balanced dijets R_B (A_J < 0.15, the median of
the PYTHIA distribution, so R_B = 50 % in pp by construction) falls rapidly
with centrality.  In the Module 9 picture the two jets start at the same
point and cross DIFFERENT path lengths.  Does the energy loss calibrated on
single hadrons (analysis/m09_quenching.py) produce that imbalance?

MODEL
-----
- pp: PYTHIA 8 hard QCD, pT-hat > 80 GeV (unweighted), 5.02 TeV; anti-kT
  R = 0.4 on all visible final particles (|eta| < 3), jets with |eta| < 2.
- Each event is placed at a random collision point of a Glauber Pb-Pb event
  of the centrality class (the pool of analysis/m09_quenching.py glauber).
- Each jet loses   Delta E_jet = f * kappa * I_zeta(x0, phi_jet) * (pT/10 GeV)^beta
  with kappa, zeta, beta of the four hadron-calibrated models (m09_model.json).
  f = the fraction of the parton's lost energy that leaves the R = 0.4 cone
  (the rest is radiated inside and stays in the jet); f = 1 and 0.5.
- After the loss the jets are re-ordered and the CMS selection applied.
- Loss fluctuations ("fluct" models): Delta E_jet multiplied by an
  exponentially distributed factor of mean 1 per jet (kappa recalibrated on
  hadrons with the same fluctuations, analysis/m09_quenching.py).
- Background fluctuations: in Pb-Pb a jet's pT after rho*A subtraction still
  carries delta pT ~ Gaussian(0, sigma_bkg) from the underlying event.  For
  independently emitted particles (Poisson number, exponential pT -- exactly
  the case Module 4 validated by embedding: 2.73 vs 2.75 GeV)
      sigma_bkg^2 = A * dN/(deta dphi) * 2 <pT>^2,   A = pi R^2,
  with dN/(deta dphi) = 1.5 x (dNch/deta) / 2 pi (charged + neutral) and
  dNch/deta, <pT> of the Angantyr Pb-Pb centrality class (Module 7,
  results/tables/m07/m07_PbPb.json).  Correlated fluctuations (flow,
  mini-jets) make the real width larger: a lower bound.  Each jet of an
  event gets an independent delta pT, also in the "no loss" baseline, so
  the background's own effect on A_J is visible.
NOT in the model: jet-by-jet f, fake jets made of background, the medium
response that is recovered at large angles.

OUTPUT
  data/generated/m09_dijet_pp.npz (not in git); results/tables/m09/m09_dijet.json;
  results/figures/m09/m09_dijet_aj.png, m09_dijet_rb.png
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from hicpy.glauber import path_integrals                       # noqa: E402

GEN, TAB, FIG = REPO / "data" / "generated", REPO / "results" / "tables" / "m09", REPO / "results" / "figures" / "m09"
CARD = REPO / "generators" / "pythia_hardqcd_5TeV.cmnd"
N_JETS = 4
SEL = {"pt1": 120.0, "pt2": 50.0, "dphi": 2 * np.pi / 3, "eta": 2.0}
GROUPS = {"0-10%": (0, 1), "10-30%": (2, 3), "30-50%": (4, 5), "50-80%": (6, 7, 8)}
AJ_EDGES = np.linspace(0, 1, 21)
F_OUT = (1.0, 0.5)
WIDTH = 0.5
R_JET = 0.4
CHARGED_TO_ALL = 1.5


def background_sigma():
    """sigma_bkg (GeV) per Module 9 centrality class from the Angantyr Pb-Pb
    classes of Module 7 (same edges for 0-80 %)."""
    m7 = json.loads((REPO / "results" / "tables" / "m07" / "m07_PbPb.json").read_text())
    out = []
    for c in m7["classes"][:9]:
        dens = CHARGED_TO_ALL * c["dndeta"] / (2 * np.pi)
        out.append(float(np.sqrt(np.pi * R_JET ** 2 * dens * 2 * c["meanpt"]["pi"] ** 2)))
    return np.array(out)


def _chunk(job):
    n_events, seed = job
    try:
        import pythia8
    except ImportError:
        import pythia8mc as pythia8
    from hicpy.jets import cluster, pseudojets, NEUTRINOS
    py = pythia8.Pythia("", False)
    py.readFile(str(CARD))
    for s in (f"Random:seed = {seed}", "PhaseSpace:pTHatMin = 80."):
        py.readString(s)
    if not py.init():
        raise RuntimeError("PYTHIA init() failed")
    out = np.zeros((n_events, N_JETS, 3))
    n = 0
    while n < n_events:
        if not py.next():
            continue
        ev = py.event
        sel = [i for i in range(ev.size()) if ev[i].isFinal() and ev[i].idAbs() not in NEUTRINOS and abs(ev[i].eta()) < 3.0]
        parts = pseudojets([ev[i].px() for i in sel], [ev[i].py() for i in sel],
                           [ev[i].pz() for i in sel], [ev[i].e() for i in sel])
        jets, cs = cluster(parts, "antikt", 0.4, 30.0, SEL["eta"])       # keep cs alive while jets are read
        for k, j in enumerate(jets[:N_JETS]):
            out[n, k] = j.pt(), j.eta(), j.phi_std()
        n += 1
    info = py.infoPython() if hasattr(py, "infoPython") else py.info
    return out, info.sigmaGen(), py.settings.parm("Pythia:versionNumber")


def make_pp(args) -> int:
    w = args.workers or os.cpu_count() or 1
    jobs = [(args.n_events // w + (i < args.n_events % w), 20261301 + i) for i in range(w)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=w) as ex:
        parts = list(ex.map(_chunk, jobs))
    jets = np.concatenate([p[0] for p in parts])
    GEN.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(GEN / "m09_dijet_pp.npz", jets=jets, sigma_mb=np.mean([p[1] for p in parts]), version=parts[0][2])
    lead = jets[:, 0, 0]
    print(f"[pp] {len(jets)} events (pT-hat > 80 GeV), leading jet > 120 GeV/c in {np.mean(lead > 120):.1%}; "
          f"PYTHIA {parts[0][2]:.3f}, {time.time() - t0:.0f} s")
    return 0


def select_aj(j):
    """j: (n, N_JETS, 3) jets after loss.  Re-order by pT, apply the CMS dijet
    selection; returns A_J of the selected events and their mask."""
    order = np.argsort(-j[:, :, 0], axis=1)
    js = np.take_along_axis(j, order[:, :, None], axis=1)
    pt1, pt2 = js[:, 0, 0], js[:, 1, 0]
    dphi = np.abs((js[:, 0, 2] - js[:, 1, 2] + np.pi) % (2 * np.pi) - np.pi)
    ok = (pt1 > SEL["pt1"]) & (pt2 > SEL["pt2"]) & (dphi > SEL["dphi"])
    return (pt1[ok] - pt2[ok]) / (pt1[ok] + pt2[ok]), ok


def model(args) -> int:
    from hicpy.style import apply_style, save_figure, OKABE_ITO, plt
    jets = np.load(GEN / "m09_dijet_pp.npz")["jets"]
    gz = np.load(GEN / "m09_glauber.npz", allow_pickle=True)
    cls, coll, wounded = gz["cls"], gz["coll"], gz["wounded"]
    mods = json.loads((TAB / "m09_model.json").read_text())["models"]
    mods = {"no loss": {"kappa": 0.0, "zeta": 0, "beta": 0.0, "fluct": False}, **mods}
    sig_bkg = background_sigma()
    rng = np.random.default_rng(12)
    real = jets[:, :, 0] > 0
    aj_pp, ok_pp = select_aj(jets)
    res = {"selection": {k: float(v) for k, v in SEL.items()}, "n_pp_events": int(len(jets)),
           "pp": {"n_dijets": int(ok_pp.sum()), "mean_aj": float(aj_pp.mean()), "rb_015": float(np.mean(aj_pp < 0.15)),
                  "median_aj": float(np.median(aj_pp)), "hist": np.histogram(aj_pp, AJ_EDGES, density=True)[0].tolist()},
           "aj_edges": AJ_EDGES.tolist(), "groups": {}}
    print(f"[pp] {ok_pp.sum()} selected dijets, <A_J> {aj_pp.mean():.3f}, median {np.median(aj_pp):.3f}, "
          f"R_B(A_J<0.15) {np.mean(aj_pp < 0.15):.2f}")
    res["sigma_bkg_GeV"] = sig_bkg.tolist()
    print("[bkg] sigma(delta pT), R = 0.4, per class 0-5 ... 70-80 %: " + " ".join(f"{v:.1f}" for v in sig_bkg))
    x_jet = rng.exponential(1.0, size=jets.shape[:2])                  # loss fluctuation per jet
    gauss = rng.normal(size=jets.shape[:2])                             # background fluctuation per jet
    for gname, members in GROUPS.items():
        pool = np.nonzero(np.isin(cls, members))[0]
        g_of = rng.choice(pool, size=len(jets))
        sb = sig_bkg[cls[g_of]][:, None] * gauss                        # delta pT of each jet
        I = np.zeros((len(jets), N_JETS, 2))
        for g in np.unique(g_of):
            ev = np.nonzero(g_of == g)[0]
            c, wd = coll[g].astype(float), wounded[g].astype(float)
            x0 = np.repeat(c[rng.integers(0, len(c), size=len(ev))], N_JETS, axis=0)
            phi = jets[ev, :, 2].ravel()
            for z in (0, 1):
                I[ev, :, z] = path_integrals(x0, phi, wd, WIDTH, z).reshape(len(ev), N_JETS)
        res["groups"][gname] = {}
        for mname, m in mods.items():
            for f in (F_OUT if mname != "no loss" else (1.0,)):
              for bkg in (False, True):
                pt = jets[:, :, 0]
                de = f * m["kappa"] * I[:, :, m["zeta"]] * (np.clip(pt, 1e-3, None) / 10.0) ** m["beta"]
                if m.get("fluct"):
                    de = de * x_jet
                ptq = np.clip(pt - de, 0, None) + (sb if bkg else 0.0)
                q = jets.copy(); q[:, :, 0] = np.where(real, np.clip(ptq, 0, None), 0.0)
                aj, ok = select_aj(q)
                # surface bias: path integral of the (post-loss) leading and subleading jet
                order = np.argsort(-q[:, :, 0], axis=1)
                Is = np.take_along_axis(I[:, :, m["zeta"]], order, axis=1)
                key = f"{mname}, f={f:g}" + (", bkg" if bkg else "")
                res["groups"][gname][key] = {
                    "n_dijets": int(ok.sum()), "mean_aj": float(aj.mean()) if len(aj) else None,
                    "rb_015": float(np.mean(aj < 0.15)) if len(aj) else None,
                    "hist": np.histogram(aj, AJ_EDGES, density=True)[0].tolist() if len(aj) else None,
                    "survival": float(ok.sum() / max(ok_pp.sum(), 1)),
                    "mean_I_leading": float(Is[ok, 0].mean()) if ok.any() else None,
                    "mean_I_subleading": float(Is[ok, 1].mean()) if ok.any() else None}
        print(f"[{gname}] " + "  ".join(f"{k.replace(', f=1', '')}: R_B {v['rb_015']:.2f}"
                                      for k, v in res["groups"][gname].items() if ", f=1" in k))
    TAB.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    (TAB / "m09_dijet.json").write_text(json.dumps(res, indent=1))

    apply_style()
    c = 0.5 * (AJ_EDGES[1:] + AJ_EDGES[:-1])
    base = [n for n in mods if n != "no loss" and not n.endswith("fluct")]
    col = dict(zip(base, OKABE_ITO[:4]))
    colr = lambda n: "0.45" if n == "no loss" else col[n.replace(", fluct", "")]        # noqa: E731
    style = lambda n: "--" if n.endswith("fluct") else "-"                              # noqa: E731
    g0 = res["groups"]["0-10%"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), sharey=True)
    for ax, bkg in zip(axes, (False, True)):
        sfx = ", bkg" if bkg else ""
        ax.step(c, res["pp"]["hist"], where="mid", color="k", lw=2, label="pp (PYTHIA)")
        for n in mods:
            if n == "no loss" and not bkg:
                continue
            v = g0[f"{n}, f=1{sfx}"]
            ax.step(c, v["hist"], where="mid", color=colr(n), ls=style(n), lw=1.4, label=n)
        ax.set_xlabel(r"$A_J = (p_{T,1} - p_{T,2}) / (p_{T,1} + p_{T,2})$"); ax.set_xlim(0, 0.8)
        ax.set_title("0-10 %, f = 1, " + ("with background fluctuations" if bkg else "no background"), fontsize=11)
    axes[0].set_ylabel("probability density"); axes[1].legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout(); save_figure(fig, FIG / "m09_dijet_aj"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    xs = np.arange(len(GROUPS))
    for n in mods:
        ax.plot(xs, [res["groups"][g][f"{n}, f=1, bkg"]["rb_015"] for g in GROUPS], style(n), marker="o", ms=4,
                color=colr(n), label=n + (" (background only)" if n == "no loss" else ""))
    ax.axhline(res["pp"]["rb_015"], color="k", lw=1.5, label="pp")
    ax.set_xticks(xs); ax.set_xticklabels(list(GROUPS)); ax.set_xlabel("centrality")
    ax.set_ylabel(r"$R_B$ = fraction with $A_J < 0.15$ (f = 1, with background)")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout(); save_figure(fig, FIG / "m09_dijet_rb"); plt.close(fig)
    print(f"[out] {TAB / 'm09_dijet.json'}, {FIG}/m09_dijet_*.png")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("pp"); a.add_argument("--n-events", type=int, default=200000); a.add_argument("--workers", type=int)
    sub.add_parser("model")
    args = ap.parse_args()
    return {"pp": make_pp, "model": model}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
