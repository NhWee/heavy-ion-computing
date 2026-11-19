#!/usr/bin/env python3
"""Module 4: jets with FastJet -- algorithms, areas, and why heavy ions need
background subtraction.

    python scripts/m02_run_pythia.py --card generators/pythia_hardqcd_13TeV.cmnd \
           --out data/generated/m04_hardqcd_13TeV.root
    python analysis/m04_jets.py [--infile ...] [--n-max N]

Input: Module 2's TTree layout (final-state particles), hard-QCD PYTHIA with
pTHat > 20 GeV.  Jets are built from all visible final-state particles
(neutrinos removed) with |eta| < 3; jets are kept for |y| < 2.

Measurements
  1. inclusive anti-kT jet cross section dsigma/dpT dy for R = 0.2, 0.4, 0.6
     (absolute, from PYTHIA's sigma_gen); jets below 40 GeV are flagged as
     biased by the pTHat cut
  2. jet AREA for anti-kT, kT and Cambridge/Aachen at R = 0.4: anti-kT jets
     should be circles (A ~ pi R^2 = 0.503), kT/CA jets irregular
  3. dijet azimuthal separation of the two leading jets (back-to-back at pi,
     with a tail from extra radiation) and the jet mass
  4. EMBEDDING: every event is overlaid with an uncorrelated toy background
     of known density rho_true; matched jets shift by delta pT, which should
     be ~ rho_true * A on average.  Event-by-event rho from the kT-jet median
     then removes the shift; what survives is the region-to-region
     FLUCTUATION of the background, which no subtraction can remove.

Verification criteria (docs/roadmap.md, Module 4)
  - anti-kT area distribution narrow around pi R^2; kT/CA broad
  - <delta pT> before subtraction = rho_true * <A> within a few %
  - <delta pT> after subtraction consistent with 0
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import awkward as ak
import numpy as np
import uproot

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from hicpy.jets import (pseudojets, area_definition, cluster, median_rho,  # noqa: E402
                        toy_background, delta_phi, delta_r)
from hicpy.style import apply_style, stamp, save_figure, OKABE_ITO, plt                 # noqa: E402

PT_EDGES = np.array([20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200, 300], float)
JET_RAP = 2.0
PART_ETA = 3.0
BIAS_PT = 40.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--infile", default=str(REPO / "data" / "generated" / "m04_hardqcd_13TeV.root"))
    ap.add_argument("--n-max", type=int, default=None)
    ap.add_argument("--bkg-density", type=float, default=30.0, help="toy dN/(deta dphi)")
    ap.add_argument("--bkg-mean-pt", type=float, default=0.5, help="toy <pT> in GeV/c")
    ap.add_argument("--seed", type=int, default=20260827)
    ap.add_argument("--figdir", default=str(REPO / "results" / "figures" / "m04"))
    ap.add_argument("--tabdir", default=str(REPO / "results" / "tables" / "m04"))
    args = ap.parse_args()
    figdir, tabdir = pathlib.Path(args.figdir), pathlib.Path(args.tabdir)
    figdir.mkdir(parents=True, exist_ok=True); tabdir.mkdir(parents=True, exist_ok=True)

    with uproot.open(args.infile) as f:
        t = f["events"]
        ev = t.arrays(["particle_px", "particle_py", "particle_pz", "particle_E", "particle_pdg"],
                      entry_stop=args.n_max)
        prov = json.loads(str(f["provenance"]))
    sigma_mb = prov["cross_section_mb"]["sigma_gen"]
    n_ev = len(ev)
    f64 = lambda a: ak.values_astype(a, np.float64)                      # noqa: E731
    px, py, pz, e = map(f64, (ev.particle_px, ev.particle_py, ev.particle_pz, ev.particle_E))
    pt = np.hypot(px, py); p = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
    eta = np.arctanh((pz / p) * (1 - 1e-12))
    apdg = abs(ev.particle_pdg)
    visible = (apdg != 12) & (apdg != 14) & (apdg != 16)       # neutrinos are invisible
    keep = visible & (abs(eta) < PART_ETA)
    px, py, pz, e = px[keep], py[keep], pz[keep], e[keep]
    print(f"[read] {pathlib.Path(args.infile).name}: {n_ev} events, sigma_gen = {sigma_mb:.4g} mb, "
          f"PYTHIA {prov.get('pythia_version')}")

    rng = np.random.default_rng(args.seed)
    rho_true = args.bkg_density * args.bkg_mean_pt
    ad = area_definition(PART_ETA, 0.01, seed=args.seed)

    spectra = {R: [] for R in (0.2, 0.4, 0.6)}
    areas = {a: [] for a in ("antikt", "kt", "ca")}
    mass04, dphi, lead_pt = [], [], []
    emb = {"pt_sig": [], "pt_emb": [], "area": [], "rho_est": [], "sigma_est": [],
           "area_sig": [], "rho_sig": []}
    for i in range(n_ev):
        parts = pseudojets(px[i], py[i], pz[i], e[i])
        for R in (0.2, 0.6):
            jets, cs = cluster(parts, "antikt", R, PT_EDGES[0], JET_RAP)
            spectra[R].extend(j.pt() for j in jets)
        jets04 = None
        for algo in ("antikt", "kt", "ca"):
            jets, cs = cluster(parts, algo, 0.4, PT_EDGES[0], JET_RAP, ad)
            areas[algo].extend(j.area() for j in jets if j.pt() > BIAS_PT)
            if algo == "antikt":
                spectra[0.4].extend(j.pt() for j in jets)
                mass04.extend(j.m() for j in jets if j.pt() > BIAS_PT)
                if len(jets) >= 2 and jets[0].pt() > BIAS_PT and jets[1].pt() > 20:
                    dphi.append(float(delta_phi(jets[0].phi(), jets[1].phi())))
                if jets:
                    lead_pt.append(jets[0].pt())
                jets04, cs04 = jets, cs
                # ---- embedding: leading signal jet in a toy background ----
                if jets and jets[0].pt() > BIAS_PT:
                    bx, by, bz, be = toy_background(rng, args.bkg_density, args.bkg_mean_pt, PART_ETA)
                    full = parts + pseudojets(bx, by, bz, be)
                    ejets, ecs = cluster(full, "antikt", 0.4, 5.0, JET_RAP + 0.4, ad)
                    match = [j for j in ejets if delta_r(j, jets[0]) < 0.2]
                    if match:
                        m = max(match, key=lambda j: j.pt())
                        rho, sig = median_rho(full, ad, 0.4, JET_RAP, 2)
                        # the PYTHIA event has its own underlying event (MPI); the
                        # median sees it too, so the signal jet is subtracted the
                        # same way before comparing (rho_sig * A_sig)
                        rho_s, _ = median_rho(parts, ad, 0.4, JET_RAP, 2)
                        emb["rho_sig"].append(rho_s); emb["area_sig"].append(jets[0].area())
                        emb["pt_sig"].append(jets[0].pt()); emb["pt_emb"].append(m.pt())
                        emb["area"].append(m.area()); emb["rho_est"].append(rho)
                        emb["sigma_est"].append(sig)
        if (i + 1) % 1000 == 0:
            print(f"  {i + 1} events", flush=True)

    # ---------------------------------------------------------------- numbers
    dy = 2 * JET_RAP
    xs = {}
    for R, v in spectra.items():
        c, _ = np.histogram(v, PT_EDGES)
        w = np.diff(PT_EDGES)
        xs[R] = {"dsigma_dpt_dy_nb": (c * sigma_mb * 1e6 / (n_ev * w * dy)).tolist(),
                 "err_nb": (np.sqrt(c) * sigma_mb * 1e6 / (n_ev * w * dy)).tolist(), "counts": c.tolist()}
    A0 = np.pi * 0.4 ** 2
    area_stats = {a: {"mean": float(np.mean(v)), "rms": float(np.std(v)),
                      "frac_within_10pct_of_piR2": float(np.mean(np.abs(np.asarray(v) / A0 - 1) < 0.1)),
                      "n": len(v)} for a, v in areas.items()}
    E = {k: np.asarray(v) for k, v in emb.items()}
    dpt_raw = E["pt_emb"] - E["pt_sig"]
    dpt_sub = (E["pt_emb"] - E["rho_est"] * E["area"]) - (E["pt_sig"] - E["rho_sig"] * E["area_sig"])
    half68 = lambda x: float(0.5 * (np.percentile(x, 84) - np.percentile(x, 16)))   # noqa: E731
    expect = rho_true * E["area"]
    embedding = {"rho_true": rho_true, "n_jets": int(len(dpt_raw)),
                 "rho_est_mean": float(E["rho_est"].mean()), "rho_est_rms": float(E["rho_est"].std()),
                 "area_mean": float(E["area"].mean()),
                 "rho_ue_pythia_mean": float(E["rho_sig"].mean()),
                 "dpt_raw_mean": float(dpt_raw.mean()), "dpt_raw_rms": float(dpt_raw.std()),
                 "dpt_raw_mean_err": float(dpt_raw.std() / np.sqrt(len(dpt_raw))),
                 "dpt_sub_half68": half68(dpt_sub),
                 "rhoA_expected_mean": float(expect.mean()),
                 "dpt_sub_mean": float(dpt_sub.mean()), "dpt_sub_rms": float(dpt_sub.std()),
                 "dpt_sub_mean_err": float(dpt_sub.std() / np.sqrt(len(dpt_sub))),
                 "sigma_est_mean": float(E["sigma_est"].mean()),
                 "fluct_expected": float(np.sqrt(A0 * args.bkg_density * 2 * args.bkg_mean_pt ** 2))}
    dphi = np.asarray(dphi)
    summary = {"input": pathlib.Path(args.infile).name, "n_events": n_ev, "sigma_gen_mb": sigma_mb,
               "pythia_version": prov.get("pythia_version"), "jet_abs_rap_max": JET_RAP,
               "pt_edges": PT_EDGES.tolist(), "bias_below_pt": BIAS_PT, "spectra": {str(k): v for k, v in xs.items()},
               "areas_R04": area_stats, "jet_mass_R04": {"mean": float(np.mean(mass04)), "n": len(mass04)},
               "dijet_dphi": {"n": int(len(dphi)), "frac_above_2.8": float(np.mean(dphi > 2.8)),
                              "mean": float(dphi.mean())},
               "embedding": embedding, "background": {"dndetadphi": args.bkg_density,
                                                       "mean_pt": args.bkg_mean_pt}}
    (tabdir / "m04_jets_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"[1] anti-kT jets (|y| < {JET_RAP}, pT > {BIAS_PT:.0f}): R=0.2/0.4/0.6 counts "
          + " / ".join(str(int(np.sum(np.asarray(spectra[R]) > BIAS_PT))) for R in (0.2, 0.4, 0.6)))
    print(f"[2] R = 0.4 areas (pi R^2 = {A0:.3f}):")
    for a, s in area_stats.items():
        print(f"      {a:7s} <A> = {s['mean']:.3f}, rms {s['rms']:.3f}, "
              f"{100 * s['frac_within_10pct_of_piR2']:.0f} % within 10 % of pi R^2")
    print(f"[3] dijet dphi: {100 * summary['dijet_dphi']['frac_above_2.8']:.0f} % above 2.8 rad; "
          f"<m_jet>(R=0.4, pT>{BIAS_PT:.0f}) = {summary['jet_mass_R04']['mean']:.1f} GeV")
    print(f"[4] embedding, rho_true = {rho_true:.1f} GeV/area, {len(dpt_raw)} matched jets:")
    print(f"      <rho_est> = {embedding['rho_est_mean']:.2f} ± {embedding['rho_est_rms']:.2f} (event RMS); "
          f"PYTHIA's own UE <rho> = {embedding['rho_ue_pythia_mean']:.2f}")
    print(f"      <dpT> raw = {embedding['dpt_raw_mean']:.2f} ± {embedding['dpt_raw_mean_err']:.2f} GeV vs rho_true <A> = "
          f"{embedding['rhoA_expected_mean']:.2f} GeV")
    print(f"      <dpT> subtracted = {embedding['dpt_sub_mean']:.2f} ± {embedding['dpt_sub_mean_err']:.2f} GeV, "
          f"RMS {embedding['dpt_sub_rms']:.2f}, 68 % half-width {embedding['dpt_sub_half68']:.2f} GeV "
          f"(Poisson-fluctuation estimate {embedding['fluct_expected']:.2f})")

    # ---------------------------------------------------------------- figures
    apply_style()
    ctr = 0.5 * (PT_EDGES[1:] + PT_EDGES[:-1])
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for k, R in enumerate((0.2, 0.4, 0.6)):
        ax.errorbar(ctr * (1 + 0.02 * (k - 1)), xs[R]["dsigma_dpt_dy_nb"], yerr=xs[R]["err_nb"],
                    fmt="os^"[k], ms=5, color=OKABE_ITO[k], label=f"anti-$k_T$, R = {R}")
    ax.axvspan(PT_EDGES[0], BIAS_PT, color="0.92", zorder=0)
    ax.text(29, 2e1, "biased by\n$\\hat{p}_T$ > 20 GeV", fontsize=10, ha="center", color="0.4")
    ax.set_xscale("log"); ax.set_yscale("log")
    from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
    ax.xaxis.set_major_locator(FixedLocator([20, 30, 50, 100, 200, 300]))
    ax.xaxis.set_major_formatter(ScalarFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"jet $p_{\rm T}$ (GeV/$c$)")
    ax.set_ylabel(r"d$^2\sigma$/d$p_{\rm T}$d$y$ (nb c/GeV)")
    stamp(ax, ["PYTHIA 8 Monash, pp 13 TeV", rf"particle level, $|y_{{\rm jet}}| < {JET_RAP:.0f}$"],
          loc="lower left")
    ax.legend(loc="upper right")
    save_figure(fig, figdir / "m04_jet_spectra"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    bins = np.linspace(0.0, 2.0, 81)
    for k, a in enumerate(("antikt", "ca", "kt")):
        lab = {"antikt": "anti-$k_T$", "kt": "$k_T$", "ca": "Cambridge/Aachen"}[a]
        ax.hist(areas[a], bins, histtype="step", lw=2, density=True, color=OKABE_ITO[k], label=lab)
    ax.axvline(A0, color="0.4", ls="--", lw=1.2)
    ax.text(A0 + 0.04, ax.get_ylim()[1] * 0.55, r"$\pi R^2$", color="0.4")
    ax.set_xlabel("jet area A (R = 0.4)"); ax.set_ylabel("probability density")
    ax.legend(loc="upper right")
    save_figure(fig, figdir / "m04_jet_areas"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    b = np.linspace(-15, 25, 81)
    ax.hist(dpt_raw, b, histtype="step", lw=2, color=OKABE_ITO[1], label="embedded − signal")
    ax.hist(dpt_sub, b, histtype="step", lw=2, color=OKABE_ITO[0], label=r"after $\rho A$ subtraction")
    ax.axvline(embedding["rhoA_expected_mean"], color=OKABE_ITO[1], ls="--", lw=1.2)
    ax.axvline(0, color="0.4", lw=1)
    ax.set_xlabel(r"$\delta p_{\rm T}$ (GeV/$c$)"); ax.set_ylabel("jets")
    ax.set_title(rf"toy background $\rho$ = {rho_true:.0f} GeV per unit area; anti-$k_T$ R = 0.4, "
                 rf"$p_{{\rm T}}^{{\rm jet}}$ > {BIAS_PT:.0f} GeV/$c$", fontsize=11)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.25)
    ax.legend(loc="upper left", fontsize=11)
    save_figure(fig, figdir / "m04_embedding"); plt.close(fig)
    print(f"[out] {figdir}/m04_*.pdf|png, {tabdir / 'm04_jets_summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
