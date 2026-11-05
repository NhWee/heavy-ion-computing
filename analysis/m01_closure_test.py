#!/usr/bin/env python3
"""Module 1, step 3: CLOSURE TEST -- does the analysis recover the generator input?

    python analysis/m01_closure_test.py [--infile data/generated/m01_toy_events.root]

Why this exists
---------------
"The script ran without crashing" is not evidence that an analysis is correct.
The strongest cheap test available in Monte Carlo work is a *closure test*: we
know the true parameters that went into the generator, so the analysis must give
them back.  Here we fit the Tsallis form species by species and compare the
extracted (T, n) with the values in config/analysis.yaml.

Expected result
---------------
Per species (pi, K, p) the fit must recover T and n within a few sigma.  Fitting
the *combined* charged-hadron spectrum with a single pion mass does NOT close --
it returns a too-large T -- because the sample is a mixture of three masses.
That failure is deliberate and instructive: it is the same reason experiments
publish identified-particle spectra rather than one inclusive Tsallis fit.

Two more closures use quantities that follow analytically from the generator
configuration, so they test the multiplicity and dN/deta observables too:

* <N_ch> in |eta| < eta_max:  nbar * f_eta * f_pt, where f_eta is the fraction
  of the generator's dN/deta shape inside the acceptance and f_pt the fraction
  of the Tsallis spectrum below pt_max (both computed by numerical integration
  of the SAME functional forms the generator samples from).  The max(N, 1)
  floor in the generator shifts nbar by P(N=0) * 1, which is included.
* dN_ch/deta averaged over |eta| < 0.5:  nbar' * pdf_eta(|eta| < 0.5) / 1.0.

Both are compared with the statistical error of the sample mean.
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import json

import awkward as ak
import numpy as np
import uproot
import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from m01_analysis_uproot import (tsallis_invariant, do_fit, make_edges,      # noqa: E402
                                 hist_weighted, centers_widths)
from hicpy.toy_generator import (ToyConfig, MASS, PDG, _eta_shape,          # noqa: E402
                                 _tsallis_dNdpt)

SPECIES = [("pi+-", 211, 0.13957039), ("K+-", 321, 0.493677),
           ("p/pbar", 2212, 0.93827208816)]
TOLERANCE_SIGMA = 4.0


def expected_multiplicities(gen: ToyConfig, eta_max: float, pt_min: float,
                            pt_max: float):
    """Analytic expectation of <N_ch>(|eta|<eta_max, pt_min<pT<pt_max) and of
    the average dN/deta over |eta| < 0.5, from the generator's own shapes."""
    # effective mean after the max(N, 1) floor: E[max(N,1)] = nbar + P(N=0)
    k, nbar = gen.nbd_k, gen.nbd_nbar
    p0 = (k / (k + nbar)) ** k
    nbar_eff = nbar + p0

    grid = np.linspace(-gen.eta_max, gen.eta_max, 400001)
    pdf = _eta_shape(grid, gen.eta_plateau, gen.eta_sigma)
    pdf /= np.trapezoid(pdf, grid)
    f_eta = np.trapezoid(pdf[np.abs(grid) < eta_max], grid[np.abs(grid) < eta_max])
    f_05 = np.trapezoid(pdf[np.abs(grid) < 0.5], grid[np.abs(grid) < 0.5])

    # pT acceptance, species-weighted (the generator samples per species)
    fr = np.array([gen.frac_pi, gen.frac_K, gen.frac_p]); fr /= fr.sum()
    f_pt = 0.0
    for w, code in zip(fr, (PDG["pi"], PDG["K"], PDG["p"])):
        g = np.geomspace(gen.pt_min, gen.pt_max, 200001)
        d = _tsallis_dNdpt(g, MASS[code], gen.tsallis_T, gen.tsallis_n)
        inside = (g > pt_min) & (g < pt_max)
        f_pt += w * np.trapezoid(d[inside], g[inside]) / np.trapezoid(d, g)

    return nbar_eff * f_eta * f_pt, nbar_eff * f_05 / 1.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(REPO / "config" / "analysis.yaml"))
    ap.add_argument("--infile", default=None)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    T_true = cfg["generator"]["tsallis_T"]
    n_true = cfg["generator"]["tsallis_n"]
    acut = cfg["analysis"]["track_cuts"]
    eta_max = acut["eta_abs_max"]
    gen_cfg = {k: v for k, v in cfg["generator"].items() if k != "name"}

    infile = pathlib.Path(args.infile) if args.infile else REPO / cfg["paths"]["generated_file"]
    with uproot.open(infile) as f:
        tree = f["events"]
        n_ev = tree.num_entries
        a = tree.arrays(["particle_px", "particle_py", "particle_pz", "particle_pdg"])
        prov = json.loads(str(f["provenance"]))
    # the truth is what the FILE says it was generated with, not the current yaml
    gen = ToyConfig(**prov["config"])
    if prov["config"] != {**gen_cfg, "n_events": prov["config"]["n_events"]}:
        print("note: generator config in the file differs from config/analysis.yaml; "
              "using the file's provenance as truth")
        T_true, n_true = gen.tsallis_T, gen.tsallis_n

    pt = np.hypot(a.particle_px, a.particle_py)
    p = np.sqrt(a.particle_px ** 2 + a.particle_py ** 2 + a.particle_pz ** 2)
    eta = np.arctanh((a.particle_pz / p) * (1.0 - 1e-12))
    apdg = abs(a.particle_pdg)

    # ---------------------------------------------------- multiplicity closure
    is_ch = (apdg == 211) | (apdg == 321) | (apdg == 2212)
    sel_acc = is_ch & (abs(eta) < eta_max) & (pt > acut["pt_min"]) & (pt < acut["pt_max"])
    nch = ak.to_numpy(ak.sum(sel_acc, axis=1))
    n_central = ak.to_numpy(ak.sum(is_ch & (pt > acut["pt_min"]) & (abs(eta) < 0.5), axis=1))
    exp_nch, exp_dndeta = expected_multiplicities(gen, eta_max, acut["pt_min"],
                                                  acut["pt_max"])
    ok = True
    print(f"closure test: {n_ev} events, |eta| < {eta_max}\n")
    print(f"{'observable':26s} {'measured':>18s} {'expected':>10s} {'pull':>6s}")
    for label, arr, expv in (("<N_ch> in acceptance", nch, exp_nch),
                             ("dN/deta, |eta|<0.5", n_central, exp_dndeta)):
        mean, err = arr.mean(), arr.std(ddof=1) / np.sqrt(len(arr))
        pull = (mean - expv) / err
        print(f"{label:26s} {mean:10.4f} +- {err:5.4f} {expv:10.4f} {pull:6.1f}")
        if abs(pull) > TOLERANCE_SIGMA:
            ok = False
    print()

    edges = make_edges(cfg["analysis"]["binning"]["pt"])
    ctr, w = centers_widths(edges)
    norm_deta = 2.0 * eta_max

    # ---------------------------------------------------- Tsallis closure
    print(f"generator truth:  T = {T_true*1000:.1f} MeV,  n = {n_true:.3f}\n")
    print(f"{'species':8s} {'N':>9s} {'T (MeV)':>16s} {'pull':>6s} "
          f"{'n':>15s} {'pull':>6s} {'chi2/ndf':>9s}")

    for name, code, m0 in SPECIES:
        sel = (apdg == code) & (abs(eta) < eta_max)
        vals = ak.to_numpy(ak.flatten(pt[sel]))
        wts = 1.0 / (2.0 * np.pi * vals * norm_deta)
        sumw, sumw_err = hist_weighted(vals, wts, edges)
        y, ye = sumw / (n_ev * w), sumw_err / (n_ev * w)
        fit = do_fit(lambda x, C, T, n, _m=m0: tsallis_invariant(x, C, T, n, _m),
                     ctr, y, ye, p0=[y[0] * 2, 0.13, 7.0],
                     bounds=([0, 0.01, 1.0], [np.inf, 1.0, 50.0]), label=name)
        if not fit["converged"]:
            print(f"{name:8s}  FIT DID NOT CONVERGE: {fit.get('error')}")
            ok = False
            continue
        T, dT = fit["params"][1], fit["errors"][1]
        nn, dn = fit["params"][2], fit["errors"][2]
        pull_T, pull_n = (T - T_true) / dT, (nn - n_true) / dn
        print(f"{name:8s} {len(vals):9d} {T*1000:8.1f} +- {dT*1000:4.1f} {pull_T:6.1f} "
              f"{nn:8.3f} +- {dn:5.3f} {pull_n:6.1f} {fit['chi2_ndf']:9.2f}")
        if abs(pull_T) > TOLERANCE_SIGMA or abs(pull_n) > TOLERANCE_SIGMA:
            ok = False

    # deliberate non-closure: inclusive charged hadrons fitted with the pion mass
    sel = ((apdg == 211) | (apdg == 321) | (apdg == 2212)) & (abs(eta) < eta_max)
    vals = ak.to_numpy(ak.flatten(pt[sel]))
    wts = 1.0 / (2.0 * np.pi * vals * norm_deta)
    sumw, sumw_err = hist_weighted(vals, wts, edges)
    y, ye = sumw / (n_ev * w), sumw_err / (n_ev * w)
    fit = do_fit(tsallis_invariant, ctr, y, ye,
                 p0=[y[0] * 2, 0.13, 7.0],
                 bounds=([0, 0.01, 1.0], [np.inf, 1.0, 50.0]), label="inclusive")
    print(f"{'incl.':8s} {len(vals):9d} {fit['params'][1]*1000:8.1f} +- "
          f"{fit['errors'][1]*1000:4.1f} {(fit['params'][1]-T_true)/fit['errors'][1]:6.1f} "
          f"{fit['params'][2]:8.3f} +- {fit['errors'][2]:5.3f} "
          f"{(fit['params'][2]-n_true)/fit['errors'][2]:6.1f} {fit['chi2_ndf']:9.2f}"
          "   <- expected NOT to close (3-mass mixture fitted with m_pi)")

    print("\nCLOSURE " + ("PASSED" if ok else "FAILED")
          + f"  (multiplicity, dN/deta and per-species Tsallis pulls within "
            f"{TOLERANCE_SIGMA:.0f} sigma)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
