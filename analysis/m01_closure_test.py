#!/usr/bin/env python3
"""Module 1, step 3: CLOSURE TEST -- does the analysis recover the generator input?

    python analysis/m01_closure_test.py

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
"""
from __future__ import annotations

import pathlib
import sys

import awkward as ak
import numpy as np
import uproot
import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from m01_analysis_uproot import (tsallis_invariant, do_fit, make_edges,      # noqa: E402
                                 hist_weighted, centers_widths)

SPECIES = [("pi+-", 211, 0.13957039), ("K+-", 321, 0.493677),
           ("p/pbar", 2212, 0.93827208816)]
TOLERANCE_SIGMA = 4.0


def main() -> int:
    cfg = yaml.safe_load(open(REPO / "config" / "analysis.yaml"))
    T_true = cfg["generator"]["tsallis_T"]
    n_true = cfg["generator"]["tsallis_n"]
    eta_max = cfg["analysis"]["track_cuts"]["eta_abs_max"]

    with uproot.open(REPO / cfg["paths"]["generated_file"]) as f:
        tree = f["events"]
        n_ev = tree.num_entries
        a = tree.arrays(["particle_px", "particle_py", "particle_pz", "particle_pdg"])

    pt = np.hypot(a.particle_px, a.particle_py)
    p = np.sqrt(a.particle_px ** 2 + a.particle_py ** 2 + a.particle_pz ** 2)
    eta = np.arctanh((a.particle_pz / p) * (1.0 - 1e-12))
    apdg = abs(a.particle_pdg)

    edges = make_edges(cfg["analysis"]["binning"]["pt"])
    ctr, w = centers_widths(edges)
    norm_deta = 2.0 * eta_max

    print(f"closure test: {n_ev} events, |eta| < {eta_max}")
    print(f"generator truth:  T = {T_true*1000:.1f} MeV,  n = {n_true:.3f}\n")
    print(f"{'species':8s} {'N':>9s} {'T (MeV)':>16s} {'pull':>6s} "
          f"{'n':>15s} {'pull':>6s} {'chi2/ndf':>9s}")

    ok = True
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
          + f"  (per-species pulls within {TOLERANCE_SIGMA:.0f} sigma)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
