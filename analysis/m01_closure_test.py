#!/usr/bin/env python3
"""Module 1, step 3: CLOSURE TEST -- does the analysis recover the generator input?

    python analysis/m01_closure_test.py [--infile data/generated/m01_toy_events.root]
                                        [--tabdir results/tables]

Why this exists
---------------
"The script ran without crashing" is not evidence that an analysis is correct.
The strongest cheap test available in Monte Carlo work is a *closure test*: we
know the true parameters that went into the generator, so the analysis must give
them back.  The truth is read from the input FILE's provenance (what it was
generated with), never from the current config/analysis.yaml.

What is closed (every generator ingredient of hicpy.toy_generator)
-----------------------------------------------------------------
1. Multiplicity, NBD mean -- <N_ch> in the acceptance against
   nbar' * f_eta * f_pt, where f_eta is the fraction of the generator's dN/deta
   shape inside |eta| < eta_max and f_pt the fraction of the Tsallis spectrum
   inside (pt_min, pt_max), both integrated from the SAME functional forms the
   generator samples; nbar' = nbar + P(N=0) accounts for the max(N, 1) floor.
2. dN_ch/deta averaged over |eta| < 0.5 against nbar' * pdf_eta(|eta| < 0.5).
3. Multiplicity, NBD shape -- k from the moments of the accepted N_ch,
   k = <N>^2 / (Var N - <N>).  Binomial thinning of an NBD keeps k, so the
   accepted distribution must return the generator's k (bootstrap error).
4. pT spectrum -- Tsallis (T, n) per species (pi, K, p), each fitted with its
   own mass by the same binned Poisson likelihood as the main analysis.
   Fitting the *combined* spectrum with the pion mass does NOT close (T comes
   out ~13 MeV, ~20 sigma, too high): that failure is deliberate and is the
   reason experiments publish identified-particle spectra.
5. Z line shape -- m_Z and Gamma_Z from the generated dimuon masses WITHOUT
   any muon cut (truth level), same Poisson Breit-Wigner fit as the analysis.
   (What the muon cuts do to the fitted peak is measured separately, with
   pseudo-experiments, in analysis/m01_fit_validation.py.)

Items 1-3 are compared with the statistical error of the sample estimate;
every pull must be within TOLERANCE_SIGMA.  Results are written to
results/tables/m01_closure_test.json so that the numbers quoted in the docs
are traceable to a file, not to a console log.
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

from hicpy import kinematics as kin                                        # noqa: E402
from hicpy.histograms import make_edges, hist_counts, select_bins          # noqa: E402
from hicpy.fitting import (fit_poisson, tsallis_expected_counts,          # noqa: E402
                           bw_expected_counts, sub_edges)
from hicpy.toy_generator import (ToyConfig, MASS, PDG, M_Z, GAMMA_Z,      # noqa: E402
                                 expected_mean_multiplicity, eta_fraction,
                                 pt_fraction, species_fractions)

SPECIES = [("pi+-", PDG["pi"]), ("K+-", PDG["K"]), ("p/pbar", PDG["p"])]
TOLERANCE_SIGMA = 4.0
N_BOOTSTRAP = 200
TSALLIS_LIMITS = {"C": (0.0, None), "T": (0.01, 1.0), "n": (1.0, 50.0)}


def expected_multiplicities(gen: ToyConfig, eta_max: float, pt_min: float,
                            pt_max: float):
    """Analytic expectation of <N_ch>(|eta|<eta_max, pt_min<pT<pt_max) and of
    the average dN/deta over |eta| < 0.5, from the generator's own shapes
    (hicpy.toy_generator: max(N,1) floor, eta shape, per-species Tsallis)."""
    nbar_eff = expected_mean_multiplicity(gen)
    f_pt = sum(w * pt_fraction(gen, MASS[code], pt_min, pt_max)
               for code, w in species_fractions(gen).items())
    return (nbar_eff * eta_fraction(gen, eta_max) * f_pt,
            nbar_eff * eta_fraction(gen, 0.5) / 1.0)


def nbd_k_moments(n):
    """Method-of-moments NBD shape: k = mean^2 / (var - mean)."""
    m, v = n.mean(), n.var(ddof=1)
    return m * m / (v - m)


def pull_row(label, value, error, truth):
    pull = (value - truth) / error if error > 0 else float("inf")
    return {"observable": label, "value": float(value), "error": float(error),
            "truth": float(truth), "pull": float(pull),
            "ok": bool(abs(pull) <= TOLERANCE_SIGMA)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(REPO / "config" / "analysis.yaml"))
    ap.add_argument("--infile", default=None)
    ap.add_argument("--tabdir", default=None, help="override paths.table_dir")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    acut = cfg["analysis"]["track_cuts"]
    eta_max = acut["eta_abs_max"]
    infile = pathlib.Path(args.infile) if args.infile else REPO / cfg["paths"]["generated_file"]
    tabdir = pathlib.Path(args.tabdir) if args.tabdir else REPO / cfg["paths"]["table_dir"]
    tabdir.mkdir(parents=True, exist_ok=True)

    with uproot.open(infile) as f:
        tree = f["events"]
        n_ev = tree.num_entries
        a = tree.arrays(["particle_px", "particle_py", "particle_pz", "particle_E",
                         "particle_pdg"])
        prov = json.loads(str(f["provenance"]))

    # The truth is what the FILE says it was generated with, not the current yaml.
    gen = ToyConfig(**prov["config"])
    yaml_gen = {k: v for k, v in cfg["generator"].items() if k != "name"}
    if {**yaml_gen, "n_events": gen.n_events} != prov["config"]:
        print("note: the generator block of the config differs from the file's "
              "provenance -- using the file's provenance as truth")

    f64 = lambda arr: ak.values_astype(arr, np.float64)          # noqa: E731
    px, py, pz, e = (f64(a.particle_px), f64(a.particle_py), f64(a.particle_pz),
                     f64(a.particle_E))
    pt = kin.pt(px, py)
    eta = kin.eta(px, py, pz)
    apdg = abs(a.particle_pdg)
    rows, ok = [], True

    print(f"closure test: {infile.name}, {n_ev} events, |eta| < {eta_max} "
          f"(truth from the file's provenance, seed {gen.seed})\n")
    print(f"{'observable':28s} {'measured':>20s} {'truth':>10s} {'pull':>6s}")

    def show(row):
        print(f"{row['observable']:28s} {row['value']:10.4f} +- {row['error']:7.4f} "
              f"{row['truth']:10.4f} {row['pull']:6.1f}")

    # ------------------------------------------------ 1-3) multiplicity closure
    # Charged HADRONS only (pi, K, p): the NBD expectation is for the hadrons
    # the generator draws from the NBD; the two muons of an embedded Z are
    # extra and are deliberately excluded here (the analysis itself counts
    # every charged particle, muons included -- see module01_root.md).
    is_h = (apdg == 211) | (apdg == 321) | (apdg == 2212)
    sel_acc = is_h & (abs(eta) < eta_max) & (pt > acut["pt_min"]) & (pt < acut["pt_max"])
    nch = ak.to_numpy(ak.sum(sel_acc, axis=1)).astype(np.float64)
    n_central = ak.to_numpy(ak.sum(is_h & (pt > acut["pt_min"]) & (abs(eta) < 0.5),
                                   axis=1)).astype(np.float64)
    exp_nch, exp_dndeta = expected_multiplicities(gen, eta_max, acut["pt_min"],
                                                  acut["pt_max"])
    for label, arr, expv in (("<N_ch> in acceptance", nch, exp_nch),
                             ("dN/deta, |eta|<0.5", n_central, exp_dndeta)):
        rows.append(pull_row(label, arr.mean(), arr.std(ddof=1) / np.sqrt(len(arr)), expv))
        show(rows[-1])

    rng = np.random.default_rng(gen.seed)
    boot = [nbd_k_moments(nch[rng.integers(0, len(nch), len(nch))])
            for _ in range(N_BOOTSTRAP)]
    rows.append(pull_row("NBD k (accepted N_ch)", nbd_k_moments(nch),
                         float(np.std(boot, ddof=1)), gen.nbd_k))
    show(rows[-1])

    # ------------------------------------------------ 4) Tsallis closure
    edges_all = make_edges(cfg["analysis"]["binning"]["pt"])
    fr = cfg["analysis"]["fits"]["tsallis_pt_range"]
    in_fr = select_bins(edges_all, fr[0], fr[1])
    edges = sub_edges(edges_all, in_fr)
    norm_deta = 2.0 * eta_max
    print(f"\nTsallis, Poisson likelihood on raw counts, {edges[0]:.3f} < pT < "
          f"{edges[-1]:.3f} GeV/c; generator truth T = {gen.tsallis_T*1000:.1f} MeV, "
          f"n = {gen.tsallis_n:.3f}")
    print(f"{'species':8s} {'N':>9s} {'T (MeV)':>16s} {'pull':>6s} "
          f"{'n':>15s} {'pull':>6s} {'-2lnL/ndf':>11s}")

    def tsallis_fit(sel, m0, label):
        vals = ak.to_numpy(ak.flatten(pt[sel]))
        k, _ = hist_counts(vals, edges_all)
        c0 = k[0] / (n_ev * norm_deta * 2 * np.pi * edges[0] * (edges[1] - edges[0]))
        fit = fit_poisson(tsallis_expected_counts(edges, n_ev, norm_deta, m0), k[in_fr],
                          p0={"C": c0, "T": 0.13, "n": 7.0}, limits=TSALLIS_LIMITS,
                          label=label, function="tsallis_invariant")
        return len(vals), fit

    species_fits = {}
    for name, code in SPECIES:
        n_sel, fit = tsallis_fit((apdg == code) & (abs(eta) < eta_max), MASS[code], name)
        species_fits[name] = fit
        if not fit["converged"]:
            print(f"{name:8s}  FIT DID NOT CONVERGE: {fit.get('error')}")
            ok = False
            continue
        rT = pull_row(f"Tsallis T [{name}]", fit["params"][1], fit["errors"][1], gen.tsallis_T)
        rn = pull_row(f"Tsallis n [{name}]", fit["params"][2], fit["errors"][2], gen.tsallis_n)
        rows += [rT, rn]
        print(f"{name:8s} {n_sel:9d} {rT['value']*1000:8.1f} +- {rT['error']*1000:4.1f} "
              f"{rT['pull']:6.1f} {rn['value']:8.3f} +- {rn['error']:5.3f} {rn['pull']:6.1f} "
              f"{fit['deviance']:6.1f}/{fit['ndf']:<3d}")

    # deliberate non-closure: inclusive charged hadrons fitted with the pion mass
    n_sel, incl = tsallis_fit(is_h & (abs(eta) < eta_max), MASS[PDG["pi"]], "inclusive")
    if incl["converged"]:
        pT_ = (incl["params"][1] - gen.tsallis_T) / incl["errors"][1]
        pn_ = (incl["params"][2] - gen.tsallis_n) / incl["errors"][2]
        print(f"{'incl.':8s} {n_sel:9d} {incl['params'][1]*1000:8.1f} +- "
              f"{incl['errors'][1]*1000:4.1f} {pT_:6.1f} {incl['params'][2]:8.3f} +- "
              f"{incl['errors'][2]:5.3f} {pn_:6.1f} {incl['deviance']:6.1f}/{incl['ndf']:<3d}"
              "   <- expected NOT to close (3-mass mixture fitted with m_pi)")

    # ------------------------------------------------ 5) Z line shape, truth level
    is_mu = apdg == 13
    two = ak.num(px[is_mu]) == 2
    s = lambda arr: ak.to_numpy(ak.sum(arr[is_mu][two], axis=1))      # noqa: E731
    m_true = kin.invariant_mass(s(e), s(px), s(py), s(pz))
    zedges = make_edges(cfg["analysis"]["binning"]["mass"])
    zc, _ = hist_counts(m_true, zedges)
    zfit = fit_poisson(bw_expected_counts(zedges), zc,
                       p0={"N": max(zc.sum(), 1.0), "m0": 91.2, "gamma": 2.5},
                       limits={"N": (0.0, None), "m0": (80.0, 100.0), "gamma": (0.1, 20.0)},
                       label="Z truth (no muon cuts)", function="rel_breit_wigner")
    print(f"\nZ line shape, all {len(m_true)} generated Z -> mu mu (no muon cuts), "
          f"{int(zc.sum())} in {zedges[0]:.0f}-{zedges[-1]:.0f} GeV")
    if zfit["converged"]:
        for lab, i, truth in (("Z mass (truth level)", 1, M_Z),
                              ("Z width (truth level)", 2, GAMMA_Z)):
            rows.append(pull_row(lab, zfit["params"][i], zfit["errors"][i], truth))
            show(rows[-1])
    else:
        print(f"  Z FIT DID NOT CONVERGE: {zfit.get('error')}")
        ok = False

    ok &= all(r["ok"] for r in rows)
    verdict = "PASSED" if ok else "FAILED"
    print(f"\nCLOSURE {verdict}  (all pulls within {TOLERANCE_SIGMA:.0f} sigma: "
          f"multiplicity mean and shape, dN/deta, per-species Tsallis, Z line shape)")

    out = tabdir / "m01_closure_test.json"
    out.write_text(json.dumps({
        "input_file": str(infile.name),
        "n_events": int(n_ev),
        "truth_from": "input file provenance",
        "generator_config": prov["config"],
        "tolerance_sigma": TOLERANCE_SIGMA,
        "passed": bool(ok),
        "checks": rows,
        "fits": {"species": species_fits,
                 "inclusive_pion_mass_expected_not_to_close": incl,
                 "z_truth": zfit},
        "generator_provenance": prov,
    }, indent=2))
    print(f"[out] {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
