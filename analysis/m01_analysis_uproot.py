#!/usr/bin/env python3
"""Module 1, step 2: read the ROOT file back and produce the first physics figures.

    python analysis/m01_analysis_uproot.py [--config config/analysis.yaml]

This is the "modern Python" arm of Module 1: uproot + awkward + numpy +
matplotlib, no ROOT installation needed.  `analysis/m01_analysis_root.C` and
`analysis/m01_analysis_pyroot.py` do the same thing with ROOT itself, so the two
ecosystems can be compared directly.

Observables produced
--------------------
1. P(N_ch)              charged multiplicity distribution in |eta| < 0.8
2. dN_ch/deta           charged-particle pseudorapidity density
3. invariant pT yield   1/(2 pi pT N_ev) d2N/(dpT deta), with a Tsallis fit
4. m(mu+ mu-)           dimuon invariant mass, with a Breit-Wigner fit

Fits (hicpy.fitting; docs/troubleshooting.md #15, #16)
------------------------------------------------------
* Physics results use a binned Poisson likelihood on the RAW counts with the
  model integrated over each bin.  The Tsallis fit therefore has no
  bin-centre approximation, and the 240-event Z peak is fitted with the right
  statistics (a chi2 with sigma = sqrt(N) biased its width low by 0.8 sigma).
* The Tsallis fit is ALSO done the Session 1-3 way -- least squares at the bin
  centres -- because that algorithm is identical in the PyROOT and C++ arms
  (TF1 chi2).  It is written to the "tsallis" block that
  m01_compare_arms.py compares; it is a cross-check of the implementations,
  not a physics number.

Every figure carries its dataset stamp (system, energy, cuts); every fit records
range, function, method, initial parameters, convergence status, goodness of
fit, covariance and correlation into results/tables/m01_fit_results.json.
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

from hicpy import kinematics as kin              # noqa: E402
from hicpy.events import read_particles, track_selection   # noqa: E402
from hicpy import pdg as pdgref                  # noqa: E402  (reference values)
from hicpy.style import apply_style, stamp, save_figure, OKABE_ITO, plt   # noqa: E402
# Shared helpers (moved to src/hicpy; re-exported here for older imports such
# as analysis/m02_compare_data.py)
from hicpy.histograms import (make_edges, centers_widths, hist_counts,  # noqa: E402
                              hist_weighted, hist_record, select_bins)
from hicpy.fitting import (tsallis_invariant, rel_breit_wigner, fit_chi2,  # noqa: E402
                           fit_poisson, tsallis_expected_counts,
                           bw_expected_counts, bin_integrals, sub_edges)

# "Charged particle" = electric charge != 0, read from the particle_charge
# branch.  NOT a list of PDG codes: a list written for the toy sample (pi, K, p)
# silently undercounts a PYTHIA sample (e, mu, Sigma, Xi, Omega ...) and no
# longer matches the experimental "charged primary" definition.


def rel(path: pathlib.Path) -> str:
    """Repo-relative path when possible, absolute otherwise (for logs/manifests)."""
    try:
        return str(pathlib.Path(path).resolve().relative_to(REPO))
    except ValueError:
        return str(path)


# ---------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(REPO / "config" / "analysis.yaml"))
    ap.add_argument("--infile", default=None)
    ap.add_argument("--figdir", default=None, help="override paths.figure_dir")
    ap.add_argument("--tabdir", default=None, help="override paths.table_dir")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    acut = cfg["analysis"]["track_cuts"]
    mcut = cfg["analysis"]["muon_cuts"]
    bins = cfg["analysis"]["binning"]

    infile = pathlib.Path(args.infile) if args.infile else REPO / cfg["paths"]["generated_file"]
    figdir = pathlib.Path(args.figdir) if args.figdir else REPO / cfg["paths"]["figure_dir"]
    tabdir = pathlib.Path(args.tabdir) if args.tabdir else REPO / cfg["paths"]["table_dir"]
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    apply_style()

    # ---------------------------------------------------------------- read
    # hicpy.events.read_particles promotes the float32 branches to float64
    # BEFORE any arithmetic (docs/troubleshooting.md #12: numpy's weighted
    # histogram accumulates a cumulative sum, and in float32 the sparse high-pT
    # bins came out quantised, their sum of w^2 exactly zero) and computes pT
    # and eta = asinh(pz/pT) with ufunc-only hicpy.kinematics, which works on
    # jagged awkward arrays (troubleshooting #3).
    ev = read_particles(infile)
    n_events, prov = ev.n_events, ev.provenance
    print(f"[read] {infile.name}: {n_events} events, "
          f"{ak.sum(ak.num(ev.px))} particles "
          f"(seed {prov['config']['seed']}, git {prov['git_commit']})")
    # The figure stamp describes the FILE that was analysed, so the collision
    # energy comes from its provenance: the toy records sqrt_s_TeV, PYTHIA
    # (scripts/m02_run_pythia.py) records eCM_GeV.  The config value is only a
    # fallback for files that record neither.
    pcfg = prov.get("config", {})
    if "sqrt_s_TeV" in pcfg:
        sqrt_s = float(pcfg["sqrt_s_TeV"])
    elif "eCM_GeV" in pcfg:
        sqrt_s = float(pcfg["eCM_GeV"]) / 1000.0
    else:
        sqrt_s = float(cfg["generator"]["sqrt_s_TeV"])

    px, py, pz, e = ev.px, ev.py, ev.pz, ev.E
    pt_a, eta_a = ev.pt, ev.eta          # eta = asinh(pz/pT), same in all three arms
    apdg = abs(ev.pdg)
    is_charged = ev.charge != 0

    # ------------------------------------------------- 1) multiplicity P(Nch)
    sel_mult = track_selection(ev, acut)      # charged, |eta|, pt_min < pT < pt_max
    nch = ak.to_numpy(ak.sum(sel_mult, axis=1))

    edges = make_edges(bins["mult"])
    c, err = hist_counts(nch, edges)
    ctr, w = centers_widths(edges)
    prob, prob_err = c / (n_events * w), err / (n_events * w)
    h_mult = hist_record(edges, prob, prob_err)

    fig, ax = plt.subplots()
    ax.errorbar(ctr, prob, yerr=prob_err, fmt="o", color=OKABE_ITO[0],
                label="toy MC (this work)")
    ax.set_yscale("log")
    ax.set_xlabel(r"$N_{\rm ch}$")
    ax.set_ylabel(r"$P(N_{\rm ch})$")
    stamp(ax, [f"toy pp, $\\sqrt{{s}}$ = {sqrt_s:.0f} TeV",
               rf"$|\eta| < {acut['eta_abs_max']}$",
               rf"$p_{{\rm T}} > {acut['pt_min']}$ GeV/$c$",
               rf"$\langle N_{{\rm ch}}\rangle$ = {nch.mean():.1f}"])
    ax.legend(loc="lower left")
    save_figure(fig, figdir / "m01_multiplicity")
    plt.close(fig)

    # ------------------------------------------------------- 2) dN_ch/deta
    eta_flat = ak.to_numpy(ak.flatten(eta_a[is_charged
                                            & (pt_a > acut["pt_min"])]))
    edges = make_edges(bins["eta"])
    c, err = hist_counts(eta_flat, edges)
    ctr, w = centers_widths(edges)
    dndeta, dndeta_err = c / (n_events * w), err / (n_events * w)
    h_eta = hist_record(edges, dndeta, dndeta_err)

    # A single "value at eta = 0" is ambiguous: with edges at 0 the two arms can
    # pick different bins.  Define it instead as the average density over
    # |eta| < 0.5, computed from the counts directly (no binning at all).
    n_central = int(np.sum(np.abs(eta_flat) < 0.5))
    dndeta_0 = n_central / (n_events * 1.0)

    fig, ax = plt.subplots()
    ax.errorbar(ctr, dndeta, yerr=dndeta_err, fmt="o", color=OKABE_ITO[0],
                label="toy MC (this work)")
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel(r"$\mathrm{d}N_{\rm ch}/\mathrm{d}\eta$")
    ax.set_ylim(0, max(dndeta) * 1.45)
    stamp(ax, [f"toy pp, $\\sqrt{{s}}$ = {sqrt_s:.0f} TeV",
               rf"$p_{{\rm T}} > {acut['pt_min']}$ GeV/$c$, charged particles",
               rf"$\mathrm{{d}}N_{{\rm ch}}/\mathrm{{d}}\eta|_{{\eta\approx0}}$ = "
               rf"{dndeta_0:.2f}"], loc="upper left")
    ax.legend(loc="lower center")
    save_figure(fig, figdir / "m01_dndeta")
    plt.close(fig)

    # ------------------------------------------- 3) invariant pT spectrum
    # Same track selection as N_ch (charged, |eta| cut, pt_min < pT < pt_max)
    # so that the three arms histogram exactly the same particles.
    pt_flat = ak.to_numpy(ak.flatten(pt_a[sel_mult]))
    edges = make_edges(bins["pt"])
    ctr, w = centers_widths(edges)
    delta_eta = 2.0 * acut["eta_abs_max"]
    # NOTE (docs/troubleshooting.md #7): the 1/(2 pi pT) Jacobian is applied
    # PER PARTICLE, with that particle's own pT -- not with the bin centre.
    # On a log pT axis the bins are wide and the two differ by several percent,
    # which is enough to make this arm and the ROOT arm disagree.  Per-particle
    # weighting is also what experiments do.
    weights = 1.0 / (2.0 * np.pi * pt_flat * delta_eta)
    sumw, sumw_err = hist_weighted(pt_flat, weights, edges)
    yld, yld_err = sumw / (n_events * w), sumw_err / (n_events * w)
    h_pt = hist_record(edges, yld, yld_err)

    # Fit range is a stated choice (config), not "wherever the bins happen to
    # be empty".  Tsallis describes the soft hadron spectrum; in this toy the
    # bins above a few GeV/c are dominated by Z-decay muons (z_fraction = 2 %
    # is ~1e5 times the real rate) and in Module 2 by the hard-scattering tail.
    fr = cfg["analysis"]["fits"]["tsallis_pt_range"]
    in_range = select_bins(edges, fr[0], fr[1])     # bin centre inside, as ROOT's "R"

    # (a) physics fit: Poisson likelihood on the RAW counts of the same
    #     particles, model = N_ev * delta_eta * integral over the bin of
    #     2 pi pT f_inv.  Same C, T, n meaning as the invariant-yield form.
    counts_pt, _ = hist_counts(pt_flat, edges)
    fit = fit_poisson(tsallis_expected_counts(sub_edges(edges, in_range), n_events,
                                              delta_eta),
                      counts_pt[in_range],
                      p0={"C": yld[0] * 2, "T": 0.13, "n": 7.0},
                      limits={"C": (0.0, None), "T": (0.01, 1.0), "n": (1.0, 50.0)},
                      label="Tsallis fit to the pT spectrum (Poisson likelihood)",
                      function="tsallis_invariant")
    fit["requested_range"] = [float(fr[0]), float(fr[1])]
    fit["fit_range"] = [float(sub_edges(edges, in_range)[0]),
                        float(sub_edges(edges, in_range)[-1])]

    # (b) arm cross-check: least squares at the bin centres, the algorithm the
    #     ROOT arms run (TF1 chi2).  Not a physics result -- see the docstring.
    fit_xc = fit_chi2(tsallis_invariant, ctr[in_range], yld[in_range], yld_err[in_range],
                      p0=[yld[0] * 2, 0.13, 7.0],
                      bounds=([0, 0.01, 1.0], [np.inf, 1.0, 50.0]),
                      label="Tsallis chi2 at bin centres (three-arm cross-check only)")
    fit_xc["requested_range"] = [float(fr[0]), float(fr[1])]

    fig, ax = plt.subplots()
    ax.errorbar(ctr, yld, yerr=yld_err, fmt="o", color=OKABE_ITO[0],
                label="toy MC (this work)")
    if fit["converged"]:
        xs = np.geomspace(*fit["fit_range"], 300)
        ax.plot(xs, tsallis_invariant(xs, *fit["params"]), "-",
                color=OKABE_ITO[1],
                label=("Tsallis fit (Poisson LL)\n"
                       rf"$T$ = {fit['params'][1]*1000:.0f} MeV, "
                       rf"$n$ = {fit['params'][2]:.2f}"))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$p_{\rm T}$ (GeV/$c$)")
    ax.set_ylabel(r"$\frac{1}{2\pi p_{\rm T}N_{\rm ev}}\,"
                  r"\frac{\mathrm{d}^2N}{\mathrm{d}p_{\rm T}\mathrm{d}\eta}$"
                  r"  ($c^2$/GeV$^2$)")
    stamp(ax, [f"toy pp, $\\sqrt{{s}}$ = {sqrt_s:.0f} TeV",
               rf"$|\eta| < {acut['eta_abs_max']}$, charged particles",
               (rf"$-2\ln\lambda$/ndf = {fit['deviance']:.1f}/{fit['ndf']}"
                if fit["converged"] else "fit failed")])
    ax.legend(loc="lower left", fontsize=10)
    save_figure(fig, figdir / "m01_pt_spectrum")
    plt.close(fig)

    # -------------------------------------------------- 4) dimuon mass
    is_mu = apdg == 13
    mu_ok = is_mu & (pt_a > mcut["pt_min"]) & (abs(eta_a) < mcut["eta_abs_max"])
    n_mu = ak.num(px[mu_ok])
    keep = n_mu == 2
    def sum_over_pair(arr):
        return ak.to_numpy(ak.sum(arr[mu_ok][keep], axis=1))

    mass = kin.invariant_mass(sum_over_pair(e), sum_over_pair(px),
                              sum_over_pair(py), sum_over_pair(pz))

    edges = make_edges(bins["mass"])
    c, err = hist_counts(mass, edges)
    ctr, w = centers_widths(edges)
    h_mass = hist_record(edges, c, err)        # raw counts, not normalised
    # Poisson likelihood with a relativistic Breit-Wigner and NO background
    # term: the toy contains no dimuon background by construction, and a
    # background term bounded at >= 0 absorbs the Breit-Wigner tails and
    # biases the width low (-0.3 sigma; m01_fit_validation.py).  Module 2
    # (Drell-Yan continuum, heavy-flavour muons) must add one.
    zfit = fit_poisson(bw_expected_counts(edges), c,
                       p0={"N": max(c.sum(), 1.0), "m0": 91.2, "gamma": 2.5},
                       limits={"N": (0.0, None), "m0": (80.0, 100.0), "gamma": (0.1, 20.0)},
                       label="Breit-Wigner fit to m(mu+mu-) (Poisson likelihood)",
                       function="rel_breit_wigner")

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.errorbar(ctr, c, yerr=err, fmt="o", color=OKABE_ITO[0],
                label="toy MC (this work)")
    if zfit["converged"]:
        n_fit, m_fit, g_fit = zfit["params"]
        xs = np.linspace(edges[0], edges[-1], 400)
        norm = bin_integrals(rel_breit_wigner, edges, m_fit, g_fit).sum()
        # expected counts per bin of width w[0], drawn as a smooth curve
        ax.plot(xs, n_fit * rel_breit_wigner(xs, m_fit, g_fit) / norm * w[0], "-",
                color=OKABE_ITO[1], label="Breit-Wigner fit (Poisson LL)")
    ax.set_ylim(0, c.max() * 2.2)           # headroom so stamp/legend never overlap data
    ax.axvline(pdgref.Z_MASS, ls=":", lw=1.2, color="0.4", label=r"RPP 2026 $m_Z$")
    ax.set_xlabel(r"$m_{\mu^+\mu^-}$ (GeV/$c^2$)")
    ax.set_ylabel(f"counts / {w[0]:.1f} GeV/$c^2$")
    stamp(ax, [f"toy pp, $\\sqrt{{s}}$ = {sqrt_s:.0f} TeV",
               rf"2 $\mu$, $p_{{\rm T}} > {mcut['pt_min']:.0f}$ GeV/$c$",
               rf"$|\eta_\mu| < {mcut['eta_abs_max']}$",
               f"{len(mass)} candidates",
               (rf"$m_{{\rm fit}}$ = {zfit['params'][1]:.2f} $\pm$ "
                rf"{zfit['errors'][1]:.2f} GeV/$c^2$") if zfit["converged"] else "",
               (rf"$\Gamma_{{\rm fit}}$ = {zfit['params'][2]:.2f} $\pm$ "
                rf"{zfit['errors'][2]:.2f} GeV") if zfit["converged"] else ""],
          loc="upper right")
    ax.legend(loc="upper left", fontsize=10)
    save_figure(fig, figdir / "m01_dimuon_mass")
    plt.close(fig)

    # ------------------------------------------------------------- outputs
    summary = {
        "arm": "uproot",
        "input_file": rel(infile),
        "n_events": int(n_events),
        "mean_Nch_in_acceptance": float(nch.mean()),
        "rms_Nch": float(nch.std(ddof=1)),
        "dNdeta_central": dndeta_0,      # averaged over |eta| < 0.5
        "mean_pt_GeV": float(pt_flat.mean()),
        "n_dimuon_candidates": int(len(mass)),
        "mean_dimuon_mass_GeV": float(mass.mean()) if len(mass) else None,
        "cuts": {"track": acut, "muon": mcut},
        "fits": [fit, zfit, fit_xc],
        # Physics result (Poisson likelihood) in flat form
        "tsallis_poisson": ({"T": fit["params"][1], "T_err": fit["errors"][1],
                             "n": fit["params"][2], "n_err": fit["errors"][2],
                             "deviance": fit["deviance"], "ndf": fit["ndf"]}
                            if fit["converged"] else None),
        "z_fit": ({"m": zfit["params"][1], "m_err": zfit["errors"][1],
                   "gamma": zfit["params"][2], "gamma_err": zfit["errors"][2],
                   "n_signal": zfit["params"][0], "deviance": zfit["deviance"],
                   "ndf": zfit["ndf"]} if zfit["converged"] else None),
        # Cross-check block: the bin-centre chi2 fit, in the flat form every arm
        # writes, so m01_compare_arms.py can compare the three implementations.
        "tsallis": ({"T": fit_xc["params"][1], "T_err": fit_xc["errors"][1],
                     "n": fit_xc["params"][2], "n_err": fit_xc["errors"][2],
                     "chi2": fit_xc["chi2"], "ndf": fit_xc["ndf"],
                     "method": "chi2 at bin centres (cross-check only)"}
                    if fit_xc["converged"] else None),
        "histograms": {"mult": h_mult, "eta": h_eta, "pt": h_pt, "mass": h_mass},
        "generator_provenance": prov,
    }
    (tabdir / "m01_fit_results.json").write_text(json.dumps(summary, indent=2))

    print(f"[res] <Nch> (|eta|<{acut['eta_abs_max']}) = {nch.mean():.2f} "
          f"+- {nch.std(ddof=1)/np.sqrt(n_events):.2f} (RMS {nch.std(ddof=1):.2f})")
    print(f"[res] <pT> = {pt_flat.mean():.3f} GeV/c")
    if fit["converged"]:
        print(f"[fit] Tsallis  T = {fit['params'][1]*1000:.1f} +- "
              f"{fit['errors'][1]*1000:.1f} MeV, n = {fit['params'][2]:.3f} +- "
              f"{fit['errors'][2]:.3f}, -2lnL/ndf = {fit['deviance']:.1f}/{fit['ndf']} "
              f"(Poisson likelihood)")
    if fit_xc["converged"]:
        print(f"[x-chk] Tsallis chi2 at bin centres: T = {fit_xc['params'][1]*1000:.1f} "
              f"MeV, n = {fit_xc['params'][2]:.3f}, chi2/ndf = {fit_xc['chi2_ndf']:.2f} "
              f"(for m01_compare_arms.py only)")
    if zfit["converged"]:
        print(f"[fit] Z peak   m = {zfit['params'][1]:.3f} +- {zfit['errors'][1]:.3f} "
              f"GeV/c2, Gamma = {zfit['params'][2]:.2f} +- {zfit['errors'][2]:.2f} GeV, "
              f"-2lnL/ndf = {zfit['deviance']:.1f}/{zfit['ndf']} (Poisson likelihood)")
    print(f"[out] figures -> {rel(figdir)}/m01_*.pdf|png")
    print(f"[out] tables  -> {rel(tabdir / 'm01_fit_results.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
