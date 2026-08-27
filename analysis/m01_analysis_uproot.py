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
4. m(mu+ mu-)           dimuon invariant mass, with a Breit-Wigner + linear fit

Every figure carries its dataset stamp (system, energy, cuts); every fit records
range, function, initial parameters, convergence status, chi2/ndf and the
covariance matrix into results/tables/m01_fit_results.json.
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
from scipy.optimize import curve_fit

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from hicpy import kinematics as kin              # noqa: E402
from hicpy.style import apply_style, stamp, OKABE_ITO, plt   # noqa: E402

CHARGED_HADRONS = {211, 321, 2212}               # |pdg| of pi+-, K+-, p/pbar


def rel(path: pathlib.Path) -> str:
    """Repo-relative path when possible, absolute otherwise (for logs/manifests)."""
    try:
        return str(pathlib.Path(path).resolve().relative_to(REPO))
    except ValueError:
        return str(path)


# ---------------------------------------------------------------- histograms
def make_edges(spec):
    if spec["type"] == "log":
        return np.geomspace(spec["lo"], spec["hi"], spec["n_bins"] + 1)
    return np.linspace(spec["lo"], spec["hi"], spec["n_bins"] + 1)


def hist_with_errors(values, edges):
    """Unweighted histogram -> (counts, sqrt(N) errors). Poisson errors are the
    correct statistical uncertainty for a raw count in a bin."""
    counts, _ = np.histogram(values, bins=edges)
    return counts.astype(float), np.sqrt(counts.astype(float))


def centers_widths(edges):
    return 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)


# ---------------------------------------------------------------- fit models
def tsallis_invariant(pt, C, T, n, m0=0.13957039):
    """Invariant yield form: C * [1 + (mT - m0)/(nT)]^(-n)."""
    mt = np.sqrt(pt * pt + m0 * m0)
    return C * (1.0 + (mt - m0) / (n * T)) ** (-n)


def bw_plus_linear(m, A, m0, gamma, a, b):
    """Relativistic Breit-Wigner peak on a linear background."""
    num = m ** 2 * gamma ** 2
    den = (m ** 2 - m0 ** 2) ** 2 + m0 ** 2 * gamma ** 2
    return A * num / den + a + b * m


def do_fit(func, x, y, yerr, p0, bounds=(-np.inf, np.inf), label=""):
    """Fit and return a full record: params, errors, chi2/ndf, covariance."""
    mask = (yerr > 0) & np.isfinite(y)
    try:
        popt, pcov = curve_fit(func, x[mask], y[mask], p0=p0,
                               sigma=yerr[mask], absolute_sigma=True,
                               bounds=bounds, maxfev=40000)
        converged = True
    except Exception as exc:                       # noqa: BLE001
        return {"label": label, "converged": False, "error": str(exc)}
    resid = (y[mask] - func(x[mask], *popt)) / yerr[mask]
    chi2 = float(np.sum(resid ** 2))
    ndf = int(mask.sum() - len(popt))
    return {
        "label": label,
        "function": func.__name__,
        "converged": converged,
        "fit_range": [float(x[mask].min()), float(x[mask].max())],
        "n_points": int(mask.sum()),
        "p0": [float(v) for v in p0],
        "params": [float(v) for v in popt],
        "errors": [float(v) for v in np.sqrt(np.diag(pcov))],
        "covariance": [[float(v) for v in row] for row in pcov],
        "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf if ndf > 0 else None,
    }


# ---------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(REPO / "config" / "analysis.yaml"))
    ap.add_argument("--infile", default=None)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    acut = cfg["analysis"]["track_cuts"]
    mcut = cfg["analysis"]["muon_cuts"]
    bins = cfg["analysis"]["binning"]
    sqrt_s = cfg["generator"]["sqrt_s_TeV"]

    infile = pathlib.Path(args.infile) if args.infile else REPO / cfg["paths"]["generated_file"]
    figdir = REPO / cfg["paths"]["figure_dir"]
    tabdir = REPO / cfg["paths"]["table_dir"]
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    apply_style()

    # ---------------------------------------------------------------- read
    with uproot.open(infile) as f:
        tree = f["events"]
        n_events = tree.num_entries
        ev = tree.arrays(["particle_px", "particle_py", "particle_pz",
                          "particle_E", "particle_pdg", "has_z"])
        prov = json.loads(str(f["provenance"]))
    print(f"[read] {infile.name}: {n_events} events, "
          f"{ak.sum(ak.num(ev.particle_px))} particles "
          f"(seed {prov['config']['seed']}, git {prov['git_commit']})")

    px, py, pz, e, pdg = (ev.particle_px, ev.particle_py, ev.particle_pz,
                          ev.particle_E, ev.particle_pdg)

    # NOTE (docs/troubleshooting.md #2): only numpy *ufuncs* (sqrt, hypot,
    # arctanh, comparisons) broadcast over awkward arrays.  np.clip and np.isin
    # are __array_function__ dispatches and try to rectangularise the jagged
    # array, which fails.  Everything below is written with ufuncs only.
    pt_a = np.hypot(px, py)
    p_a = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
    # |pz/|p|| <= 1 by construction; the shrink factor keeps arctanh finite
    eta_a = np.arctanh((pz / p_a) * (1.0 - 1e-12))
    apdg = abs(pdg)
    is_charged_hadron = (apdg == 211) | (apdg == 321) | (apdg == 2212)

    # ------------------------------------------------- 1) multiplicity P(Nch)
    sel_mult = (is_charged_hadron
                & (abs(eta_a) < acut["eta_abs_max"])
                & (pt_a > acut["pt_min"]) & (pt_a < acut["pt_max"]))
    nch = ak.to_numpy(ak.sum(sel_mult, axis=1))

    edges = make_edges(bins["mult"])
    c, err = hist_with_errors(nch, edges)
    ctr, w = centers_widths(edges)
    prob, prob_err = c / (n_events * w), err / (n_events * w)

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
    for extn in ("pdf", "png"):
        fig.savefig(figdir / f"m01_multiplicity.{extn}")
    plt.close(fig)

    # ------------------------------------------------------- 2) dN_ch/deta
    eta_flat = ak.to_numpy(ak.flatten(eta_a[is_charged_hadron
                                            & (pt_a > acut["pt_min"])]))
    edges = make_edges(bins["eta"])
    c, err = hist_with_errors(eta_flat, edges)
    ctr, w = centers_widths(edges)
    dndeta, dndeta_err = c / (n_events * w), err / (n_events * w)

    fig, ax = plt.subplots()
    ax.errorbar(ctr, dndeta, yerr=dndeta_err, fmt="o", color=OKABE_ITO[0],
                label="toy MC (this work)")
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel(r"$\mathrm{d}N_{\rm ch}/\mathrm{d}\eta$")
    ax.set_ylim(0, max(dndeta) * 1.45)
    stamp(ax, [f"toy pp, $\\sqrt{{s}}$ = {sqrt_s:.0f} TeV",
               rf"$p_{{\rm T}} > {acut['pt_min']}$ GeV/$c$, charged hadrons",
               rf"$\mathrm{{d}}N_{{\rm ch}}/\mathrm{{d}}\eta|_{{\eta\approx0}}$ = "
               rf"{dndeta[np.argmin(abs(ctr))]:.2f}"], loc="upper left")
    ax.legend(loc="lower center")
    for extn in ("pdf", "png"):
        fig.savefig(figdir / f"m01_dndeta.{extn}")
    plt.close(fig)

    # ------------------------------------------- 3) invariant pT spectrum
    m = is_charged_hadron & (abs(eta_a) < acut["eta_abs_max"])
    pt_flat = ak.to_numpy(ak.flatten(pt_a[m]))
    edges = make_edges(bins["pt"])
    c, err = hist_with_errors(pt_flat, edges)
    ctr, w = centers_widths(edges)
    delta_eta = 2.0 * acut["eta_abs_max"]
    norm = 2.0 * np.pi * ctr * n_events * w * delta_eta
    yld, yld_err = c / norm, err / norm

    fit = do_fit(tsallis_invariant, ctr, yld, yld_err,
                 p0=[yld[0] * 2, 0.13, 7.0],
                 bounds=([0, 0.01, 1.0], [np.inf, 1.0, 50.0]),
                 label="Tsallis fit to invariant pT yield")

    fig, ax = plt.subplots()
    ax.errorbar(ctr, yld, yerr=yld_err, fmt="o", color=OKABE_ITO[0],
                label="toy MC (this work)")
    if fit["converged"]:
        xs = np.geomspace(ctr[0], ctr[-1], 300)
        ax.plot(xs, tsallis_invariant(xs, *fit["params"]), "-",
                color=OKABE_ITO[1],
                label=(rf"Tsallis: $T$ = {fit['params'][1]*1000:.0f} MeV, "
                       rf"$n$ = {fit['params'][2]:.2f}"))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$p_{\rm T}$ (GeV/$c$)")
    ax.set_ylabel(r"$\frac{1}{2\pi p_{\rm T}N_{\rm ev}}\,"
                  r"\frac{\mathrm{d}^2N}{\mathrm{d}p_{\rm T}\mathrm{d}\eta}$"
                  r"  ($c^2$/GeV$^2$)")
    stamp(ax, [f"toy pp, $\\sqrt{{s}}$ = {sqrt_s:.0f} TeV",
               rf"$|\eta| < {acut['eta_abs_max']}$, charged hadrons",
               rf"$\chi^2/$ndf = {fit.get('chi2_ndf', float('nan')):.2f}"])
    ax.legend(loc="lower left")
    for extn in ("pdf", "png"):
        fig.savefig(figdir / f"m01_pt_spectrum.{extn}")
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
    c, err = hist_with_errors(mass, edges)
    ctr, w = centers_widths(edges)
    zfit = do_fit(bw_plus_linear, ctr, c, np.maximum(err, 1.0),
                  p0=[c.max() * 1e4, 91.2, 2.5, 0.0, 0.0],
                  bounds=([0, 80, 0.1, -1e3, -1e3], [np.inf, 100, 20, 1e3, 1e3]),
                  label="Breit-Wigner + linear fit to m(mu+mu-)")

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.errorbar(ctr, c, yerr=err, fmt="o", color=OKABE_ITO[0],
                label="toy MC (this work)")
    if zfit["converged"]:
        xs = np.linspace(edges[0], edges[-1], 400)
        ax.plot(xs, bw_plus_linear(xs, *zfit["params"]), "-", color=OKABE_ITO[1],
                label="Breit-Wigner + linear fit")
    ax.set_ylim(0, c.max() * 2.2)           # headroom so stamp/legend never overlap data
    ax.axvline(91.1876, ls=":", lw=1.2, color="0.4")
    ax.text(92.5, c.max() * 0.30, r"PDG $m_Z$", fontsize=10, color="0.4")
    ax.set_xlabel(r"$m_{\mu^+\mu^-}$ (GeV/$c^2$)")
    ax.set_ylabel(f"counts / {w[0]:.1f} GeV/$c^2$")
    stamp(ax, [f"toy pp, $\\sqrt{{s}}$ = {sqrt_s:.0f} TeV",
               rf"2 $\mu$, $p_{{\rm T}} > {mcut['pt_min']:.0f}$ GeV/$c$",
               rf"$|\eta_\mu| < {mcut['eta_abs_max']}$",
               f"{len(mass)} candidates",
               (rf"$m_{{\rm fit}}$ = {zfit['params'][1]:.2f} $\pm$ "
                rf"{zfit['errors'][1]:.2f} GeV/$c^2$") if zfit["converged"] else ""],
          loc="upper right")
    ax.legend(loc="upper left", fontsize=10)
    for extn in ("pdf", "png"):
        fig.savefig(figdir / f"m01_dimuon_mass.{extn}")
    plt.close(fig)

    # ------------------------------------------------------------- outputs
    summary = {
        "input_file": rel(infile),
        "n_events": int(n_events),
        "mean_Nch_in_acceptance": float(nch.mean()),
        "rms_Nch": float(nch.std(ddof=1)),
        "dNdeta_at_eta0": float(dndeta[np.argmin(abs(ctr))]) if len(ctr) else None,
        "mean_pt_GeV": float(pt_flat.mean()),
        "n_dimuon_candidates": int(len(mass)),
        "mean_dimuon_mass_GeV": float(mass.mean()) if len(mass) else None,
        "cuts": {"track": acut, "muon": mcut},
        "fits": [fit, zfit],
        "generator_provenance": prov,
    }
    (tabdir / "m01_fit_results.json").write_text(json.dumps(summary, indent=2))

    print(f"[res] <Nch> (|eta|<{acut['eta_abs_max']}) = {nch.mean():.2f} "
          f"+- {nch.std(ddof=1)/np.sqrt(n_events):.2f} (RMS {nch.std(ddof=1):.2f})")
    print(f"[res] <pT> = {pt_flat.mean():.3f} GeV/c")
    if fit["converged"]:
        print(f"[fit] Tsallis  T = {fit['params'][1]*1000:.1f} +- "
              f"{fit['errors'][1]*1000:.1f} MeV, n = {fit['params'][2]:.3f} +- "
              f"{fit['errors'][2]:.3f}, chi2/ndf = {fit['chi2_ndf']:.2f}")
    if zfit["converged"]:
        print(f"[fit] Z peak   m = {zfit['params'][1]:.3f} +- {zfit['errors'][1]:.3f} "
              f"GeV/c2, Gamma = {zfit['params'][2]:.2f} +- {zfit['errors'][2]:.2f} GeV, "
              f"chi2/ndf = {zfit['chi2_ndf']:.2f}")
    print(f"[out] figures -> {rel(figdir)}/m01_*.pdf|png")
    print(f"[out] tables  -> {rel(tabdir / 'm01_fit_results.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
