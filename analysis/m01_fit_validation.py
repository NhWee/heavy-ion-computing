#!/usr/bin/env python3
"""Module 1: VALIDATE THE FIT METHODS with pseudo-experiments.

    python analysis/m01_fit_validation.py [--n-toys-pt 300] [--n-toys-z 1000]
                                          [--tabdir results/tables] [--figdir results/figures]

Why this exists
---------------
The closure test asks "does the analysis give back the truth on THIS sample?"
-- one sample, so a 1-2 sigma deviation cannot be told apart from a
fluctuation.  This script asks the METHOD question: on many independent
samples drawn from the generator's own model, with the statistics of the real
sample, is each fit statistic unbiased (pull mean 0) and are its errors right
(pull width 1)?  pull = (fitted - true) / fitted error.

Studies
-------
A. Tsallis (T, n) per species, fit range and binning from the config
     chi2_centre  least squares at bin centres on the invariant yield
                  (Session 1-3 method; still the three-arm cross-check)
     poisson      Poisson likelihood on raw counts, model integrated over
                  each bin (the physics fit since 2026-11-07)
B. Z line shape (m, Gamma), N_Z ~ Binomial(n_events, z_fraction) per toy
     neyman_cuts      chi2 with sigma = max(sqrt(N), 1), BW + linear
                      (Session 1-3 method), muon cuts applied
     poisson_cuts     Poisson likelihood, BW only (physics fit), muon cuts
     poisson_truth    same, all generated Z, no muon cuts
     poisson_bkg_cuts Poisson likelihood, BW + linear background >= 0, cuts
C. Muon acceptance versus dimuon mass, and the shift of the fitted peak that
   the cuts cause (paired per toy: poisson_cuts - poisson_truth).

Outputs: results/tables/m01_fit_validation.json,
         results/figures/m01_fit_validation.pdf|png
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
import warnings

import numpy as np
import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from hicpy import kinematics as kin                                        # noqa: E402
from hicpy import toy_generator as tg                                      # noqa: E402
from hicpy.histograms import (make_edges, centers_widths, hist_counts,     # noqa: E402
                              hist_weighted, select_bins)
from hicpy.fitting import (tsallis_invariant, bw_plus_linear, fit_chi2,    # noqa: E402
                           fit_poisson, tsallis_expected_counts,
                           bw_expected_counts, bin_integrals, rel_breit_wigner,
                           sub_edges)
from hicpy.style import apply_style, save_figure, OKABE_ITO, plt           # noqa: E402

SPECIES = [("pi+-", tg.PDG["pi"]), ("K+-", tg.PDG["K"]), ("p/pbar", tg.PDG["p"])]
TSALLIS_LIMITS = {"C": (0.0, None), "T": (0.01, 1.0), "n": (1.0, 50.0)}
Z_LIMITS = {"N": (0.0, None), "m0": (80.0, 100.0), "gamma": (0.1, 20.0)}


def pull_stats(values, errors, truth):
    """Bias and pull summary of one parameter over the converged toys."""
    v, e = np.asarray(values, float), np.asarray(errors, float)
    n = len(v)
    if n < 2:
        return {"n_converged": n}
    pull = (v - truth) / e
    return {"n_converged": int(n), "truth": float(truth),
            "bias": float(v.mean() - truth), "bias_err": float(v.std(ddof=1) / np.sqrt(n)),
            "pull_mean": float(pull.mean()), "pull_mean_err": float(pull.std(ddof=1) / np.sqrt(n)),
            "pull_sd": float(pull.std(ddof=1)), "mean_fit_error": float(e.mean()),
            "spread_of_values": float(v.std(ddof=1))}


# --------------------------------------------------------------- study A
def study_tsallis(cfg, gen, rng, n_toys):
    acut = cfg["analysis"]["track_cuts"]
    deta = 2.0 * acut["eta_abs_max"]
    n_ev = gen.n_events
    edges_all = make_edges(cfg["analysis"]["binning"]["pt"])
    ctr, w = centers_widths(edges_all)
    fr = cfg["analysis"]["fits"]["tsallis_pt_range"]
    in_fr = select_bins(edges_all, fr[0], fr[1])
    edges = sub_edges(edges_all, in_fr)
    nbar = tg.expected_mean_multiplicity(gen) * tg.eta_fraction(gen, acut["eta_abs_max"])
    out = {}
    for name, code in SPECIES:
        m0 = tg.MASS[code]
        n_exp = nbar * tg.species_fractions(gen)[code] * n_ev
        res = {"chi2_centre": {"T": [], "dT": [], "n": [], "dn": []},
               "poisson": {"T": [], "dT": [], "n": [], "dn": []}}
        for _ in range(n_toys):
            pts = tg.sample_pt(rng, int(rng.poisson(n_exp)), m0, gen)
            sumw, sumw_err = hist_weighted(pts, 1.0 / (2.0 * np.pi * pts * deta), edges_all)
            y, ye = sumw / (n_ev * w), sumw_err / (n_ev * w)
            f1 = fit_chi2(lambda x, C, T, n, _m=m0: tsallis_invariant(x, C, T, n, _m),
                          ctr[in_fr], y[in_fr], ye[in_fr], p0=[y[0] * 2, 0.13, 7.0],
                          bounds=([0, 0.01, 1.0], [np.inf, 1.0, 50.0]))
            k, _ = hist_counts(pts, edges_all)
            f2 = fit_poisson(tsallis_expected_counts(edges, n_ev, deta, m0), k[in_fr],
                             p0={"C": y[0] * 2, "T": 0.13, "n": 7.0},
                             limits=TSALLIS_LIMITS, n_starts=1)
            for key, f in (("chi2_centre", f1), ("poisson", f2)):
                if f["converged"]:
                    r = res[key]
                    r["T"].append(f["params"][1]); r["dT"].append(f["errors"][1])
                    r["n"].append(f["params"][2]); r["dn"].append(f["errors"][2])
        out[name] = {"n_toys": int(n_toys), "particles_per_toy": float(n_exp)}
        for key, r in res.items():
            out[name][key] = {"T": pull_stats(r["T"], r["dT"], gen.tsallis_T),
                              "n": pull_stats(r["n"], r["dn"], gen.tsallis_n)}
        print(f"  {name:7s} T pull: chi2_centre {out[name]['chi2_centre']['T']['pull_mean']:+.2f}"
              f" / poisson {out[name]['poisson']['T']['pull_mean']:+.2f}   n pull: chi2_centre "
              f"{out[name]['chi2_centre']['n']['pull_mean']:+.2f} / poisson "
              f"{out[name]['poisson']['n']['pull_mean']:+.2f}")
    return out, [float(edges[0]), float(edges[-1])]


# --------------------------------------------------------------- study B, C
def z_masses(rng, n, gen, mcut):
    """True dimuon masses of n generated Z and the mask of those passing the cuts."""
    m_true, mup, mum = tg.sample_z_decays(rng, n, gen)
    ok = np.ones(n, dtype=bool)
    for mu in (mup, mum):
        ok &= (kin.pt(mu[0], mu[1]) > mcut["pt_min"]) & \
              (np.abs(kin.eta(mu[0], mu[1], mu[2])) < mcut["eta_abs_max"])
    s = [a + b for a, b in zip(mup, mum)]
    m = kin.invariant_mass(s[3], s[0], s[1], s[2])          # = m_true up to rounding
    return m, ok


def bw_plus_bkg_counts(edges):
    """N * normalised BW + linear background with non-negative densities
    (counts per GeV) b_lo, b_hi at the window edges."""
    bw = bw_expected_counts(edges)
    ctr, w = centers_widths(edges)
    frac = (ctr - edges[0]) / (edges[-1] - edges[0])

    def expected(N, m0, gamma, b_lo, b_hi):
        return bw(N, m0, gamma) + (b_lo + (b_hi - b_lo) * frac) * w
    return expected


def study_z(cfg, gen, rng, n_toys):
    mcut = cfg["analysis"]["muon_cuts"]
    edges = make_edges(cfg["analysis"]["binning"]["mass"])
    ctr, _ = centers_widths(edges)
    keys = ("neyman_cuts", "poisson_cuts", "poisson_truth", "poisson_bkg_cuts")
    res = {k: {"m": [], "dm": [], "g": [], "dg": []} for k in keys}
    shift = []
    for _ in range(n_toys):
        nz = int(rng.binomial(gen.n_events, gen.z_fraction))
        m, ok = z_masses(rng, nz, gen, mcut)
        c_cut, e_cut = hist_counts(m[ok], edges)
        c_all, _ = hist_counts(m, edges)
        fits = {
            "neyman_cuts": fit_chi2(bw_plus_linear, ctr, c_cut, np.maximum(e_cut, 1.0),
                                    p0=[c_cut.max() * 1e4, 91.2, 2.5, 0.0, 0.0],
                                    bounds=([0, 80, 0.1, -1e3, -1e3],
                                            [np.inf, 100, 20, 1e3, 1e3])),
            "poisson_cuts": fit_poisson(bw_expected_counts(edges), c_cut,
                                        p0={"N": max(c_cut.sum(), 1.0), "m0": 91.2, "gamma": 2.5},
                                        limits=Z_LIMITS, n_starts=1),
            "poisson_truth": fit_poisson(bw_expected_counts(edges), c_all,
                                         p0={"N": max(c_all.sum(), 1.0), "m0": 91.2, "gamma": 2.5},
                                         limits=Z_LIMITS, n_starts=1),
            "poisson_bkg_cuts": fit_poisson(bw_plus_bkg_counts(edges), c_cut,
                                            p0={"N": max(c_cut.sum(), 1.0), "m0": 91.2,
                                                "gamma": 2.5, "b_lo": 0.1, "b_hi": 0.1},
                                            limits={**Z_LIMITS, "b_lo": (0.0, None),
                                                    "b_hi": (0.0, None)}, n_starts=1),
        }
        for k, f in fits.items():
            if f["converged"]:
                i_m, i_g = (1, 2)
                res[k]["m"].append(f["params"][i_m]); res[k]["dm"].append(f["errors"][i_m])
                res[k]["g"].append(f["params"][i_g]); res[k]["dg"].append(f["errors"][i_g])
        if fits["poisson_cuts"]["converged"] and fits["poisson_truth"]["converged"]:
            shift.append(fits["poisson_cuts"]["params"][1] - fits["poisson_truth"]["params"][1])
    out = {"n_toys": int(n_toys)}
    for k, r in res.items():
        out[k] = {"m": pull_stats(r["m"], r["dm"], tg.M_Z),
                  "gamma": pull_stats(r["g"], r["dg"], tg.GAMMA_Z)}
        print(f"  {k:17s} m: bias {out[k]['m']['bias']:+.4f} GeV, pull {out[k]['m']['pull_mean']:+.2f}"
              f" (sd {out[k]['m']['pull_sd']:.2f})   Gamma: bias {out[k]['gamma']['bias']:+.4f} GeV,"
              f" pull {out[k]['gamma']['pull_mean']:+.2f} (sd {out[k]['gamma']['pull_sd']:.2f})")
    sh = np.array(shift)
    out["peak_shift_from_cuts"] = {"mean_GeV": float(sh.mean()),
                                   "error_GeV": float(sh.std(ddof=1) / np.sqrt(len(sh))),
                                   "n_pairs": int(len(sh)),
                                   "definition": "per toy: m(poisson_cuts) - m(poisson_truth), "
                                                 "same generated events"}
    print(f"  peak shift caused by the muon cuts: {sh.mean():+.4f} +- "
          f"{sh.std(ddof=1)/np.sqrt(len(sh)):.4f} GeV")
    return out


def study_acceptance(cfg, gen, rng, n):
    mcut = cfg["analysis"]["muon_cuts"]
    m, ok = z_masses(rng, n, gen, mcut)
    edges = np.arange(80.0, 103.0, 1.0)
    ctr, _ = centers_widths(edges)
    n_all, _ = np.histogram(m, edges)
    n_ok, _ = np.histogram(m[ok], edges)
    acc = n_ok / np.maximum(n_all, 1)
    err = np.sqrt(acc * (1 - acc) / np.maximum(n_all, 1))
    sel = (ctr > 85) & (ctr < 98)
    x = ctr[sel] - tg.M_Z
    W = 1.0 / err[sel] ** 2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    a, b = cov @ (A.T @ (W * acc[sel]))
    window = (m > edges[0]) & (m < edges[-1])
    print(f"  acceptance at m_Z = {a:.4f}, slope = {b*1e3:+.3f} +- {np.sqrt(cov[1,1])*1e3:.3f} "
          f"per mille / GeV ({b/a*100:+.3f} % of itself per GeV)")
    return {"n_decays": int(n), "bin_edges": edges.tolist(), "acceptance": acc.tolist(),
            "acceptance_err": err.tolist(),
            "linear_fit_85_98": {"acc_at_mZ": float(a), "slope_per_GeV": float(b),
                                 "slope_err": float(np.sqrt(cov[1, 1])),
                                 "relative_slope_per_GeV": float(b / a)},
            "mean_mass_80_102": {"all": float(m[window].mean()),
                                 "accepted": float(m[window & ok].mean())}}


# --------------------------------------------------------------- figure
SPECIES_TEX = {"pi+-": r"$\pi^\pm$", "K+-": r"$K^\pm$", "p/pbar": r"$p\,(\bar{p})$"}


def make_figure(res, figdir):
    """Three panels: Tsallis pull means, Z pull means, acceptance vs mass.

    Markers show the mean pull with its error (bias of the method); the
    text "sd" is the pull width (1 = correctly sized errors)."""
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 5.4),
                             gridspec_kw={"width_ratios": [1.1, 1.0, 1.0]})
    style = {"chi2_centre": ("s", OKABE_ITO[1], True, r"$\chi^2$ at bin centres (old)"),
             "poisson": ("o", OKABE_ITO[0], False, "Poisson LL, bin-integrated (new)"),
             "neyman_cuts": ("s", OKABE_ITO[1], True, r"Neyman $\chi^2$, BW+linear, cuts (old)"),
             "poisson_cuts": ("o", OKABE_ITO[0], False, "Poisson LL, BW, cuts (new)"),
             "poisson_truth": ("D", OKABE_ITO[2], False, "Poisson LL, BW, no cuts"),
             "poisson_bkg_cuts": ("v", OKABE_ITO[3], False,
                                  r"Poisson LL, BW + bkg $\geq 0$, cuts")}

    def point(ax, s, y, meth, label):
        mk, col, hollow, lab = style[meth]
        ax.errorbar(s["pull_mean"], y, xerr=s["pull_mean_err"], fmt=mk, color=col,
                    mfc="white" if hollow else col, label=lab if label else None)

    def finish(ax, ticks, xlim, sd_x, title):
        ax.set_yticks([t[0] for t in ticks], [t[1] for t in ticks], fontsize=12)
        ax.tick_params(axis="y", which="minor", left=False, right=False)
        ax.invert_yaxis()
        ax.axvline(0, color="0.5", lw=1)
        ax.set_xlim(*xlim)
        ax.set_xlabel("mean pull  (fit $-$ true) / error")
        ax.set_title(title, fontsize=12, pad=16)
        # x in data units, y in axes units: a column header above the sd values
        ax.text(sd_x, 1.01, "pull sd", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=9, color="0.35")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), fontsize=9.5,
                  frameon=False)

    ax = axes[0]
    ticks, y = [], 0.0
    for par in ("T", "n"):
        for name, _ in SPECIES:
            for meth in ("chi2_centre", "poisson"):
                s = res["tsallis"][name][meth][par]
                point(ax, s, y, meth, par == "T" and name == "pi+-")
                ax.text(1.75, y, f"{s['pull_sd']:.2f}", va="center", ha="center",
                        fontsize=8.5, color="0.35")
                y += 1
            ticks.append((y - 1.5, f"{par}, {SPECIES_TEX[name]}"))
        y += 0.6
    finish(ax, ticks, (-0.6, 2.0), 1.75,
           "Tsallis fit, %d toys per species" % res["tsallis"]["pi+-"]["n_toys"])

    ax = axes[1]
    ticks, y = [], 0.0
    for par, plab in (("m", r"$m_Z$"), ("gamma", r"$\Gamma_Z$")):
        for meth in ("neyman_cuts", "poisson_cuts", "poisson_truth", "poisson_bkg_cuts"):
            s = res["z"][meth][par]
            point(ax, s, y, meth, par == "m")
            ax.text(0.7, y, f"{s['pull_sd']:.2f}", va="center", ha="center",
                    fontsize=8.5, color="0.35")
            y += 1
        ticks.append((y - 2.5, plab))
        y += 0.6
    finish(ax, ticks, (-1.1, 0.9), 0.7, "Z line shape, %d toys" % res["z"]["n_toys"])

    ax = axes[2]
    a = res["acceptance"]
    e = np.array(a["bin_edges"])
    c = 0.5 * (e[1:] + e[:-1])
    lo, hi = min(a["acceptance"]), max(a["acceptance"])
    xs = np.linspace(80.5, 101.5, 300)
    bw = rel_breit_wigner(xs, tg.M_Z, tg.GAMMA_Z)
    ax.fill_between(xs, lo - 0.012, lo - 0.012 + (hi - lo) * 0.55 * bw / bw.max(),
                    color="0.85", lw=0, label="Z line shape (arb. units)")
    ax.errorbar(c, a["acceptance"], yerr=a["acceptance_err"], fmt="o", color=OKABE_ITO[0],
                label="toy Z: both muons accepted")
    lf = a["linear_fit_85_98"]
    xl = np.linspace(85, 98, 50)
    ax.plot(xl, lf["acc_at_mZ"] + lf["slope_per_GeV"] * (xl - tg.M_Z), "-",
            color=OKABE_ITO[1],
            label=rf"linear fit 85-98 GeV: {lf['relative_slope_per_GeV']*100:+.2f} %/GeV")
    ax.set_xlabel(r"$m_{\mu^+\mu^-}$ (GeV/$c^2$)")
    ax.set_ylabel("muon-pair acceptance")
    ax.set_title(r"$p_{\rm T}^\mu > 20$ GeV/$c$, $|\eta_\mu| < 2.4$", fontsize=12, pad=16)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), fontsize=9.5, frameon=False)
    # explicit margins: tight_layout would shrink the axes to make room for the
    # legends placed below them
    fig.subplots_adjust(left=0.08, right=0.985, top=0.90, bottom=0.36, wspace=0.42)
    save_figure(fig, figdir / "m01_fit_validation")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(REPO / "config" / "analysis.yaml"))
    ap.add_argument("--n-toys-pt", type=int, default=None)
    ap.add_argument("--n-toys-z", type=int, default=None)
    ap.add_argument("--n-z-acceptance", type=int, default=None)
    ap.add_argument("--tabdir", default=None)
    ap.add_argument("--figdir", default=None)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    val = dict(cfg["validation"])
    for k in ("n_toys_pt", "n_toys_z", "n_z_acceptance"):
        if getattr(args, k) is not None:
            val[k] = getattr(args, k)
    gen = tg.ToyConfig(**{k: v for k, v in cfg["generator"].items() if k != "name"})
    tabdir = pathlib.Path(args.tabdir) if args.tabdir else REPO / cfg["paths"]["table_dir"]
    figdir = pathlib.Path(args.figdir) if args.figdir else REPO / cfg["paths"]["figure_dir"]
    tabdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    warnings.filterwarnings("ignore", category=RuntimeWarning)
    # independent, reproducible streams for the three studies
    rng_a, rng_b, rng_c = (np.random.default_rng(s)
                           for s in np.random.SeedSequence(val["seed"]).spawn(3))
    t0 = time.time()
    print(f"A. Tsallis methods, {val['n_toys_pt']} toys per species")
    tsallis, fit_range = study_tsallis(cfg, gen, rng_a, val["n_toys_pt"])
    print(f"B. Z line-shape methods, {val['n_toys_z']} toys")
    z = study_z(cfg, gen, rng_b, val["n_toys_z"])
    print(f"C. muon acceptance vs mass, {val['n_z_acceptance']} decays")
    acc = study_acceptance(cfg, gen, rng_c, val["n_z_acceptance"])
    res = {"validation_config": val, "generator_config": gen.to_dict(),
           "tsallis_fit_range_GeV": fit_range, "tsallis": tsallis, "z": z,
           "acceptance": acc,
           "pull_definition": "(fitted - true) / fitted error; truth = generator inputs"}
    (tabdir / "m01_fit_validation.json").write_text(json.dumps(res, indent=2))
    make_figure(res, figdir)
    print(f"[out] {tabdir / 'm01_fit_validation.json'}  ({time.time() - t0:.0f} s)")
    print(f"[out] {figdir / 'm01_fit_validation'}.pdf|png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
