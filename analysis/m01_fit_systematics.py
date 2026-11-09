#!/usr/bin/env python3
"""Module 1: SYSTEMATIC uncertainty of the Tsallis fit from analysis choices.

    python analysis/m01_fit_systematics.py [--infile ...] [--tabdir ...] [--figdir ...]

The statistical error of a fit says how much (T, n) would scatter if the
experiment were repeated with the SAME analysis.  It says nothing about how
much they depend on choices the analyst made: where the fit starts and ends,
how the spectrum is binned.  This script varies those choices one at a time
(config: systematics.tsallis) and keeps statistical and systematic
uncertainties separate, as the final numbers must.

Procedure
---------
For each variation v of the nominal fit, delta_v = x_v - x_nom.  The two fits
share most of their data, so delta_v is NOT compared with their statistical
errors but with its own statistical spread sigma_delta,v: how much delta_v
moves when the experiment is repeated.  That spread is measured by bootstrap
-- resample the particles with replacement, redo the nominal and every
variation on each replica -- which is valid for any kind of variation.
A variation counts as a systematic effect only if
|delta_v| > threshold * sigma_delta,v; within each group (fit-range low edge,
fit-range high edge, binning) the largest significant |delta| is taken, and
the groups are added in quadrature.

The shortcut sigma_B = sqrt(|sigma_v^2 - sigma_nom^2|) (Barlow,
hep-ex/0207026) is also reported.  It holds when one data set is a subset of
the other (fit-range changes) but not for re-binning the same data, where it
goes to ~0 and turns statistical noise into a "significant" shift -- the
pion control below shows exactly that (troubleshooting #18).

Two spectra are analysed
------------------------
* inclusive charged particles with the pion-mass Tsallis form -- the physics
  result of m01_analysis_uproot.py.  The form is NOT the true model of this
  three-species mixture, so its parameters are expected to depend on the fit
  range: that dependence IS the systematic uncertainty.
* pions alone with the pion mass -- the TRUE generator model.  A control:
  here no variation should be significant.  If one is, the procedure (not the
  physics) is producing false systematics.

Output: results/tables/m01_tsallis_systematics.json and
        results/figures/m01_tsallis_systematics.pdf|png
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

from hicpy.events import read_particles, track_selection                   # noqa: E402
from hicpy.histograms import make_edges, hist_counts, select_bins          # noqa: E402
from hicpy.fitting import fit_poisson, tsallis_expected_counts, sub_edges  # noqa: E402
from hicpy.toy_generator import M_PI                                       # noqa: E402
from hicpy.style import apply_style, save_figure, OKABE_ITO, plt           # noqa: E402

LIMITS = {"C": (0.0, None), "T": (0.01, 1.0), "n": (1.0, 50.0)}


def fit_spectrum(pts, n_ev, deta, binning, fit_range, n_starts=3):
    edges_all = make_edges(binning)
    sel = select_bins(edges_all, *fit_range)
    edges = sub_edges(edges_all, sel)
    k, _ = hist_counts(pts, edges_all)
    c0 = max(k[sel][0], 1.0) / (n_ev * deta * 2 * np.pi * edges[0] * (edges[1] - edges[0]))
    f = fit_poisson(tsallis_expected_counts(edges, n_ev, deta, M_PI), k[sel],
                    p0={"C": c0, "T": 0.13, "n": 7.0}, limits=LIMITS, n_starts=n_starts)
    f["fit_range"] = [float(edges[0]), float(edges[-1])]
    f["n_bins_total"] = int(binning["n_bins"])
    return f


def variations(cfg):
    nom_range = list(cfg["analysis"]["fits"]["tsallis_pt_range"])
    nom_bins = dict(cfg["analysis"]["binning"]["pt"])
    sy = cfg["systematics"]["tsallis"]
    out = [("nominal", "nominal", nom_range, nom_bins)]
    for lo in sy["range_low"]:
        out.append(("range_low", f"fit from {lo} GeV/c", [lo, nom_range[1]], nom_bins))
    for hi in sy["range_high"]:
        out.append(("range_high", f"fit up to {hi} GeV/c", [nom_range[0], hi], nom_bins))
    for nb in sy["n_bins"]:
        out.append(("binning", f"{nb} log bins", nom_range, {**nom_bins, "n_bins": nb}))
    return out


def analyse(pts, n_ev, deta, cfg, threshold, rng_boot, n_boot):
    rows, nom = [], None
    var = variations(cfg)
    for group, label, rng, binning in var:
        f = fit_spectrum(pts, n_ev, deta, binning, rng)
        row = {"group": group, "label": label, "fit_range_requested": rng,
               "fit_range": f.get("fit_range"), "n_bins_total": binning["n_bins"],
               "converged": f["converged"]}
        if f["converged"]:
            row.update(T=f["params"][1], T_err=f["errors"][1],
                       n=f["params"][2], n_err=f["errors"][2],
                       deviance=f["deviance"], ndf=f["ndf"])
        if group == "nominal":
            nom = row
        rows.append(row)
    # bootstrap: the statistical spread of every delta_v, correlations included
    boot = {par: [[] for _ in var] for par in ("T", "n")}
    for _ in range(n_boot):
        rep_pts = pts[rng_boot.integers(0, len(pts), len(pts))]
        fits = [fit_spectrum(rep_pts, n_ev, deta, b, r, n_starts=1) for _, _, r, b in var]
        if not fits[0]["converged"]:
            continue
        for i, f in enumerate(fits[1:], start=1):
            if f["converged"]:
                boot["T"][i].append(f["params"][1] - fits[0]["params"][1])
                boot["n"][i].append(f["params"][2] - fits[0]["params"][2])
    summary = {}
    for par in ("T", "n"):
        groups = {}
        for i, r in enumerate(rows[1:], start=1):
            if not r["converged"]:
                continue
            d = r[par] - nom[par]
            sb = np.sqrt(abs(r[f"{par}_err"] ** 2 - nom[f"{par}_err"] ** 2))
            sd = float(np.std(boot[par][i], ddof=1)) if len(boot[par][i]) > 2 else float("nan")
            r[f"{par}_delta"] = float(d)
            r[f"{par}_sigma_delta_bootstrap"] = sd
            r[f"{par}_n_bootstrap"] = len(boot[par][i])
            r[f"{par}_barlow_sigma"] = float(sb)
            r[f"{par}_significant"] = bool(np.isfinite(sd) and abs(d) > threshold * sd)
            r[f"{par}_significant_by_barlow"] = bool(abs(d) > threshold * sb)
            if r[f"{par}_significant"]:
                groups[r["group"]] = max(groups.get(r["group"], 0.0), abs(d))
        syst = float(np.sqrt(sum(v * v for v in groups.values())))
        summary[par] = {"nominal": nom[par], "stat": nom[f"{par}_err"], "syst": syst,
                        "syst_by_group": groups,
                        "combination": "max significant |delta| per group, groups in quadrature"}
    return rows, summary


def make_figure(results, figdir):
    """delta_v = x_v - x_nominal with its bootstrap statistical spread (the
    right error bar for a shift between overlapping fits), the nominal
    statistical error as a grey band for scale, stars on significant shifts."""
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))
    series = ((("inclusive", r"inclusive, $m_\pi$ form (physics result)"), OKABE_ITO[1], "s", -0.12),
              (("pion_control", r"$\pi^\pm$ only, true model (control)"), OKABE_ITO[0], "o", 0.12))
    for ax, par, unit, scale in ((axes[0], "T", "MeV", 1e3), (axes[1], "n", "", 1.0)):
        stat = results["inclusive"]["summary"][par]["stat"] * scale
        ax.axhspan(-stat, stat, color="0.9", lw=0, label="nominal stat. error (inclusive)")
        for (key, lab), col, mk, off in series:
            rows = results[key]["variations"][1:]
            x = np.arange(1, len(rows) + 1) + off
            d = np.array([r.get(f"{par}_delta", np.nan) for r in rows]) * scale
            e = np.array([r.get(f"{par}_sigma_delta_bootstrap", np.nan) for r in rows]) * scale
            ax.errorbar(x, d, yerr=e, fmt=mk, color=col, mfc=col if mk == "o" else "white",
                        label=lab)
            sig = np.array([bool(r.get(f"{par}_significant")) for r in rows])
            ax.plot(x[sig], d[sig], "*", color="k", ms=10, zorder=5,
                    label="significant (> 2 bootstrap $\\sigma$)" if key == "inclusive" else None)
        ax.axhline(0, color="0.5", lw=1)
        labels = [r["label"] for r in results["inclusive"]["variations"][1:]]
        ax.set_xticks(range(1, len(labels) + 1), labels, rotation=35, ha="right", fontsize=10)
        ax.tick_params(axis="x", which="minor", bottom=False, top=False)
        ax.set_xlim(0.4, len(labels) + 0.6)
        ax.set_ylabel(rf"${par} - {par}_{{\rm nominal}}$" + (f" ({unit})" if unit else ""))
    axes[0].legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    save_figure(fig, figdir / "m01_tsallis_systematics")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(REPO / "config" / "analysis.yaml"))
    ap.add_argument("--infile", default=None)
    ap.add_argument("--tabdir", default=None)
    ap.add_argument("--figdir", default=None)
    ap.add_argument("--n-bootstrap", type=int, default=None,
                    help="override systematics.tsallis.n_bootstrap (tests)")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    infile = pathlib.Path(args.infile) if args.infile else REPO / cfg["paths"]["generated_file"]
    tabdir = pathlib.Path(args.tabdir) if args.tabdir else REPO / cfg["paths"]["table_dir"]
    figdir = pathlib.Path(args.figdir) if args.figdir else REPO / cfg["paths"]["figure_dir"]
    tabdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    sy = cfg["systematics"]["tsallis"]
    threshold = float(sy["significance_threshold"])
    n_boot = int(args.n_bootstrap if args.n_bootstrap is not None else sy["n_bootstrap"])

    ev = read_particles(infile)
    acut = cfg["analysis"]["track_cuts"]
    deta = 2.0 * acut["eta_abs_max"]
    incl = ak.to_numpy(ak.flatten(ev.pt[track_selection(ev, acut)]))
    pions = ak.to_numpy(ak.flatten(ev.pt[(abs(ev.pdg) == 211) & (abs(ev.eta) < acut["eta_abs_max"])]))

    results = {}
    seeds = np.random.SeedSequence(int(sy["seed"])).spawn(2)
    for (key, pts), ss in zip((("inclusive", incl), ("pion_control", pions)), seeds):
        rows, summary = analyse(pts, ev.n_events, deta, cfg, threshold,
                                np.random.default_rng(ss), n_boot)
        results[key] = {"n_particles": int(len(pts)), "variations": rows, "summary": summary}
        print(f"\n{key}: {len(pts)} particles")
        print(f"  {'variation':22s} {'T (MeV)':>16s} {'dT':>7s} {'sig':>5s} {'n':>14s} {'dn':>7s} {'sig':>5s}"
              f"   (sig = delta / bootstrap spread; Barlow shortcut in brackets)")
        for r in rows:
            if not r["converged"]:
                print(f"  {r['label']:22s} FIT DID NOT CONVERGE")
                continue
            dT = r.get("T_delta", 0.0) * 1e3
            dn = r.get("n_delta", 0.0)
            sig = lambda p, k: (r[f"{p}_delta"] / r[k]) if r.get(k) else 0.0   # noqa: E731
            sT, sn = sig("T", "T_sigma_delta_bootstrap"), sig("n", "n_sigma_delta_bootstrap")
            bT, bn = sig("T", "T_barlow_sigma"), sig("n", "n_barlow_sigma")
            print(f"  {r['label']:22s} {r['T']*1e3:8.2f} +- {r['T_err']*1e3:4.2f} {dT:+7.2f} {sT:+5.1f} "
                  f"{r['n']:7.3f} +- {r['n_err']:5.3f} {dn:+7.3f} {sn:+5.1f}"
                  + (f"   [{bT:+5.1f} {bn:+5.1f}]" if r["group"] != "nominal" else ""))
        s = summary
        print(f"  => T = {s['T']['nominal']*1e3:.1f} +- {s['T']['stat']*1e3:.1f} (stat) "
              f"+- {s['T']['syst']*1e3:.1f} (syst) MeV,  n = {s['n']['nominal']:.3f} +- "
              f"{s['n']['stat']:.3f} (stat) +- {s['n']['syst']:.3f} (syst)")

    out = {"input_file": infile.name, "significance_threshold": threshold,
           "n_bootstrap": n_boot,
           "procedure": "significant if |delta| > threshold * bootstrap sd of delta "
                        "(particles resampled, nominal and variation refitted on each "
                        "replica); max significant |delta| per group, groups in quadrature. "
                        "Barlow's sqrt(|sigma_v^2 - sigma_nom^2|) reported for comparison.",
           **results, "generator_provenance": ev.provenance}
    (tabdir / "m01_tsallis_systematics.json").write_text(json.dumps(out, indent=2))
    make_figure(results, figdir)
    print(f"\n[out] {tabdir / 'm01_tsallis_systematics.json'}")
    print(f"[out] {figdir / 'm01_tsallis_systematics'}.pdf|png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
