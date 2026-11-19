#!/usr/bin/env python3
"""Module 2, step 5: two-parameter scan (MPI pT0Ref x colour-reconnection range)
against ALICE dN/deta AND the pT spectrum at the same time.

    python analysis/m02_tune_scan.py generate [--n-events 10000] [--workers 4]
    python analysis/m02_tune_scan.py fit

Why two parameters and two observables
--------------------------------------
The one-at-a-time study (m02_sensitivity.py) showed that pT0Ref moves the
NUMBER of particles but hardly the SHAPE of the pT spectrum, while colour
reconnection moves the shape.  One observable cannot separate them; two can.
This is a toy version of what a tune is -- not a new tune: two parameters, one
energy, one experiment, and the systematic-correlation model below is an
assumption.

generate
    For every grid point run scripts/m02_run_pythia.py with the two --set
    overrides, compute the observables in the DATA's binning with the same
    functions the comparison scripts use (INEL>0 from the particles,
    dN/deta in 18 eta bins, invariant pT yield in 47 bins), write one JSON per
    point under results/tables/m02/scan/, and delete the ROOT file.

fit
    chi2 per point for dN/deta, for the pT spectrum and for both, under two
    models of the systematic uncertainties, because their bin-to-bin
    correlations are not given in the HEPData tables:
      uncorr   stat and sys added in quadrature, independent bins
      corr     sys fully correlated within each observable (one nuisance
               normalisation per observable, profiled analytically), stat
               independent
    The truth lies between; a conclusion that holds under both is robust.
    MC statistical errors are added to the diagonal in both models.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import pathlib
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "analysis"))
sys.path.insert(0, str(REPO / "src"))

SCAN_DIR = REPO / "results" / "tables" / "m02" / "scan"
GEN_DIR = REPO / "data" / "generated" / "scan"
PT0REF = [2.28, 2.40, 2.52, 2.64, 2.76]          # Monash default 2.28 GeV
CR_RANGE = [1.0, 1.4, 1.8, 2.6]                  # Monash default 1.8
DNDETA_JSON = REPO / "data" / "raw" / "alice_pp13TeV_dndeta.json"
PT_JSON = REPO / "data" / "raw" / "alice_pp13TeV_pt.json"


def tag(pt0, cr):
    return f"pt0_{pt0:.2f}_cr_{cr:.1f}"


def observables(root_file):
    from m02_compare_data import load_dataset, read_mc, select_events, dndeta
    from m02_compare_pt import invariant_yield
    _, d1 = load_dataset(DNDETA_JSON)
    _, d2 = load_dataset(PT_JSON)
    mc = read_mc(pathlib.Path(root_file))
    mask = select_events(mc, "INEL>0")
    e1 = np.concatenate([d1["xlo"], [d1["xhi"][-1]]])
    e2 = np.concatenate([d2["xlo"], [d2["xhi"][-1]]])
    y1, s1 = dndeta(mc, mask, e1, 0.0)
    y2, s2 = invariant_yield(mc, mask, e2)
    return {"n_events": int(len(mask)), "f_inel0": float(mask.mean()),
            "dndeta": y1.tolist(), "dndeta_err": s1.tolist(),
            "pt_yield": y2.tolist(), "pt_yield_err": s2.tolist(),
            "pythia_version": mc["provenance"].get("pythia_version"),
            "settings": mc["provenance"].get("extra_settings")}


def run_point(args):
    pt0, cr, n_events = args
    out_json = SCAN_DIR / f"{tag(pt0, cr)}.json"
    if out_json.exists():
        return f"{tag(pt0, cr)}: exists"
    GEN_DIR.mkdir(parents=True, exist_ok=True)
    root = GEN_DIR / f"{tag(pt0, cr)}.root"
    cmd = [sys.executable, str(REPO / "scripts" / "m02_run_pythia.py"), "--n-events", str(n_events),
           "--out", str(root),
           "--set", f"MultipartonInteractions:pT0Ref = {pt0}",
           "--set", f"ColourReconnection:range = {cr}"]
    log = GEN_DIR / f"{tag(pt0, cr)}.log"
    with open(log, "w") as fh:
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        return f"{tag(pt0, cr)}: generation FAILED (see {log})"
    obs = observables(root)
    obs.update({"pT0Ref": pt0, "cr_range": cr})
    out_json.write_text(json.dumps(obs, indent=1))
    root.unlink()
    root.with_suffix(".json").unlink(missing_ok=True)
    return f"{tag(pt0, cr)}: done"


def generate(n_events, workers):
    SCAN_DIR.mkdir(parents=True, exist_ok=True)
    points = [(p, c, n_events) for p, c in itertools.product(PT0REF, CR_RANGE)]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for msg in ex.map(run_point, points):
            print(msg, flush=True)


def chi2(mc, mc_err, data, stat, sys_, model):
    """chi2 for one observable.  'corr': sys = one fully correlated relative
    normalisation (profiled); 'uncorr': sys added per bin in quadrature."""
    mc, mc_err, data, stat, sys_ = map(np.asarray, (mc, mc_err, data, stat, sys_))
    if model == "uncorr":
        var = stat ** 2 + sys_ ** 2 + mc_err ** 2
        return float(np.sum((mc - data) ** 2 / var)), 1.0
    # data -> data * (1 + eps); eps ~ N(0, s), s = mean relative sys
    var = stat ** 2 + mc_err ** 2
    s = float(np.mean(sys_ / data))
    # minimise sum (mc - d(1+eps))^2/var + eps^2/s^2 over eps (linear)
    a = np.sum(data ** 2 / var) + 1.0 / s ** 2
    b = np.sum((mc - data) * data / var)
    eps = b / a
    c2 = float(np.sum((mc - data * (1 + eps)) ** 2 / var) + (eps / s) ** 2)
    return c2, 1.0 + eps


def fit():
    from m02_compare_data import load_dataset
    _, d1 = load_dataset(DNDETA_JSON)
    _, d2 = load_dataset(PT_JSON)
    pts = [json.loads(p.read_text()) for p in sorted(SCAN_DIR.glob("pt0_*.json"))]
    if not pts:
        print("no scan points; run `generate` first"); return 1
    rows = []
    for p in pts:
        r = {"pT0Ref": p["pT0Ref"], "cr_range": p["cr_range"], "f_inel0": p["f_inel0"]}
        for model in ("uncorr", "corr"):
            c1, n1 = chi2(p["dndeta"], p["dndeta_err"], d1["y"], d1["stat"], d1["sys"], model)
            c2, n2 = chi2(p["pt_yield"], p["pt_yield_err"], d2["y"], d2["stat"], d2["sys"], model)
            r[model] = {"dndeta": c1, "pt": c2, "total": c1 + c2,
                        "norm_dndeta": n1, "norm_pt": n2}
        central = np.abs(0.5 * (d1["xlo"] + d1["xhi"])) < 0.5
        r["dndeta_central_ratio"] = float(np.mean(np.asarray(p["dndeta"])[central] / d1["y"][central]))
        e2 = np.concatenate([d2["xlo"], [d2["xhi"][-1]]]); c = 0.5 * (e2[1:] + e2[:-1])
        w = np.asarray(p["pt_yield"]) * 2 * np.pi * c * np.diff(e2)
        wd = d2["y"] * 2 * np.pi * c * np.diff(e2)
        r["mean_pt_ratio"] = float(((w * c).sum() / w.sum()) / ((wd * c).sum() / wd.sum()))
        rows.append(r)
    n1, n2 = len(d1["y"]), len(d2["y"])
    print(f"{len(rows)} points; bins: dN/deta {n1}, pT {n2}")
    for model in ("uncorr", "corr"):
        print(f"\n== systematics model: {model}  (chi2 dN/deta | pT | total)")
        best = min(rows, key=lambda r: r[model]["total"])
        best1 = min(rows, key=lambda r: r[model]["dndeta"])
        hdr = "pT0Ref \\ CR " + "".join(f"{c:>18.1f}" for c in CR_RANGE)
        print(hdr)
        for p0 in PT0REF:
            line = f"{p0:11.2f} "
            for cr in CR_RANGE:
                r = next(x for x in rows if x["pT0Ref"] == p0 and x["cr_range"] == cr)
                m = r[model]
                line += f"{m['dndeta']:6.0f}|{m['pt']:5.0f}|{m['total']:5.0f}"
            print(line)
        for lab, b in (("best (dN/deta only)", best1), ("best (both)", best)):
            print(f"{lab:20s}: pT0Ref {b['pT0Ref']:.2f}, CR range {b['cr_range']:.1f}; chi2 dN/deta "
                  f"{b[model]['dndeta']:.0f}/{n1}, pT {b[model]['pt']:.0f}/{n2}; dN/deta(|eta|<0.4) "
                  f"MC/data {b['dndeta_central_ratio']:.3f}, <pT> MC/data {b['mean_pt_ratio']:.3f}")
    monash = next(x for x in rows if x["pT0Ref"] == 2.28 and x["cr_range"] == 1.8)
    out = {"grid": {"pT0Ref": PT0REF, "cr_range": CR_RANGE}, "n_bins": {"dndeta": n1, "pt": n2},
           "monash": monash, "points": rows,
           "note": "toy two-parameter scan; systematic correlations bracketed by two models; not a tune"}
    path = REPO / "results" / "tables" / "m02" / "m02_tune_scan.json"
    path.write_text(json.dumps(out, indent=2))

    # heat maps: chi2/ndf for dN/deta, pT and both, uncorrelated-systematics model
    from hicpy.style import apply_style, save_figure, plt
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9), sharey=True)
    for ax, (key, ndf, title) in zip(axes, (("dndeta", n1, r"d$N$/d$\eta$"), ("pt", n2, r"$p_{\rm T}$ spectrum"),
                                             ("total", n1 + n2, "both"))):
        z = np.array([[next(r for r in rows if r["pT0Ref"] == p0 and r["cr_range"] == cr)["uncorr"][key] / ndf
                       for cr in CR_RANGE] for p0 in PT0REF])
        im = ax.imshow(np.log10(z), origin="lower", cmap="viridis_r", aspect="auto")
        for i in range(len(PT0REF)):
            for j in range(len(CR_RANGE)):
                ax.text(j, i, f"{z[i, j]:.1f}", ha="center", va="center", fontsize=9,
                        color="white" if np.log10(z[i, j]) > np.log10(z).mean() else "black")
        ax.set_xticks(range(len(CR_RANGE))); ax.set_xticklabels([f"{c:.1f}" for c in CR_RANGE])
        ax.set_yticks(range(len(PT0REF))); ax.set_yticklabels([f"{p:.2f}" for p in PT0REF])
        ax.tick_params(which="minor", left=False, bottom=False, top=False, right=False)
        ax.set_xlabel("ColourReconnection:range"); ax.set_title(rf"$\chi^2$/ndf, {title}", fontsize=13)
        i0, j0 = PT0REF.index(2.28), CR_RANGE.index(1.8)
        ax.add_patch(plt.Rectangle((j0 - 0.5, i0 - 0.5), 1, 1, fill=False, ec="#D63F34", lw=2))
    axes[0].set_ylabel("MPI pT0Ref (GeV)")
    figdir = REPO / "results" / "figures" / "m02"; figdir.mkdir(parents=True, exist_ok=True)
    save_figure(fig, figdir / "m02_tune_scan"); plt.close(fig)
    print(f"\n[out] {path}, {figdir / 'm02_tune_scan.pdf'} (red box = Monash)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["generate", "fit"])
    ap.add_argument("--n-events", type=int, default=10000)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = ap.parse_args()
    if args.mode == "generate":
        generate(args.n_events, args.workers)
        return 0
    return fit()


if __name__ == "__main__":
    raise SystemExit(main())
