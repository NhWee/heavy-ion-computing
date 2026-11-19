#!/usr/bin/env python3
"""Module 2, step 6: does the string-fragmentation pT width fix the SHAPE?

    python analysis/m02_stringpt_scan.py generate [--n-events 20000] [--workers 4]
    python analysis/m02_stringpt_scan.py fit

After the pT0Ref x CR scan (m02_tune_scan.py) the dN/deta normalisation is
fixed by pT0Ref = 2.40 GeV, but the pT spectrum keeps its shape: MC/data 0.85
at 0.15 GeV/c rising to 1.23 near 0.5 GeV/c.  Neither MPI nor colour
reconnection controls the transverse momentum a hadron gets when a string
breaks.  In the Lund model a q-qbar pair created in a string breaking gets
opposite transverse momenta drawn from a Gaussian of width

    StringPT:sigma   (Monash 0.335 GeV)

so a smaller sigma makes every primary hadron softer.  This scan varies that
one parameter with pT0Ref = 2.40 and Monash colour reconnection held fixed,
and reports the pT-spectrum chi2 (0.15 < pT < 5 GeV/c, both systematic
models), the ratio in a soft and an intermediate pT region, <pT>, and
dN/deta -- sigma also moves the multiplicity a little (softer hadrons, more
of them for the same string energy).
"""
from __future__ import annotations

import argparse
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

from m02_tune_scan import observables, chi2, PT_FIT_MAX, DNDETA_JSON, PT_JSON   # noqa: E402

SIGMAS = [0.26, 0.29, 0.31, 0.335, 0.36, 0.39]     # GeV; Monash 0.335
PT0REF = 2.40
OUT_DIR = REPO / "results" / "tables" / "m02" / "stringpt"
GEN_DIR = REPO / "data" / "generated" / "scan"


def tag(sig):
    return f"stringpt_{sig:.3f}"


def run_point(args):
    sig, n_events = args
    out = OUT_DIR / f"{tag(sig)}.json"
    if out.exists():
        return f"{tag(sig)}: exists"
    GEN_DIR.mkdir(parents=True, exist_ok=True)
    root = GEN_DIR / f"{tag(sig)}.root"
    cmd = [sys.executable, str(REPO / "scripts" / "m02_run_pythia.py"), "--n-events", str(n_events),
           "--out", str(root), "--set", f"MultipartonInteractions:pT0Ref = {PT0REF}",
           "--set", f"StringPT:sigma = {sig}"]
    with open(GEN_DIR / f"{tag(sig)}.log", "w") as fh:
        if subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode != 0:
            return f"{tag(sig)}: generation FAILED"
    obs = observables(root)
    obs.update({"StringPT:sigma": sig, "pT0Ref": PT0REF})
    out.write_text(json.dumps(obs, indent=1))
    root.unlink(); root.with_suffix(".json").unlink(missing_ok=True)
    return f"{tag(sig)}: done"


def fit():
    from m02_compare_data import load_dataset
    _, d1 = load_dataset(DNDETA_JSON)
    _, d2 = load_dataset(PT_JSON)
    pts = sorted((json.loads(p.read_text()) for p in OUT_DIR.glob("stringpt_*.json")),
                 key=lambda r: r["StringPT:sigma"])
    if not pts:
        print("no points; run `generate` first"); return 1
    e2 = np.concatenate([d2["xlo"], [d2["xhi"][-1]]]); c = 0.5 * (e2[1:] + e2[:-1])
    k = d2["xhi"] <= PT_FIT_MAX + 1e-9
    central = np.abs(0.5 * (d1["xlo"] + d1["xhi"])) < 0.5
    wd = d2["y"] * 2 * np.pi * c * np.diff(e2)
    rows = []
    for p in pts:
        y = np.asarray(p["pt_yield"]); ye = np.asarray(p["pt_yield_err"]); r = y / d2["y"]
        w = y * 2 * np.pi * c * np.diff(e2)
        row = {"sigma": p["StringPT:sigma"],
               "dndeta_ratio": float(np.mean(np.asarray(p["dndeta"])[central] / d1["y"][central])),
               "mean_pt_ratio": float(((w * c).sum() / w.sum()) / ((wd * c).sum() / wd.sum())),
               "ratio_0.15-0.25": float(r[(c > 0.15) & (c < 0.25)].mean()),
               "ratio_0.4-0.8": float(r[(c > 0.4) & (c < 0.8)].mean()),
               "ratio_1-4": float(r[(c > 1) & (c < 4)].mean()),
               "shape_spread": float(r[k].max() / r[k].min())}
        for model in ("uncorr", "corr"):
            row[f"chi2_pt_{model}"] = chi2(y[k], ye[k], d2["y"][k], d2["stat"][k], d2["sys"][k], model)[0]
            row[f"chi2_dndeta_{model}"] = chi2(p["dndeta"], p["dndeta_err"], d1["y"], d1["stat"],
                                               d1["sys"], model)[0]
        rows.append(row)
    n2 = int(k.sum())
    print(f"pT0Ref = {PT0REF}, CR as Monash; pT chi2 over {n2} bins (pT < {PT_FIT_MAX:g} GeV/c)")
    print(f"{'sigma':>6s} {'dN/deta':>8s} {'<pT>':>6s} {'0.15-0.25':>9s} {'0.4-0.8':>8s} {'1-4':>6s} "
          f"{'max/min':>7s} {'chi2 pT unc':>11s} {'corr':>6s} {'chi2 dNdeta unc':>15s}")
    for r in rows:
        print(f"{r['sigma']:6.3f} {r['dndeta_ratio']:8.3f} {r['mean_pt_ratio']:6.3f} {r['ratio_0.15-0.25']:9.3f} "
              f"{r['ratio_0.4-0.8']:8.3f} {r['ratio_1-4']:6.3f} {r['shape_spread']:7.3f} "
              f"{r['chi2_pt_uncorr']:11.0f} {r['chi2_pt_corr']:6.0f} {r['chi2_dndeta_uncorr']:15.1f}")
    out = REPO / "results" / "tables" / "m02" / "m02_stringpt_scan.json"
    out.write_text(json.dumps({"pT0Ref": PT0REF, "pt_fit_max_gev": PT_FIT_MAX, "points": rows,
                               "pythia_version": pts[0].get("pythia_version")}, indent=2))

    from hicpy.style import apply_style, save_figure, OKABE_ITO, plt
    apply_style()
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for i, r in enumerate(rows):
        p = pts[i]; rr = np.asarray(p["pt_yield"]) / d2["y"]
        ax.plot(c[k], rr[k], "-", lw=1.6, color=plt.cm.viridis(i / max(len(rows) - 1, 1)),
                label=rf"$\sigma$ = {r['sigma']:.3f}" + (" (Monash)" if abs(r["sigma"] - 0.335) < 1e-6 else ""))
    ax.axhline(1, color="0.4", lw=1)
    ax.set_xscale("log"); ax.set_xlabel(r"$p_{\rm T}$ (GeV/$c$)"); ax.set_ylabel("MC / ALICE")
    from matplotlib.ticker import FixedLocator, NullFormatter, FormatStrFormatter
    ax.xaxis.set_major_locator(FixedLocator([0.2, 0.5, 1, 2, 5]))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g")); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(0.7, 1.4); ax.legend(fontsize=9, ncol=2, loc="upper left")
    ax.set_title(f"pT0Ref = {PT0REF} GeV, |η| < 0.8, INEL>0", fontsize=11)
    figdir = REPO / "results" / "figures" / "m02"; figdir.mkdir(parents=True, exist_ok=True)
    save_figure(fig, figdir / "m02_stringpt_scan"); plt.close(fig)
    print(f"[out] {out}, {figdir / 'm02_stringpt_scan.pdf'}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["generate", "fit"])
    ap.add_argument("--n-events", type=int, default=20000)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    if a.mode == "generate":
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for msg in ex.map(run_point, [(s, a.n_events) for s in SIGMAS]):
                print(msg, flush=True)
        return 0
    return fit()


if __name__ == "__main__":
    raise SystemExit(main())
