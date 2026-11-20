#!/usr/bin/env python3
"""Module 2, step 8b: compare the pion-origin runs of several PYTHIA settings.

    python analysis/m02_pion_origin_summary.py

Reads results/tables/m02/m02_pion_origin*.json (written by
analysis/m02_pion_origin.py) and puts the pi, K, p MC/data of all settings
side by side.

The number to watch is the pion SHAPE RATIO
    S = [MC/data at 0.4-0.8 GeV/c] / [MC/data at 0.1-0.25 GeV/c]
S = 1 means the MC has the data's shape there whatever its normalisation.
Monash has S = 1.38: MC pions are too few at the lowest pT and too many
around 0.5 GeV/c.  A setting "explains" the pion shape only if it moves S.

chi2 is over 0.1-3 GeV/c with data stat+sys in quadrature treated as
uncorrelated -- a ranking between settings, not a goodness of fit.

Output: results/tables/m02/m02_pion_origin_summary.json,
        results/figures/m02/m02_pion_origin_summary.{pdf,png}
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from hicpy.style import apply_style, save_figure, OKABE_ITO, plt   # noqa: E402

TAB = REPO / "results" / "tables" / "m02"
ORDER = ["", "_udv025", "_udv10", "_udv15", "_bose", "_rescatter", "_pt0240_sig031", "_pt0240_sig031_rescatter"]
SHOW = ["", "_udv15", "_rescatter", "_pt0240_sig031", "_pt0240_sig031_rescatter"]


def main() -> int:
    runs = {}
    for tag in ORDER:
        f = TAB / f"m02_pion_origin{tag}.json"
        if f.exists():
            runs[tag] = json.loads(f.read_text())
    if "" not in runs:
        print("no Monash run (results/tables/m02/m02_pion_origin.json) -- run analysis/m02_pion_origin.py first")
        return 1
    rows = []
    print(f"{'setting':40s} {'pi .1-.25':>9s} {'pi .4-.8':>8s} {'S':>5s} {'chi2':>5s} "
          f"{'K .8-1.5':>8s} {'p .25-.4':>8s} {'p .8-1.5':>8s}")
    for tag, d in runs.items():
        pi, K, p = (d["species"][s]["bands"] for s in ("pi", "K", "p"))
        a, b = pi["0.1-0.25"]["mc_over_data"], pi["0.4-0.8"]["mc_over_data"]
        row = {"tag": tag, "label": d.get("label", tag), "pythia_version": d["pythia_version"],
               "extra_settings": d["extra_settings"], "pi_0.1-0.25": a, "pi_0.4-0.8": b, "shape_ratio": b / a,
               "chi2_pi_0.1-3": d["pion_rescale_diagnostic"]["chi2_unscaled"],
               "n_bins": d["pion_rescale_diagnostic"]["n_bins"],
               "K_0.8-1.5": K["0.8-1.5"]["mc_over_data"], "p_0.25-0.4": p["0.25-0.4"]["mc_over_data"],
               "p_0.8-1.5": p["0.8-1.5"]["mc_over_data"]}
        rows.append(row)
        print(f"{row['label']:40s} {a:9.3f} {b:8.3f} {b / a:5.2f} {row['chi2_pi_0.1-3']:5.0f} "
              f"{row['K_0.8-1.5']:8.3f} {row['p_0.25-0.4']:8.3f} {row['p_0.8-1.5']:8.3f}")

    apply_style()
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8), sharey=True)
    colors = dict(zip(SHOW, [OKABE_ITO[0], OKABE_ITO[2], OKABE_ITO[1], OKABE_ITO[3], OKABE_ITO[4]]))
    for ax, sp in zip(axes, ("pi", "K", "p")):
        ref = runs[""]["species"][sp]
        e = np.array(ref["edges"]); c = 0.5 * (e[1:] + e[:-1])
        dy, dt = np.array(ref["data"]), np.array(ref["data_tot"])
        ax.fill_between(c, 1 - dt / dy, 1 + dt / dy, step="mid", color="0.88", lw=0, label="ALICE unc.")
        for tag in SHOW:
            if tag not in runs:
                continue
            r = runs[tag]["species"][sp]; k = c < 3.0
            ax.step(c[k], (np.array(r["mc_total"]) / dy)[k], where="mid", lw=1.7, color=colors[tag],
                    label=runs[tag].get("label", tag))
        ax.axhline(1, color="k", lw=0.8); ax.axvspan(0.4, 0.8, color=OKABE_ITO[4], alpha=0.10, lw=0)
        ax.set_xscale("log"); ax.set_xlim(0.1, 3.0); ax.set_ylim(0.6, 1.6)
        ticks = [0.1, 0.2, 0.5, 1, 2, 3]
        ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
        ax.set_xlabel(r"$p_{\rm T}$ (GeV/$c$)")
        ax.set_title({"pi": r"$\pi^\pm$", "K": r"K$^\pm$", "p": r"p, $\bar{\rm p}$"}[sp])
    axes[0].set_ylabel("PYTHIA / ALICE, pp 13 TeV, |y| < 0.5")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=9, frameon=False)
    fig.subplots_adjust(bottom=0.30)
    save_figure(fig, REPO / "results" / "figures" / "m02" / "m02_pion_origin_summary"); plt.close(fig)
    (TAB / "m02_pion_origin_summary.json").write_text(json.dumps(rows, indent=1))
    print(f"[out] {TAB / 'm02_pion_origin_summary.json'}, results/figures/m02/m02_pion_origin_summary.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
