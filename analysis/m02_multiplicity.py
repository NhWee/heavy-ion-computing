#!/usr/bin/env python3
"""Module 2, step 9: particle ratios and <pT> versus multiplicity, PYTHIA vs ALICE.

    python analysis/m02_multiplicity.py generate [--n-events 1000000] [--workers N]
                                                 [--set "Key = value" ...] [--tag _x] [--label "..."]
    python analysis/m02_multiplicity.py compare  [--tags "" _rescatter ...]

QUESTION
--------
Two observations of section 13 of docs/modules/module02_pythia.md need a
multiplicity-dependent test:
  - K/pi is 7-9 % low in minimum bias.  ALICE sees strangeness ENHANCEMENT:
    K/pi, Lambda/pi, Xi/pi, Omega/pi rise with multiplicity, the more the more
    strange quarks the hadron has.  Is PYTHIA's deficit the same at every
    multiplicity (a normalisation) or growing (missing physics)?
  - hadronic rescattering improved the pion shape.  Its effect must GROW with
    the hadron density, i.e. with multiplicity; <pT>(dNch/deta) per species is
    where that shows.

DATA (HEPData ins1784041, ALICE arXiv:2003.02394; tables in data/raw/ins1784041/)
  Table 17  (K+ + K-)/(pi+ + pi-), (p + pbar)/(pi+ + pi-), pT-integrated
  Table 18  2 K0S / pi, (Lambda + Lambdabar)/pi, (Xi- + Xibar+)/pi
  Table 19  (Omega- + Omegabar+)/pi          (5 merged classes)
  Table 20  <pT> of pi, K, p
  x axis: <dNch/deta> in |eta| < 0.5 of each V0M class (the table gives a
  low-high band; the centre is used).  All in |y| < 0.5, INEL>0.

MC EVENT CLASSES, as ALICE defines them
  - INEL>0: >= 1 charged final-state particle in |eta| < 1
  - V0M estimator: charged final-state particles in the V0 acceptance,
    V0A 2.8 < eta < 5.1 and V0C -3.7 < eta < -1.7 (forward, so the class
    selection does not bias the mid-rapidity yields directly)
  - classes I..X = percentiles of the INEL>0 sample in V0M: 0-0.92, 0.92-4.6,
    4.6-9.2, 9.2-13.8, 13.8-18.4, 18.4-27.6, 27.6-36.8, 36.8-46, 46-64.5,
    64.5-100 % (ties in the integer V0M count broken at random)
  Each MC class is then placed at ITS OWN <dNch/deta>, so the comparison is
  at equal mid-rapidity multiplicity even though PYTHIA's overall dN/deta
  is 6 % high.  Strange hadrons: K0S, Lambda, Xi, Omega all have c*tau > 1 cm
  and are not decayed (card: tau0Max = 10 mm) -- like ALICE's yields, which
  are corrected for feed-down.

OUTPUT
  data/generated/m02_mult{tag}.npz                 per-event counts (not in git)
  results/tables/m02/m02_mult{tag}.json            per-class MC values
  results/tables/m02/m02_mult_compare.json         data and all MC sets
  results/figures/m02/m02_mult_ratios.{pdf,png}, m02_mult_meanpt.{pdf,png}
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

GEN, TAB, FIG = REPO / "data" / "generated", REPO / "results" / "tables" / "m02", REPO / "results" / "figures" / "m02"
RAW = REPO / "data" / "raw" / "ins1784041"
CLASS_EDGES = [0, 0.92, 4.6, 9.2, 13.8, 18.4, 27.6, 36.8, 46, 64.5, 100]
CLASSES = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
SPECIES = {"pi": (211,), "K": (321,), "p": (2212,), "K0S": (310,), "Lambda": (3122,), "Xi": (3312,), "Omega": (3334,)}
PT_SPECIES = ("pi", "K", "p")
SEED0 = 20261120
# ratio name -> (numerator species, factor, data table, dependent index)
RATIOS = {"K/pi": ("K", 1.0, 17, 0), "p/pi": ("p", 1.0, 17, 1), "2K0S/pi": ("K0S", 2.0, 18, 0),
          "Lambda/pi": ("Lambda", 1.0, 18, 1), "Xi/pi": ("Xi", 1.0, 18, 2), "Omega/pi": ("Omega", 1.0, 19, 0)}
MEANPT = {"pi": (20, 0), "K": (20, 1), "p": (20, 2)}


# ----------------------------------------------------------------------------- generation
def _chunk(args):
    n_events, seed, settings, card = args
    try:
        import pythia8
    except ImportError:
        import pythia8mc as pythia8
    py = pythia8.Pythia("", False)            # quiet construction
    py.readFile(card)
    for s in ("Random:setSeed = on", f"Random:seed = {seed}", "Next:numberCount = 0", "Print:quiet = on"):
        py.readString(s)
    for s in settings:
        if not py.readString(s):
            raise RuntimeError(f"PYTHIA rejected {s!r}")
    if not py.init():
        raise RuntimeError("PYTHIA init() failed")
    code = {c: k for k, v in SPECIES.items() for c in v}
    rows = []
    n_acc = n_tried = 0
    while n_acc < n_events and n_tried < 2 * n_events + 100:
        n_tried += 1
        if not py.next():
            continue
        n_acc += 1
        nv0 = nch1 = nch05 = 0
        cnt = dict.fromkeys(SPECIES, 0); spt = dict.fromkeys(PT_SPECIES, 0.0)
        ev = py.event
        for i in range(ev.size()):
            p = ev[i]
            if not p.isFinal():
                continue
            if p.isCharged():
                eta = p.eta()
                if abs(eta) < 1.0:
                    nch1 += 1
                    if abs(eta) < 0.5:
                        nch05 += 1
                if 2.8 < eta < 5.1 or -3.7 < eta < -1.7:
                    nv0 += 1
            sp = code.get(abs(p.id()))
            if sp is not None and abs(p.y()) < 0.5:
                cnt[sp] += 1
                if sp in spt:
                    spt[sp] += p.pT()
        rows.append([nv0, nch1, nch05] + [cnt[s] for s in SPECIES] + [spt[s] for s in PT_SPECIES])
    return np.array(rows, float), py.settings.parm("Pythia:versionNumber")


def generate(args) -> int:
    card = str(REPO / "generators" / "pythia_mb_13TeV.cmnd")
    w = args.workers or os.cpu_count() or 1
    chunks = [(args.n_events // w + (i < args.n_events % w), SEED0 + i, args.set, card) for i in range(w)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=w) as ex:
        parts = list(ex.map(_chunk, chunks))
    data = np.vstack([p[0] for p in parts]); version = parts[0][1]
    cols = ["nv0", "nch1", "nch05"] + [f"n_{s}" for s in SPECIES] + [f"sumpt_{s}" for s in PT_SPECIES]
    GEN.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(GEN / f"m02_mult{args.tag}.npz", data=data, columns=np.array(cols),
                        settings=np.array(args.set, dtype=object), version=version, allow_pickle=True)
    print(f"[gen] {len(data)} events in {time.time() - t0:.0f} s on {w} workers, PYTHIA {version:.3f}, "
          f"settings {args.set or 'none'}")
    res = classify(data, cols)
    res.update({"label": args.label, "pythia_version": round(version, 3), "extra_settings": args.set,
                "n_events": int(len(data)), "seeds": [c[1] for c in chunks]})
    TAB.mkdir(parents=True, exist_ok=True)
    (TAB / f"m02_mult{args.tag}.json").write_text(json.dumps(res, indent=1))
    for c in res["classes"]:
        print(f"  {c['class']:>4s}: <dNch/deta> {c['dndeta']:6.2f}  K/pi {c['ratios']['K/pi']:.4f}  "
              f"p/pi {c['ratios']['p/pi']:.4f}  Xi/pi {c['ratios']['Xi/pi']:.5f}  <pT>pi {c['meanpt']['pi']:.3f}")
    print(f"[out] {TAB / f'm02_mult{args.tag}.json'}")
    return 0


def classify(data, cols):
    """V0M percentile classes of the INEL>0 events; per class pooled ratios."""
    col = {c: i for i, c in enumerate(cols)}
    inel0 = data[data[:, col["nch1"]] > 0]
    rng = np.random.default_rng(1)
    key = inel0[:, col["nv0"]] + rng.random(len(inel0))          # random tie-break
    order = np.argsort(-key, kind="stable")                         # most active first
    n = len(order); out = []
    for name, lo, hi in zip(CLASSES, CLASS_EDGES[:-1], CLASS_EDGES[1:]):
        sel = inel0[order[int(round(n * lo / 100)):int(round(n * hi / 100))]]
        npi = sel[:, col["n_pi"]].sum()
        rec = {"class": name, "percentile": [lo, hi], "n_events": int(len(sel)),
               "dndeta": float(sel[:, col["nch05"]].mean()),          # |eta| < 0.5 -> width 1
               "ratios": {}, "ratio_stat": {}, "meanpt": {}, "counts": {}}
        for r, (sp, fac, _, _) in RATIOS.items():
            nx = sel[:, col[f"n_{sp}"]].sum()
            rec["ratios"][r] = float(fac * nx / npi) if npi else None
            rec["ratio_stat"][r] = float(fac * np.sqrt(nx) / npi) if npi else None
            rec["counts"][sp] = int(nx)
        rec["counts"]["pi"] = int(npi)
        for sp in PT_SPECIES:
            nsp = sel[:, col[f"n_{sp}"]].sum()
            rec["meanpt"][sp] = float(sel[:, col[f"sumpt_{sp}"]].sum() / nsp) if nsp else None
        out.append(rec)
    return {"n_inel0": int(n), "classes": out,
            "dndeta_inel0": float(inel0[:, col["nch05"]].mean())}


# ----------------------------------------------------------------------------- data
def read_table(n, dep):
    doc = yaml.safe_load((RAW / f"Table_{n}.yaml").read_text())
    x = np.array([0.5 * (v["low"] + v["high"]) for v in doc["independent_variables"][0]["values"]])
    vals = doc["dependent_variables"][dep]["values"]
    y = np.array([float(v["value"]) for v in vals])
    err = {lab: np.array([next((float(e["symerror"]) for e in v.get("errors", []) if e["label"] == lab), 0.0)
                          for v in vals]) for lab in ("stat", "full syst.", "unc. syst.")}
    return {"x": x, "y": y, "stat": err["stat"], "sys": err["full syst."], "sys_unc": err["unc. syst."]}


def mc_points(mc, quantity, merge_pairs=False):
    """(dN/deta, value, relative stat. error) per MC class.  merge_pairs: classes
    I+II, III+IV, ... as in ALICE's Omega table (Table 19), for statistics."""
    cls = mc["classes"]
    if not merge_pairs:
        x = np.array([c["dndeta"] for c in cls])
        y = np.array([c[quantity[0]][quantity[1]] for c in cls], float)
        if quantity[0] == "ratios":
            e = np.array([c["ratio_stat"][quantity[1]] for c in cls], float) / np.where(y > 0, y, np.nan)
        else:
            e = np.zeros_like(y)
        return x, y, e
    sp, fac = RATIOS[quantity[1]][:2]
    xs, ys, es = [], [], []
    for a, b in zip(cls[0::2], cls[1::2]):
        n = a["n_events"] + b["n_events"]
        nx = a["counts"][sp] + b["counts"][sp]; npi = a["counts"]["pi"] + b["counts"]["pi"]
        xs.append((a["dndeta"] * a["n_events"] + b["dndeta"] * b["n_events"]) / n)
        ys.append(fac * nx / npi); es.append(1 / np.sqrt(nx) if nx else np.nan)
    return np.array(xs), np.array(ys), np.array(es)


def interp_mc(mc, quantity, x, merge_pairs=False, max_extrap=0.10):
    """MC value (and relative stat. error) at the data's <dNch/deta>: linear in
    dNch/deta between MC points; outside the MC range linear extrapolation from
    the two end points, allowed up to max_extrap (relative) beyond the range."""
    xs, ys, es = mc_points(mc, quantity, merge_pairs)
    o = np.argsort(xs); xs, ys, es = xs[o], ys[o], es[o]
    x = np.atleast_1d(np.asarray(x, float))
    val = np.interp(x, xs, ys); err = np.interp(x, xs, es)
    lo, hi = x < xs[0], x > xs[-1]
    val[lo] = ys[0] + (ys[1] - ys[0]) / (xs[1] - xs[0]) * (x[lo] - xs[0])
    val[hi] = ys[-1] + (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]) * (x[hi] - xs[-1])
    err[lo], err[hi] = es[0], es[-1]
    bad = (x < xs[0] * (1 - max_extrap)) | (x > xs[-1] * (1 + max_extrap))
    val[bad] = np.nan
    return val, err


# ----------------------------------------------------------------------------- compare
def compare(args) -> int:
    from hicpy.style import apply_style, save_figure, OKABE_ITO, plt
    mcs = {}
    for t in args.tags:
        f = TAB / f"m02_mult{t}.json"
        if f.exists():
            mcs[t] = json.loads(f.read_text())
    if not mcs:
        print("no MC sets: run `generate` first"); return 1
    colors = [OKABE_ITO[i % len(OKABE_ITO)] for i in (0, 1, 2, 3, 4, 5)]
    summary = {"data": "HEPData ins1784041 (ALICE arXiv:2003.02394)", "ratios": {}, "meanpt": {}}

    def one(ax, data, key, title, scale=1.0):
        ax.errorbar(data["x"], data["y"] * scale, yerr=data["stat"] * scale, fmt="o", color="k", ms=4, zorder=5,
                    label="ALICE")
        for x, y, s in zip(data["x"], data["y"] * scale, data["sys"] * scale):
            ax.add_patch(plt.Rectangle((x * 0.96, y - s), x * 0.08, 2 * s, fc="none", ec="k", lw=0.8))
        res = {}
        merge = key == ("ratios", "Omega/pi")
        ih, il = np.argmax(data["x"]), np.argmin(data["x"])
        d_rise = data["y"][ih] / data["y"][il]
        d_rise_err = d_rise * np.hypot(np.hypot(data["stat"][ih], data["sys_unc"][ih]) / data["y"][ih],
                                       np.hypot(data["stat"][il], data["sys_unc"][il]) / data["y"][il])
        for (t, mc), col in zip(mcs.items(), colors):
            xs, ys, es = mc_points(mc, key, merge)
            ax.errorbar(xs, ys * scale, yerr=ys * es * scale, fmt="-s", color=col, lw=1.6, ms=3, capsize=0,
                        label=mc.get("label", t))
            v, e = interp_mc(mc, key, data["x"], merge)
            at = v / data["y"]
            rise = v[ih] / v[il]
            res[t] = {"mc_over_data_lowest": float(at[il]), "mc_over_data_highest": float(at[ih]),
                      "data_rise_high_over_low": float(d_rise), "data_rise_err": float(d_rise_err),
                      "mc_rise_high_over_low": float(rise), "mc_rise_stat_err": float(rise * np.hypot(e[ih], e[il])),
                      "mc_rel_stat_highest": float(e[ih])}
        ax.set_xscale("log"); ax.set_title(title, fontsize=12)
        from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
        ticks = [2, 3, 5, 10, 20, 30]
        ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FixedFormatter([str(t) for t in ticks])); ax.set_xlim(2, 32)
        ax.set_xlabel(r"$\langle dN_{\rm ch}/d\eta\rangle_{|\eta|<0.5}$")
        return res

    apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.6))
    titles = {"K/pi": r"(K$^+$+K$^-$)/($\pi^+$+$\pi^-$)", "p/pi": r"(p+$\bar{\rm p}$)/($\pi^+$+$\pi^-$)",
              "2K0S/pi": r"2K$^0_S$/($\pi^+$+$\pi^-$)", "Lambda/pi": r"($\Lambda$+$\bar\Lambda$)/($\pi^+$+$\pi^-$) $\times 10$",
              "Xi/pi": r"($\Xi$+$\bar\Xi$)/($\pi^+$+$\pi^-$) $\times 10^2$",
              "Omega/pi": r"($\Omega$+$\bar\Omega$)/($\pi^+$+$\pi^-$) $\times 10^3$"}
    scale = {"Lambda/pi": 10, "Xi/pi": 1e2, "Omega/pi": 1e3}
    for ax, (r, (_, _, tab, dep)) in zip(axes.flat, RATIOS.items()):
        d = read_table(tab, dep)
        summary["ratios"][r] = one(ax, d, ("ratios", r), titles[r], scale.get(r, 1.0))
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.legend(h, l, loc="lower center", ncol=5, fontsize=9, frameon=False)
    save_figure(fig, FIG / "m02_mult_ratios"); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8))
    for ax, (sp, (tab, dep)) in zip(axes, MEANPT.items()):
        d = read_table(tab, dep)
        summary["meanpt"][sp] = one(ax, d, ("meanpt", sp), {"pi": r"$\langle p_{\rm T}\rangle$ $\pi^\pm$",
                                                             "K": r"$\langle p_{\rm T}\rangle$ K$^\pm$",
                                                             "p": r"$\langle p_{\rm T}\rangle$ p, $\bar{\rm p}$"}[sp])
        ax.set_ylabel("GeV/$c$")
    h, l = axes[0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.13, 1, 1)); fig.legend(h, l, loc="lower center", ncol=5, fontsize=9, frameon=False)
    save_figure(fig, FIG / "m02_mult_meanpt"); plt.close(fig)

    summary["mc_sets"] = {t: {"label": m.get("label"), "pythia_version": m["pythia_version"],
                              "extra_settings": m["extra_settings"], "n_events": m["n_events"],
                              "dndeta_inel0": m["dndeta_inel0"]} for t, m in mcs.items()}
    (TAB / "m02_mult_compare.json").write_text(json.dumps(summary, indent=1))
    print(f"{'':12s}" + "".join(f"{mcs[t].get('label', t)[:22]:>24s}" for t in mcs))
    for group in ("ratios", "meanpt"):
        for q, res in summary[group].items():
            cells = "".join(f"{res[t]['mc_over_data_lowest']:7.3f} {res[t]['mc_over_data_highest']:6.3f} "
                            f"{res[t]['mc_rise_high_over_low']:5.2f}±{res[t]['mc_rise_stat_err']:.2f} " for t in mcs)
            dr = next(iter(res.values()))
            print(f"{(q if group == 'ratios' else '<pT> ' + q):12s}{cells}  data rise {dr['data_rise_high_over_low']:.2f}"
                  f"±{dr['data_rise_err']:.2f}")
    print("columns per set: MC/data at lowest and highest multiplicity, MC rise high/low ± MC stat;"
          " data rise error = stat (+) multiplicity-uncorrelated syst.")
    print(f"[out] {TAB / 'm02_mult_compare.json'}, {FIG / 'm02_mult_ratios.pdf'}, {FIG / 'm02_mult_meanpt.pdf'}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate")
    g.add_argument("--n-events", type=int, default=1000000)
    g.add_argument("--workers", type=int, default=None)
    g.add_argument("--set", action="append", default=[])
    g.add_argument("--tag", default="")
    g.add_argument("--label", default="Monash")
    c = sub.add_parser("compare")
    c.add_argument("--tags", nargs="+", default=["", "_rescatter", "_pt0240_sig031_rescatter", "_rope"])
    args = ap.parse_args()
    return generate(args) if args.cmd == "generate" else compare(args)


if __name__ == "__main__":
    raise SystemExit(main())
