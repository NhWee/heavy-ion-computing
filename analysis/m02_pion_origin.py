#!/usr/bin/env python3
"""Module 2, step 8: where do PYTHIA's pions (and kaons, protons) come from?

    python analysis/m02_pion_origin.py [--n-events 20000] [--set "Key = value" ...] [--tag _x]
                                       [--label "Monash"]

QUESTION
--------
Section 12 of docs/modules/module02_pythia.md: the PYTHIA pi+- spectrum is
-18 % at 0.1-0.25 GeV/c and +17 % at 0.4-0.8 GeV/c with respect to ALICE --
its peak sits at higher pT.  A final-state pion is either made directly in
the string break-up ("direct") or is a decay product of a heavier hadron
(rho, omega, eta, K*, Delta, ...).  Decay pions have a characteristic pT
shape (rho -> pi pi gives each pion p* = 0.36 GeV/c in the rho frame), so
"too many decay pions" and "direct pions too hard" look different.  This
script splits the MC spectrum by origin and asks which part has the shape
of the excess.

HOW THE ORIGIN IS READ
----------------------
PYTHIA status codes (the |status| of a final-state particle):
  81-89  produced by hadronization (string or ministring break-up) = direct
  91-99  produced in a decay; mother1() is the decaying hadron
The immediate mother is used (pi from rho from a1 counts as "rho").  Decays
happen only for c*tau < 10 mm (card: ParticleDecays:tau0Max), so K0S and
Lambda are stable, as in ALICE's primary-particle definition.

CLOSURE
-------
Same card, same seed, same event loop as scripts/m02_run_pythia.py, so the
events are THE SAME as in data/generated/m02_pythia_mb_13TeV.root (for the
same PYTHIA build and --set options).  The sum over origins must reproduce
the total MC spectrum of analysis/m02_compare_identified.py; the script
prints that comparison when the summary JSON of that run exists.

OUTPUT
------
results/tables/m02/m02_pion_origin{tag}.json   spectra per species x origin
results/figures/m02/m02_pion_origin{tag}.{pdf,png}
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))
sys.path.insert(0, str(REPO / "scripts"))

from hicpy.style import apply_style, save_figure, plt   # noqa: E402
from m02_compare_data import load_dataset                # noqa: E402
from m02_run_pythia import card_value, git_hash          # noqa: E402

Y_MAX = 0.5
SPECIES = {"pi": 211, "K": 321, "p": 2212}
DELTA = {1114, 2114, 2214, 2224}
ORIGINS = {
    "pi": ["direct", "rho", "omega", "eta/eta'", "K*", "Delta", "other baryon", "other meson", "heavy flavour"],
    "K": ["direct", "K*", "phi", "other meson", "baryon", "heavy flavour"],
    "p": ["direct", "Delta", "other baryon", "heavy flavour", "other"],
}
COLORS = {"direct": "#0072B2", "rho": "#D55E00", "omega": "#E69F00", "eta/eta'": "#CC79A7",
          "K*": "#009E73", "Delta": "#56B4E9", "other baryon": "#8E7CC3", "other meson": "#999999",
          "heavy flavour": "#000000", "phi": "#D55E00", "baryon": "#56B4E9", "other": "#999999"}


def is_heavy(aid: int) -> bool:
    """Charm or bottom hadron (PDG numbering: quark content in the digits)."""
    if aid > 1000:
        return (aid // 1000) % 10 in (4, 5)
    return (aid // 100) % 10 in (4, 5) or aid in (441, 443, 553)


def classify(species: str, status: int, mother_id: int) -> str:
    st, m = abs(status), abs(mother_id)
    if 81 <= st <= 89:
        return "direct"
    if not 91 <= st <= 99:
        return "other meson" if species != "p" else "other"
    if is_heavy(m):
        return "heavy flavour"
    if species == "pi":
        if m in (113, 213):
            return "rho"
        if m == 223:
            return "omega"
        if m in (221, 331):
            return "eta/eta'"
        if m in (313, 323):
            return "K*"
        if m in DELTA:
            return "Delta"
        return "other baryon" if m > 1000 else "other meson"
    if species == "K":
        if m in (313, 323):
            return "K*"
        if m == 333:
            return "phi"
        return "baryon" if m > 1000 else "other meson"
    if m in DELTA:
        return "Delta"
    return "other baryon" if m > 1000 else "other"


def generate(args, edges):
    try:
        import pythia8
    except ImportError:
        import pythia8mc as pythia8
    card = pathlib.Path(args.card)
    text = card.read_text()
    py = pythia8.Pythia()                     # same steps as scripts/m02_run_pythia.py
    if not py.readFile(str(card)):
        raise SystemExit(f"PYTHIA rejected the card {card}")
    n_events = args.n_events or int(card_value(text, "Main:numberOfEvents", 1000))
    seed = int(card_value(text, "Random:seed", 0))
    for s in (f"Main:numberOfEvents = {n_events}", "Random:setSeed = on", f"Random:seed = {seed}",
              "Next:numberCount = 0"):
        py.readString(s)
    for extra in args.set:
        if not py.readString(extra):
            raise SystemExit(f"PYTHIA rejected --set {extra!r}")
    if not py.init():
        raise SystemExit("PYTHIA init() failed")

    sumw = {sp: {o: np.zeros(len(edges[sp]) - 1) for o in ORIGINS[sp]} for sp in SPECIES}
    sumpt = {sp: {o: [0.0, 0] for o in ORIGINS[sp]} for sp in SPECIES}
    code_to_sp = {v: k for k, v in SPECIES.items()}
    n_acc = n_tried = n_inel0 = 0
    t0 = time.time()
    while n_acc < n_events:
        n_tried += 1
        if n_tried > 2 * n_events + 100:
            break
        if not py.next():
            continue
        n_acc += 1
        ev = py.event
        inel0, found = False, []
        for i in range(ev.size()):
            p = ev[i]
            if not p.isFinal():
                continue
            if not inel0 and p.isCharged() and abs(p.eta()) < 1.0:
                inel0 = True
            sp = code_to_sp.get(abs(p.id()))
            if sp is None or abs(p.y()) >= Y_MAX:
                continue
            found.append((sp, p.pT(), p.status(), ev[p.mother1()].id() if p.mother1() > 0 else 0))
        if not inel0:
            continue
        n_inel0 += 1
        for sp, pt, st, mid in found:
            o = classify(sp, st, mid)
            k = np.searchsorted(edges[sp], pt, side="right") - 1
            if 0 <= k < len(edges[sp]) - 1:
                sumw[sp][o][k] += 1.0
            sumpt[sp][o][0] += pt; sumpt[sp][o][1] += 1
    info = py.infoPython() if hasattr(py, "infoPython") else py.info
    print(f"[gen] {n_acc} events ({n_tried} tried), {n_inel0} INEL>0, {time.time() - t0:.0f} s; "
          f"PYTHIA {py.settings.parm('Pythia:versionNumber'):.3f}, sigma_gen {info.sigmaGen():.4f} mb")
    return sumw, sumpt, n_acc, n_inel0, py.settings.parm("Pythia:versionNumber")


def plot(out, stem):
    """Per species: origin contributions / data (stacked), and origin fractions."""
    apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.2), sharex="col", gridspec_kw={"height_ratios": [1.3, 1]})
    for j, sp in enumerate(SPECIES):
        r = out["species"][sp]; e = np.array(r["edges"]); c = 0.5 * (e[1:] + e[:-1])
        dy = np.array(r["data"]); dt = np.array(r["data_tot"])
        ax = axes[0, j]; base = np.zeros_like(c)
        for o in ORIGINS[sp]:
            v = np.array(r["origins"][o]) / dy
            ax.fill_between(c, base, base + v, step="mid", color=COLORS[o], alpha=0.85, lw=0, label=o)
            base = base + v
        ax.fill_between(c, 1 - dt / dy, 1 + dt / dy, step="mid", facecolor="none", hatch="///", edgecolor="0.35", lw=0, label="ALICE unc.")
        ax.axhline(1, color="k", lw=1); ax.axvspan(0.4, 0.8, color="0.85", alpha=0.35, lw=0, zorder=-1)
        ax.set_ylim(0, 2.0); ax.set_title({"pi": r"$\pi^\pm$", "K": r"K$^\pm$", "p": r"p, $\bar{\rm p}$"}[sp])
        ax.legend(fontsize=7.5, ncol=3, loc="upper center", frameon=False)
        if not out.get("origins_valid", True):
            ax.text(0.03, 0.04, "origins not valid: status codes\nrewritten by BE / rescattering",
                    transform=ax.transAxes, fontsize=8, color="crimson")
        ax2 = axes[1, j]; tot = np.array(r["mc_total"])
        for o in ORIGINS[sp]:
            f = np.divide(np.array(r["origins"][o]), tot, out=np.zeros_like(tot), where=tot > 0)
            if f.max() > 0.02:
                ax2.plot(c, f, color=COLORS[o], lw=1.6, label=o)
        ax2.set_ylim(0, 1); ax2.set_xscale("log"); ax2.set_xlabel(r"$p_{\rm T}$ (GeV/$c$)")
        ax2.axvspan(0.4, 0.8, color="0.85", alpha=0.35, lw=0, zorder=-1)
    axes[0, 0].set_ylabel(f"{out.get('label', 'PYTHIA')}: origin / ALICE (stacked)")
    axes[1, 0].set_ylabel("fraction of MC")
    for ax in axes[0]:
        ax.set_xscale("log")
    save_figure(fig, stem); plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--card", default=str(REPO / "generators" / "pythia_mb_13TeV.cmnd"))
    ap.add_argument("--n-events", type=int, default=None)
    ap.add_argument("--set", action="append", default=[], metavar='"Key = value"')
    ap.add_argument("--tag", default="")
    ap.add_argument("--label", default="Monash")
    ap.add_argument("--closure-ref", default=None,
                    help="m02_identified summary JSON of the same events (default: the one for --tag)")
    ap.add_argument("--replot", action="store_true", help="redraw the figure from the JSON only")
    ap.add_argument("--data", nargs="+", default=[f"{s}=data/raw/alice_pp13TeV_{s.lower()}.json"
                                                  for s in SPECIES])
    args = ap.parse_args()
    figdir, tabdir = REPO / "results" / "figures" / "m02", REPO / "results" / "tables" / "m02"
    figdir.mkdir(parents=True, exist_ok=True); tabdir.mkdir(parents=True, exist_ok=True)
    if args.replot:
        out = json.loads((tabdir / f"m02_pion_origin{args.tag}.json").read_text())
        plot(out, figdir / f"m02_pion_origin{args.tag}")
        print(f"[out] {figdir / f'm02_pion_origin{args.tag}.pdf'} (replotted)")
        return 0

    data, edges = {}, {}
    for spec in args.data:
        sp, path = spec.split("=", 1)
        _, d = load_dataset(REPO / path)
        data[sp] = d; edges[sp] = np.concatenate([d["xlo"], [d["xhi"][-1]]])

    sumw, sumpt, n_acc, n_inel0, version = generate(args, edges)

    out = {"pythia_version": round(version, 3), "extra_settings": args.set, "n_events": n_acc,
           "n_inel0": n_inel0, "y_max": Y_MAX, "git_commit": git_hash(REPO), "species": {},
           # Bose-Einstein shifts and hadronic rescattering give final particles new status
           # codes / mothers: the totals stay right, the origin split does not.
           "origins_valid": not any(("BoseEinstein" in x or "Rescatter" in x) and "on" in x.split("=")[-1].lower()
                                    for x in args.set)}
    bands = [(0.1, 0.25), (0.25, 0.4), (0.4, 0.8), (0.8, 1.5), (1.5, 3.0)]
    for sp in SPECIES:
        e = edges[sp]; c = 0.5 * (e[1:] + e[:-1]); w = np.diff(e)
        norm = n_inel0 * w * (2 * Y_MAX)
        spec = {o: sumw[sp][o] / norm for o in ORIGINS[sp]}
        err2 = {o: sumw[sp][o] / norm ** 2 for o in ORIGINS[sp]}
        tot = sum(spec.values()); tot_err = np.sqrt(sum(err2.values()))
        dy = data[sp]["y"]; dtot = np.hypot(data[sp]["stat"], data[sp]["sys"])
        rec = {"edges": e.tolist(), "data": dy.tolist(), "data_tot": dtot.tolist(),
               "mc_total": tot.tolist(), "mc_total_err": tot_err.tolist(),
               "origins": {o: spec[o].tolist() for o in ORIGINS[sp]},
               "mean_pt": {o: (sumpt[sp][o][0] / sumpt[sp][o][1] if sumpt[sp][o][1] else None) for o in ORIGINS[sp]},
               "yield_share": {o: float(sumpt[sp][o][1] / max(1, sum(v[1] for v in sumpt[sp].values())))
                               for o in ORIGINS[sp]},
               "bands": {}}
        for lo, hi in bands:
            k = (c > lo) & (c < hi)
            if not k.any():
                continue
            d_int = float(np.sum(dy[k] * w[k])); m_int = float(np.sum(tot[k] * w[k]))
            rec["bands"][f"{lo}-{hi}"] = {
                "mc_over_data": m_int / d_int,
                "fraction": {o: float(np.sum(spec[o][k] * w[k]) / m_int) for o in ORIGINS[sp]},
                "contribution_to_ratio": {o: float(np.sum(spec[o][k] * w[k]) / d_int) for o in ORIGINS[sp]}}
        out["species"][sp] = rec

    # --- diagnostic: rescale direct and decay pions separately (pi, 0.1-3 GeV/c)
    r = out["species"]["pi"]; e = np.array(r["edges"]); c = 0.5 * (e[1:] + e[:-1])
    k = (c > 0.1) & (c < 3.0)
    direct = np.array(r["origins"]["direct"])[k]
    decay = np.array(r["mc_total"])[k] - direct
    dy, dtot = np.array(r["data"])[k], np.array(r["data_tot"])[k]
    A = np.vstack([direct / dtot, decay / dtot]).T
    coef, *_ = np.linalg.lstsq(A, dy / dtot, rcond=None)
    chi2_fit = float(np.sum((A @ coef - dy / dtot) ** 2))
    chi2_raw = float(np.sum(((direct + decay - dy) / dtot) ** 2))
    out["pion_rescale_diagnostic"] = {
        "range_gev": [0.1, 3.0], "scale_direct": float(coef[0]), "scale_decay": float(coef[1]),
        "chi2_unscaled": chi2_raw, "chi2_rescaled": chi2_fit, "n_bins": int(k.sum()),
        "note": "data stat+sys in quadrature, treated as uncorrelated: a shape diagnostic, not a tune"}

    # --- closure against the ROOT-based comparison
    ref = pathlib.Path(args.closure_ref) if args.closure_ref else tabdir / f"m02_identified{args.tag}_summary.json"
    if ref.exists():
        refj = json.loads(ref.read_text())
        same = (refj.get("extra_settings", []) == args.set
                and str(refj.get("pythia_version")) == f"{version:.3f}")
        if not same:
            print(f"[closure] skipped: {ref.name} was made with PYTHIA {refj.get('pythia_version')} "
                  f"{refj.get('extra_settings')}, this run {version:.3f} {args.set}")
        else:
            mc_ref = np.array(refj["species"]["pi"]["mc"])
            dev = np.max(np.abs(np.array(r["mc_total"]) / np.where(mc_ref > 0, mc_ref, np.nan) - 1)[k])
            out["closure_vs_identified_max_rel_dev_pi"] = float(dev)
            print(f"[closure] sum over origins vs {ref.name}: max |rel. dev.| {dev:.2e} (pi, 0.1-3 GeV/c)")

    # --- print
    if not out["origins_valid"]:
        print("[warn] Bose-Einstein / rescattering on: origin split NOT valid, only totals are")
    print(f"[pi] mean pT by origin: " + ", ".join(f"{o} {v:.3f}" for o, v in r["mean_pt"].items() if v))
    print(f"[pi] share of all |y|<0.5 pions: " + ", ".join(f"{o} {100 * v:.1f}%" for o, v in r["yield_share"].items()))
    for sp in SPECIES:
        for band, b in out["species"][sp]["bands"].items():
            top = sorted(b["fraction"].items(), key=lambda x: -x[1])[:4]
            print(f"[{sp:2s}] {band:>9s} GeV/c: MC/data {b['mc_over_data']:.3f}; " +
                  ", ".join(f"{o} {100 * f:.0f}%" for o, f in top))
    dg = out["pion_rescale_diagnostic"]
    print(f"[pi] rescale diagnostic 0.1-3 GeV/c: direct x{dg['scale_direct']:.3f}, decay x{dg['scale_decay']:.3f}; "
          f"chi2 {dg['chi2_unscaled']:.0f} -> {dg['chi2_rescaled']:.0f} / {dg['n_bins']} bins")

    out["label"] = args.label
    plot(out, figdir / f"m02_pion_origin{args.tag}")
    path = tabdir / f"m02_pion_origin{args.tag}.json"
    path.write_text(json.dumps(out, indent=1))
    print(f"[out] {path}, {figdir / f'm02_pion_origin{args.tag}.pdf'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
