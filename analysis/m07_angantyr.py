#!/usr/bin/env python3
"""Module 7: PYTHIA 8 Angantyr for p-Pb and Pb-Pb at sqrt(s_NN) = 5.02 TeV.

    python analysis/m07_angantyr.py generate --system PbPb [--n-events 5000] [--workers N] [--tag _x]
    python analysis/m07_angantyr.py generate --system pPb  [--n-events 200000]
    python analysis/m07_angantyr.py generate --system pp   [--n-events 200000]     # 5.02 TeV reference
    python analysis/m07_angantyr.py compare

QUESTIONS
---------
1. Angantyr contains its own Glauber model (fluctuating nucleon sizes and
   cross sections, fitted to PYTHIA's pp cross sections).  How do its
   N_part, N_coll and sigma_inel compare with Module 6's black-disk
   Glauber and with the measured sigma_pPb and sigma_PbPb?
2. With NO collective effects, how much charged-particle multiplicity does
   it make: (2/<N_part>) dN_ch/deta versus <N_part> for pp, p-Pb and Pb-Pb
   -- the "wounded-nucleon scaling" question.
3. Strangeness: Module 2 showed PYTHIA's pp strange/pi ratios are flat in
   multiplicity while data rise.  Does stacking nucleon-nucleon collisions
   (Angantyr, still no ropes) change that, up to dN_ch/deta ~ 2000?

FRAMES AND SELECTIONS (as ALICE)
--------------------------------
- p-Pb: ALICE quotes eta_lab with the PROTON moving to negative eta and the
  nucleon-nucleon centre of mass moving with rapidity -0.465 in the lab.
  PYTHIA's beam A (the proton) moves along +z in the NN frame, hence
  eta_lab = -eta_NN - 0.465 (exact for rapidity, the eta shift is the usual
  approximation).  Pb-Pb and pp: lab = NN frame.
- V0A: 2.8 < eta_lab < 5.1, V0C: -3.7 < eta_lab < -1.7 (charged particles).
  Pb-Pb centrality: V0M = V0A + V0C, all events (Angantyr generates
  inelastic collisions).  p-Pb: "NSD" = at least one charged particle in
  BOTH V0A and V0C (ALICE's trigger), centrality from V0A (Pb-going side).
- Identified/strange counts in |y_lab| < 0.5; K0S, Lambda, Xi, Omega are not
  decayed (tau0Max = 10 mm), as in Module 2.

DATA USED (checked values; full centrality tables need HEPData ins1410589)
  Pb-Pb 0-5 %: dN_ch/deta (|eta| < 0.5) = 1943 +- 54  (ALICE, arXiv:1512.06104)
  p-Pb NSD:    dN_ch/deta (|eta_lab| < 0.5) = 16.81 +- 0.71, 2.14 +- 0.17 per participant
               (ALICE, arXiv:1210.3615)
  sigma_PbPb = 7.67 +- 0.25 b (ALICE, arXiv:2204.10148); sigma_pPb = 2.061 +- 0.080 b
  (CMS, arXiv:1509.03893); Glauber reference: Module 6.
  Pending: the full Pb-Pb centrality dependence (HEPData ins1410589 Table 1),
  to be downloaded by hand -- hepdata.net is not reachable from the cloud session.
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

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
GEN, TAB, FIG = (REPO / "data" / "generated", REPO / "results" / "tables" / "m07",
                 REPO / "results" / "figures" / "m07")
CARDS = {"pPb": REPO / "generators" / "angantyr_pPb_5TeV.cmnd",
         "PbPb": REPO / "generators" / "angantyr_PbPb_5TeV.cmnd",
         "pp": REPO / "generators" / "pythia_mb_13TeV.cmnd"}
Y_SHIFT = {"pPb": -0.465, "PbPb": 0.0, "pp": 0.0}
FLIP = {"pPb": True, "PbPb": False, "pp": False}
SPECIES = {"pi": 211, "K": 321, "p": 2212, "K0S": 310, "Lambda": 3122, "Xi": 3312, "Omega": 3334}
PT_SPECIES = ("pi", "K", "p")
ETA_EDGES = np.linspace(-5, 5, 21)
SEED0 = 20261121
RATIOS = {"K/pi": ("K", 1.0), "p/pi": ("p", 1.0), "2K0S/pi": ("K0S", 2.0), "Lambda/pi": ("Lambda", 1.0),
          "Xi/pi": ("Xi", 1.0), "Omega/pi": ("Omega", 1.0)}
COLS = (["b", "npart_proj", "npart_targ", "ncoll", "ncoll_nd", "nv0a", "nv0c", "nch05", "nch1"]
        + [f"eta{i}" for i in range(len(ETA_EDGES) - 1)]
        + [f"n_{s}" for s in SPECIES] + [f"sumpt_{s}" for s in PT_SPECIES])

DATA = {"PbPb_0-5_dndeta": (1943.0, 54.0, "ALICE, arXiv:1512.06104, |eta| < 0.5"),
        "pPb_NSD_dndeta": (16.81, 0.71, "ALICE, arXiv:1210.3615, |eta_lab| < 0.5, NSD"),
        "pPb_NSD_per_participant": (2.14, 0.17, "ALICE, arXiv:1210.3615"),
        "sigma_PbPb_b": (7.67, 0.25, "ALICE, arXiv:2204.10148"),
        "sigma_pPb_b": (2.061, 0.080, "CMS, arXiv:1509.03893 (stat+syst+lumi)"),
        "ALICE_PbPb_glauber": {"0-1%": (401.9, 1949.0), "4-5%": (362.9, 1593.0),
                               "source": "ALICE-PUBLIC-2018-011 (TGlauberMC)"}}


# ----------------------------------------------------------------------------- generation
def _chunk(job):
    system, n_events, seed, settings = job
    try:
        import pythia8
    except ImportError:
        import pythia8mc as pythia8
    py = pythia8.Pythia("", False)
    py.readFile(str(CARDS[system]))
    base = ["Random:setSeed = on", f"Random:seed = {seed}", "Print:quiet = on", "Next:numberCount = 0"]
    if system == "pp":
        base += ["Beams:eCM = 5020."]
    for s in base + list(settings):
        if not py.readString(s):
            raise RuntimeError(f"PYTHIA rejected {s!r}")
    if not py.init():
        raise RuntimeError("PYTHIA init() failed")
    hi = py.getHeavyIonsPtr().hiInfo if system != "pp" else None
    flip, shift = FLIP[system], Y_SHIFT[system]
    code = {v: k for k, v in SPECIES.items()}
    sidx = {s: COLS.index(f"n_{s}") for s in SPECIES}
    pidx = {s: COLS.index(f"sumpt_{s}") for s in PT_SPECIES}
    e0 = COLS.index("eta0"); neta = len(ETA_EDGES) - 1
    rows = []
    n_acc = n_tried = 0
    while n_acc < n_events and n_tried < 2 * n_events + 100:
        n_tried += 1
        if not py.next():
            continue
        n_acc += 1
        row = np.zeros(len(COLS))
        if hi is not None:
            row[0:5] = hi.b(), hi.nPartProj(), hi.nPartTarg(), hi.nCollTot(), hi.nCollND()
        else:
            row[0:5] = np.nan, 1, 1, 1, 1
        ev = py.event
        for i in range(ev.size()):
            p = ev[i]
            if not p.isFinal():
                continue
            if p.isCharged():
                eta = p.eta()
                eta = (-eta if flip else eta) + shift
                if 2.8 < eta < 5.1:
                    row[5] += 1
                elif -3.7 < eta < -1.7:
                    row[6] += 1
                if abs(eta) < 1.0:
                    row[8] += 1
                    if abs(eta) < 0.5:
                        row[7] += 1
                k = int((eta + 5.0) * 2.0)
                if 0 <= k < neta and -5.0 <= eta < 5.0:
                    row[e0 + k] += 1
            sp = code.get(abs(p.id()))
            if sp is not None:
                y = p.y()
                y = (-y if flip else y) + shift
                if abs(y) < 0.5:
                    row[sidx[sp]] += 1
                    if sp in pidx:
                        row[pidx[sp]] += p.pT()
        rows.append(row)
    sig = (hi.glauberINEL(), hi.glauberINELErr()) if hi is not None else (np.nan, np.nan)
    return np.array(rows), sig, py.settings.parm("Pythia:versionNumber")


def percentile_classes(key, edges, rng):
    """Index arrays of percentile classes, largest key first (random tie break)."""
    k = key + rng.random(len(key))
    order = np.argsort(-k, kind="stable"); n = len(order)
    return [(f"{lo:g}-{hi:g}%", order[int(round(n * lo / 100)):int(round(n * hi / 100))])
            for lo, hi in zip(edges[:-1], edges[1:])]


def class_record(name, d, col):
    npart = d[:, col["npart_proj"]] + d[:, col["npart_targ"]]
    rec = {"class": name, "n_events": int(len(d)), "b": float(np.nanmean(d[:, col["b"]])) if len(d) else None,
           "npart": float(npart.mean()), "ncoll": float(d[:, col["ncoll"]].mean()),
           "dndeta": float(d[:, col["nch05"]].mean()), "dndeta_err": float(d[:, col["nch05"]].std() / np.sqrt(len(d))),
           "ratios": {}, "ratio_stat": {}, "meanpt": {}, "counts": {}}
    rec["dndeta_per_part_pair"] = rec["dndeta"] / (rec["npart"] / 2)
    npi = d[:, col["n_pi"]].sum()
    for r, (sp, fac) in RATIOS.items():
        nx = d[:, col[f"n_{sp}"]].sum()
        rec["ratios"][r] = float(fac * nx / npi) if npi else None
        rec["ratio_stat"][r] = float(fac * np.sqrt(nx) / npi) if npi else None
        rec["counts"][sp] = int(nx)
    rec["counts"]["pi"] = int(npi)
    for sp in PT_SPECIES:
        n = d[:, col[f"n_{sp}"]].sum()
        rec["meanpt"][sp] = float(d[:, col[f"sumpt_{sp}"]].sum() / n) if n else None
    return rec


def summarise(system, data, sigmas, version, settings):
    col = {c: i for i, c in enumerate(COLS)}
    rng = np.random.default_rng(5)
    out = {"system": system, "pythia_version": round(version, 3), "extra_settings": settings,
           "n_events": int(len(data)), "frame": "lab: eta_lab = -eta_NN - 0.465" if system == "pPb" else "NN = lab"}
    s = np.array([x for x in sigmas if np.isfinite(x[0]) and x[1] > 0])
    if len(s):
        w = 1 / s[:, 1] ** 2
        out["sigma_inel_mb"] = float(np.sum(w * s[:, 0]) / w.sum()); out["sigma_inel_err_mb"] = float(1 / np.sqrt(w.sum()))
    if system == "PbPb":
        sel = data
        out["selection"] = "all generated (inelastic) events; centrality = V0A + V0C"
        key = sel[:, col["nv0a"]] + sel[:, col["nv0c"]]
        edges = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        glauber_edges = [0, 1, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    elif system == "pPb":
        nsd = (data[:, col["nv0a"]] > 0) & (data[:, col["nv0c"]] > 0)
        sel = data[nsd]
        out["selection"] = "NSD: >= 1 charged in V0A and in V0C; centrality = V0A (Pb side)"
        out["nsd_fraction"] = float(nsd.mean())
        key = sel[:, col["nv0a"]]
        edges = [0, 5, 10, 20, 40, 60, 80, 100]; glauber_edges = None
    else:
        sel = data[data[:, col["nch1"]] > 0]
        out["selection"] = "INEL>0 (>= 1 charged in |eta| < 1); classes = V0A + V0C as in Module 2"
        key = sel[:, col["nv0a"]] + sel[:, col["nv0c"]]
        edges = [0, 0.92, 4.6, 9.2, 13.8, 18.4, 27.6, 36.8, 46, 64.5, 100]; glauber_edges = None
    out["all"] = class_record("all", sel, col)
    out["dndeta_vs_eta"] = {"edges": ETA_EDGES.tolist(),
                            "all": (sel[:, col["eta0"]:col["eta0"] + 20].mean(axis=0) / 0.5).tolist()}
    out["classes"] = [class_record(n, sel[i], col) for n, i in percentile_classes(key, edges, rng)]
    for c, (n, i) in zip(out["classes"], percentile_classes(key, edges, np.random.default_rng(5))):
        out["dndeta_vs_eta"][n] = (sel[i][:, col["eta0"]:col["eta0"] + 20].mean(axis=0) / 0.5).tolist()
    if glauber_edges:
        b = sel[:, col["b"]]
        out["glauber_b_order"] = [class_record(n, sel[i], col) for n, i in percentile_classes(-b, glauber_edges, rng)]
        out["glauber_v0m_order"] = [class_record(n, sel[i], col) for n, i in percentile_classes(key, glauber_edges, rng)]
    return out


def generate(args) -> int:
    w = args.workers or os.cpu_count() or 1
    w = min(w, args.n_events)
    jobs = [(args.system, args.n_events // w + (i < args.n_events % w), SEED0 + 1000 * i, args.set) for i in range(w)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=w) as ex:
        parts = list(ex.map(_chunk, jobs))
    data = np.vstack([p[0] for p in parts]); version = parts[0][2]
    GEN.mkdir(parents=True, exist_ok=True); TAB.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(GEN / f"m07_{args.system}{args.tag}.npz", data=data, columns=np.array(COLS))
    out = summarise(args.system, data, [p[1] for p in parts], version, args.set)
    out["label"] = args.label or {"pPb": "Angantyr p-Pb", "PbPb": "Angantyr Pb-Pb", "pp": "PYTHIA pp 5.02 TeV"}[args.system]
    out["seeds"] = [j[2] for j in jobs]
    path = TAB / f"m07_{args.system}{args.tag}.json"
    path.write_text(json.dumps(out, indent=1))
    a = out["all"]
    print(f"[gen] {args.system}: {len(data)} events in {time.time() - t0:.0f} s on {w} workers, PYTHIA {version:.3f}")
    if "sigma_inel_mb" in out:
        print(f"[sigma] Angantyr sigma_inel = {out['sigma_inel_mb'] / 1000:.3f} +- {out['sigma_inel_err_mb'] / 1000:.3f} b")
    print(f"[all]  {out['selection']}: {a['n_events']} events, <Npart> {a['npart']:.2f}, <Ncoll> {a['ncoll']:.2f}, "
          f"<dNch/deta> {a['dndeta']:.2f}")
    for c in out["classes"]:
        print(f"  {c['class']:>8s}: <b> {c['b'] if c['b'] is not None else float('nan'):5.2f}  <Npart> {c['npart']:6.1f}  "
              f"<Ncoll> {c['ncoll']:7.1f}  <dNch/deta> {c['dndeta']:8.1f}  per pair {c['dndeta_per_part_pair']:5.2f}")
    print(f"[out] {path}")
    return 0


# ----------------------------------------------------------------------------- compare
def compare(args) -> int:
    sys.path.insert(0, str(REPO / "analysis"))
    from hicpy.style import apply_style, save_figure, OKABE_ITO, plt
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    runs = {s: json.loads((TAB / f"m07_{s}.json").read_text()) for s in ("pp", "pPb", "PbPb")
            if (TAB / f"m07_{s}.json").exists()}
    g6 = json.loads((REPO / "results" / "tables" / "m06" / "m06_glauber_summary.json").read_text())
    FIG.mkdir(parents=True, exist_ok=True)
    summary = {"data": {k: v for k, v in DATA.items()}, "systems": {}}
    apply_style()

    # 1. Glauber: Angantyr vs Module 6
    if "PbPb" in runs:
        pb = runs["PbPb"]
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
        for ax, q, lab in zip(axes, ("npart", "ncoll"), (r"$\langle N_{\rm part}\rangle$", r"$\langle N_{\rm coll}\rangle$")):
            def mid(rows):
                return np.array([np.mean([float(x) for x in r["class"].rstrip("%").split("-")]) for r in rows])
            for rows, style, name in ((g6["pbpb"]["centrality_b_order"], dict(color=OKABE_ITO[0], ls="-"), "Module 6, b order"),
                                      (g6["pbpb"]["centrality_npart_order"], dict(color=OKABE_ITO[0], ls=":"), "Module 6, N$_{part}$ order"),
                                      (pb["glauber_b_order"], dict(color=OKABE_ITO[1], ls="-"), "Angantyr, b order"),
                                      (pb["glauber_v0m_order"], dict(color=OKABE_ITO[1], ls="--"), "Angantyr, V0M order")):
                ax.plot(mid(rows), [r[q] for r in rows], marker="o", ms=3, label=name, **style)
            ref = DATA["ALICE_PbPb_glauber"]
            ax.plot([0.5, 4.5], [ref["0-1%"][0 if q == "npart" else 1], ref["4-5%"][0 if q == "npart" else 1]], "k*", ms=10,
                    label="ALICE-PUBLIC-2018-011")
            ax.set_yscale("log"); ax.set_xlabel("centrality (%)"); ax.set_ylabel(lab)
        axes[0].legend(fontsize=8, frameon=False)
        fig.tight_layout(); save_figure(fig, FIG / "m07_glauber_vs_m06"); plt.close(fig)
        summary["glauber"] = {
            "sigma_PbPb_b": {"angantyr": pb.get("sigma_inel_mb", np.nan) / 1000, "module6": g6["pbpb"]["sigma_inel_b"],
                             "data": DATA["sigma_PbPb_b"][:2]},
            "classes": [{"class": a["class"], "angantyr_b": [a["npart"], a["ncoll"], a["b"]],
                         "angantyr_v0m": [v["npart"], v["ncoll"]],
                         "module6_b": [m["npart"], m["ncoll"], m["b"]]}
                        for a, v, m in zip(pb["glauber_b_order"], pb["glauber_v0m_order"], g6["pbpb"]["centrality_b_order"])]}
    if "pPb" in runs:
        summary["sigma_pPb_b"] = {"angantyr": runs["pPb"].get("sigma_inel_mb", np.nan) / 1000,
                                  "module6": g6["ppb"]["sigma_inel_b"], "data": DATA["sigma_pPb_b"][:2]}

    # 2. multiplicity per participant pair vs Npart
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    colors = {"pp": OKABE_ITO[2], "pPb": OKABE_ITO[3], "PbPb": OKABE_ITO[1]}
    for s, r in runs.items():
        cls = r["classes"] if s != "pp" else [r["all"]]
        ax.plot([c["npart"] for c in cls], [c["dndeta_per_part_pair"] for c in cls], "s-", color=colors[s], ms=4,
                label=r["label"] + (" (all INEL>0)" if s == "pp" else ""))
        if s == "pPb":
            a = r["all"]; ax.plot([a["npart"]], [a["dndeta_per_part_pair"]], "D", mfc="none", color=colors[s], ms=8,
                                  label="Angantyr p-Pb NSD (all)")
    v, e, _ = DATA["pPb_NSD_per_participant"]
    ax.errorbar([16.81 / 2.14], [2 * v], yerr=[2 * e], fmt="D", color="k", ms=6, label="ALICE p-Pb NSD")
    np05 = g6["pbpb"]["centrality_npart_order"]
    npart05 = sum(r["npart"] * r["n"] for r in np05[:3]) / sum(r["n"] for r in np05[:3])
    v, e, _ = DATA["PbPb_0-5_dndeta"]
    ax.errorbar([npart05], [v / (npart05 / 2)], yerr=[e / (npart05 / 2)], fmt="o", color="k", ms=6,
                label=r"ALICE Pb-Pb 0-5% ($N_{\rm part}$ from Module 6)")
    # The full Pb-Pb centrality dependence (HEPData ins1410589 Table 1) is added
    # once that table is in data/raw/ and its column layout has been checked.
    ax.set_xscale("log"); ax.set_xlabel(r"$\langle N_{\rm part}\rangle$")
    ax.set_ylabel(r"$(2/\langle N_{\rm part}\rangle)\,\langle dN_{\rm ch}/d\eta\rangle$")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout(); save_figure(fig, FIG / "m07_dndeta_per_participant"); plt.close(fig)

    # 3. strangeness vs multiplicity, pp (data + MC) -> p-Pb -> Pb-Pb (Angantyr)
    sys.path.insert(0, str(REPO / "analysis"))
    import m02_multiplicity as mm
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.6))
    pp13 = json.loads((REPO / "results" / "tables" / "m02" / "m02_mult.json").read_text())
    for ax, (r, (_, _, tab, dep)) in zip(axes.flat, mm.RATIOS.items()):
        d = mm.read_table(tab, dep)
        ax.errorbar(d["x"], d["y"], yerr=np.hypot(d["stat"], d["sys"]), fmt="o", color="k", ms=4, label="ALICE pp 13 TeV")
        ax.plot([c["dndeta"] for c in pp13["classes"]], [c["ratios"][r] for c in pp13["classes"]], "-", color=OKABE_ITO[0],
                lw=1.5, label="PYTHIA pp 13 TeV (Monash)")
        for s, run in runs.items():
            ax.plot([c["dndeta"] for c in run["classes"]], [c["ratios"][r] for c in run["classes"]], "s-", ms=3,
                    color=colors[s], lw=1.5, label=run["label"])
        ax.set_xscale("log"); ax.set_title(r, fontsize=12); ax.set_xlabel(r"$\langle dN_{\rm ch}/d\eta\rangle$")
        ticks = [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
        ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FixedFormatter([str(t) for t in ticks])); ax.set_xlim(2, 3000)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.legend(h, l, loc="lower center", ncol=5, fontsize=9, frameon=False)
    save_figure(fig, FIG / "m07_strangeness_vs_multiplicity"); plt.close(fig)

    # 4. dN/deta(eta)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, s in zip(axes, ("pPb", "PbPb")):
        if s not in runs:
            continue
        r = runs[s]; e = np.array(r["dndeta_vs_eta"]["edges"]); c = 0.5 * (e[1:] + e[:-1])
        for name in (["all"] + [k["class"] for k in r["classes"]][:: (2 if s == "PbPb" else 1)]):
            ax.step(c, r["dndeta_vs_eta"][name], where="mid", lw=1.4, label=name)
        if s == "pPb":
            v, err, _ = DATA["pPb_NSD_dndeta"]; ax.errorbar([0], [v], yerr=[err], fmt="D", color="k", label="ALICE NSD")
        else:
            ax.set_yscale("log"); v, err, _ = DATA["PbPb_0-5_dndeta"]
            ax.errorbar([0], [v], yerr=[err], fmt="o", color="k", label="ALICE 0-5%")
        ax.set_xlabel(r"$\eta_{\rm lab}$" + (" (Pb going to +)" if s == "pPb" else "")); ax.set_ylabel(r"$dN_{\rm ch}/d\eta$")
        ax.set_title(r["label"]); ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout(); save_figure(fig, FIG / "m07_dndeta_vs_eta"); plt.close(fig)

    # printed summary
    for s, r in runs.items():
        a = r["all"]
        summary["systems"][s] = {"label": r["label"], "pythia_version": r["pythia_version"], "n_events": r["n_events"],
                                 "selection": r["selection"], "all": {k: a[k] for k in ("npart", "ncoll", "dndeta",
                                                                                         "dndeta_per_part_pair")},
                                 "sigma_inel_mb": r.get("sigma_inel_mb"),
                                 "classes": [{k: c[k] for k in ("class", "npart", "ncoll", "dndeta", "dndeta_per_part_pair",
                                                                 "ratios", "meanpt")} for c in r["classes"]]}
        print(f"[{s}] {r['label']}: {a['n_events']} events ({r['selection'][:40]}...), <Npart> {a['npart']:.2f}, "
              f"<dNch/deta> {a['dndeta']:.2f}" + (f", sigma {r['sigma_inel_mb'] / 1000:.3f} b" if r.get('sigma_inel_mb') else ""))
    if "pPb" in runs:
        a = runs["pPb"]["all"]
        print(f"[pPb] NSD dNch/deta: MC {a['dndeta']:.2f} vs ALICE 16.81 +- 0.71 -> {a['dndeta'] / 16.81:.3f}; "
              f"per participant {a['dndeta'] / a['npart']:.2f} vs 2.14 +- 0.17")
    if "PbPb" in runs:
        c05 = runs["PbPb"]["classes"][0]
        print(f"[PbPb] 0-5% dNch/deta: MC {c05['dndeta']:.0f} vs ALICE 1943 +- 54 -> {c05['dndeta'] / 1943:.3f}; "
              f"<Npart> {c05['npart']:.1f} (Module 6: {npart05:.1f})")
        gl = summary["glauber"]
        print(f"[sigma] Pb-Pb: Angantyr {gl['sigma_PbPb_b']['angantyr']:.3f} b, Module 6 {gl['sigma_PbPb_b']['module6']:.3f} b, "
              f"ALICE 7.67 +- 0.25 b")
        for cl in gl["classes"]:
            print(f"   {cl['class']:>8s}: Npart Angantyr(b) {cl['angantyr_b'][0]:6.1f}  (V0M) {cl['angantyr_v0m'][0]:6.1f}  "
                  f"Module6 {cl['module6_b'][0]:6.1f} | Ncoll {cl['angantyr_b'][1]:7.1f} {cl['angantyr_v0m'][1]:7.1f} "
                  f"{cl['module6_b'][1]:7.1f} | b {cl['angantyr_b'][2]:5.2f} vs {cl['module6_b'][2]:5.2f}")
    if "sigma_pPb_b" in summary:
        sp = summary["sigma_pPb_b"]
        print(f"[sigma] p-Pb: Angantyr {sp['angantyr']:.3f} b, Module 6 {sp['module6']:.3f} b, CMS 2.061 +- 0.080 b")
    (TAB / "m07_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    print(f"[out] {TAB / 'm07_summary.json'}, {FIG}/m07_*.pdf")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate")
    g.add_argument("--system", choices=["pp", "pPb", "PbPb"], required=True)
    g.add_argument("--n-events", type=int, default=None)
    g.add_argument("--workers", type=int, default=None)
    g.add_argument("--set", action="append", default=[])
    g.add_argument("--tag", default="")
    g.add_argument("--label", default=None)
    sub.add_parser("compare")
    args = ap.parse_args()
    if args.cmd == "generate":
        args.n_events = args.n_events or {"pp": 200000, "pPb": 200000, "PbPb": 5000}[args.system]
        return generate(args)
    return compare(args)


if __name__ == "__main__":
    raise SystemExit(main())
