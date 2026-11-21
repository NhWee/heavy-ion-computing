#!/usr/bin/env python3
"""Module 8: the CMS dimuon spectrum from CERN Open Data (pp 8 TeV, 2012).

    bash scripts/m08_fetch_cms_dimuon.sh                        # WSL2: download + provenance (2 GB)
    python analysis/m08_dimuon.py histogram [--max-events N]    # WSL2: ROOT file -> histogram JSON
    python analysis/m08_dimuon.py fit [--tag _x]                # anywhere: histogram JSON -> fits, figures
    python analysis/m08_dimuon.py toy                           # closure test input (no data needed)

DATA
----
CERN Open Data record 12341, "DoubleMuParked dataset from 2012 in NanoAOD
format reduced on muons" (DOI 10.7483/OPENDATA.CMS.LVG5.QT81, CC0): ~61 M
events of CMS Run 2012 B and C, muon variables only (nMuon, Muon_pt, eta,
phi, mass, charge).  Analysis-code record: 12342.  The ROOT file is NOT in
git (2 GB); scripts/m08_fetch_cms_dimuon.sh downloads it and writes the
provenance JSON with its sha256.  The histogram JSON (small) IS committed,
so the fit step runs anywhere.

WHY THIS, FOR HEAVY IONS
------------------------
Quarkonia (J/psi, Upsilon) and the Z are the standard "thermometer" and
"standard candle" of heavy-ion physics: the loosely bound Upsilon(2S, 3S)
melt before the tightly bound Upsilon(1S) in a quark-gluon plasma (sequential
suppression), while the Z does not interact strongly and must scale with
N_coll (Module 6).  Both claims are made relative to a pp reference --
this module builds that reference from real data.

SELECTION (as in the CMS outreach analysis, record 12342)
  exactly two muons, opposite charge; mass from (pT, eta, phi) with the muon
  mass.  A second set requires both muons in the barrel, |eta| < 1.2, where
  the momentum resolution is best and the three Upsilon states separate.

OUTPUT
  results/tables/m08/m08_hist{tag}.json      histograms (committed)
  results/tables/m08/m08_fit{tag}.json       masses, resolutions, yields, ratios
  results/figures/m08/m08_spectrum{tag}.png, m08_fits{tag}.png
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys
import time

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from hicpy.dimuon import PDG, pair_mass, fit_window, toy_events     # noqa: E402

TAB, FIG = REPO / "results" / "tables" / "m08", REPO / "results" / "figures" / "m08"
RAW = REPO / "data" / "raw" / "opendata" / "Run2012BC_DoubleMuParked_Muons.root"
PROV = REPO / "data" / "raw" / "cms_2012_doublemuparked.json"
TOY = REPO / "data" / "generated" / "m08_toy_dimuon.root"
BARREL = 1.2
SPECTRUM_EDGES = np.geomspace(0.25, 300, 601)
WINDOWS = {"jpsi": (2.6, 4.2, 0.01, ["J/psi", "psi(2S)"], "cb"),
           "upsilon": (8.5, 11.5, 0.02, ["Upsilon(1S)", "Upsilon(2S)", "Upsilon(3S)"], "cb"),
           "z": (70.0, 110.0, 0.25, ["Z"], "voigt")}


def window_edges(w):
    lo, hi, step, _, _ = WINDOWS[w]
    return np.round(np.arange(lo, hi + step / 2, step), 6)


# ----------------------------------------------------------------------------- histogram
def histogram(args) -> int:
    import awkward as ak
    import uproot
    path = pathlib.Path(args.input)
    if not path.exists():
        print(f"{path} not found -- run scripts/m08_fetch_cms_dimuon.sh (WSL2) or `toy` for a test file")
        return 1
    if path.resolve() == RAW.resolve():
        prov = json.loads(PROV.read_text()) if PROV.exists() else {}
        if not prov.get("verified") or prov.get("size_bytes") != path.stat().st_size:
            print(f"{path.name}: no verified provenance for this file (incomplete download?) -- "
                  "run scripts/m08_fetch_cms_dimuon.sh first")
            return 1
    hists = {sel: {"spectrum": np.zeros(len(SPECTRUM_EDGES) - 1),
                   **{w: np.zeros(len(window_edges(w)) - 1) for w in WINDOWS}} for sel in ("all", "barrel")}
    n_read = n_pairs = 0
    t0 = time.time()
    branches = ["nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge"]
    with uproot.open(path) as f:
        tree = f["Events"]
        n_total = tree.num_entries
        stop = min(n_total, args.max_events) if args.max_events else n_total
        for arr in tree.iterate(branches, step_size=args.step, entry_stop=stop, library="ak"):
            n = ak.to_numpy(arr["nMuon"]); n_read += len(n)
            two = n == 2
            if not two.any():
                continue
            # events with exactly two muons -> regular (N, 2) numpy arrays
            reg = lambda b: ak.to_numpy(ak.to_regular(arr[b][two], axis=1)).astype(float)  # noqa: E731
            pt, eta, phi, q = reg("Muon_pt"), reg("Muon_eta"), reg("Muon_phi"), reg("Muon_charge")
            os_ = q[:, 0] != q[:, 1]
            m = pair_mass(pt[os_, 0], eta[os_, 0], phi[os_, 0], pt[os_, 1], eta[os_, 1], phi[os_, 1])
            barrel = (np.abs(eta[os_, 0]) < BARREL) & (np.abs(eta[os_, 1]) < BARREL)
            n_pairs += len(m)
            for sel, mm in (("all", m), ("barrel", m[barrel])):
                hists[sel]["spectrum"] += np.histogram(mm, SPECTRUM_EDGES)[0]
                for w in WINDOWS:
                    hists[sel][w] += np.histogram(mm, window_edges(w))[0]
            print(f"\r[read] {n_read:,} / {stop:,} events, {n_pairs:,} opposite-sign pairs, "
                  f"{time.time() - t0:.0f} s", end="", flush=True)
    print()
    prov = json.loads(PROV.read_text()) if (PROV.exists() and path.resolve() == RAW.resolve()) else None
    out = {"input": path.name, "provenance": prov, "n_events_read": int(n_read), "n_events_in_file": int(n_total),
           "n_opposite_sign_pairs": int(n_pairs),
           "selection": "nMuon == 2, opposite charge; 'barrel': both |eta| < 1.2",
           "spectrum_edges": SPECTRUM_EDGES.tolist(),
           "windows": {w: {"edges": window_edges(w).tolist(), "states": WINDOWS[w][3], "model": WINDOWS[w][4]}
                       for w in WINDOWS},
           "hist": {sel: {k: v.astype(int).tolist() for k, v in h.items()} for sel, h in hists.items()}}
    TAB.mkdir(parents=True, exist_ok=True)
    p = TAB / f"m08_hist{args.tag}.json"
    p.write_text(json.dumps(out))
    print(f"[out] {p}  ({n_pairs:,} pairs from {n_read:,} events)")
    return 0


# ----------------------------------------------------------------------------- fit
def fit(args) -> int:
    from hicpy.style import apply_style, save_figure, OKABE_ITO, plt
    h = json.loads((TAB / f"m08_hist{args.tag}.json").read_text())
    res = {"input": h["input"], "n_events_read": h["n_events_read"], "provenance": h.get("provenance"),
           "fits": {}}
    for sel in ("all", "barrel"):
        res["fits"][sel] = {}
        for w, spec in h["windows"].items():
            e = np.array(spec["edges"]); c = np.array(h["hist"][sel][w], float)
            if c.sum() < 50:
                continue
            r = fit_window(e, c, spec["states"], spec["model"])
            if not r["converged"]:
                r = fit_window(e, c, spec["states"], spec["model"], start=r["params"])
            res["fits"][sel][w] = r
    # summary numbers
    summ = {}
    for sel, fits in res["fits"].items():
        s = {}
        for w, r in fits.items():
            p, er = r["params"], r["errors"]
            ref = r["states"][0]
            for st in r["states"]:
                s[st] = {"mass": r[f"m_{st}"], "mass_err": r[f"m_{st}_err"],
                         "shift_MeV": 1000 * (r[f"m_{st}"] - PDG[st]), "mass_scale": p["scale"],
                         "mass_scale_err": er["scale"],
                         "sigma": p["sigma"] * PDG[st] / PDG[ref],
                         "sigma_over_m": p["sigma"] / PDG[ref],
                         "yield": p[f"N_{st}"], "yield_err": er[f"N_{st}"], "S/B_2sigma": r[f"S/B_{st}"]}
            s[f"chi2/ndf {w}"] = r["chi2_poisson"] / max(r["ndf"], 1)
        for num, den in (("Upsilon(2S)", "Upsilon(1S)"), ("Upsilon(3S)", "Upsilon(1S)"), ("psi(2S)", "J/psi")):
            if num in s and den in s:
                rr = s[num]["yield"] / s[den]["yield"]
                s[f"{num}/{den}"] = {"ratio": rr, "err": rr * np.hypot(s[num]["yield_err"] / s[num]["yield"],
                                                                       s[den]["yield_err"] / s[den]["yield"])}
        summ[sel] = s
    res["summary"] = summ
    FIG.mkdir(parents=True, exist_ok=True)
    out = TAB / f"m08_fit{args.tag}.json"
    out.write_text(json.dumps(res, indent=1, default=float))
    for sel, s in summ.items():
        print(f"== {sel}")
        for st in PDG:
            if st in s:
                v = s[st]
                print(f"  {st:12s} m = {v['mass']:.4f} +- {v['mass_err']:.4f} GeV (PDG {PDG[st]:.4f}, scale "
                      f"{v['mass_scale']:.5f})  sigma/m {100 * v['sigma_over_m']:.2f} %  "
                      f"N = {v['yield']:.0f} +- {v['yield_err']:.0f}  S/B "
                      + (f"{v['S/B_2sigma']:.2f}" if v['S/B_2sigma'] < 1e3 else "> 1000"))
        for k, v in s.items():
            if isinstance(v, dict) and "ratio" in v:
                print(f"  {k:24s} = {v['ratio']:.4f} +- {v['err']:.4f}  (raw yields: no acceptance correction)")
            elif k.startswith("chi2"):
                print(f"  {k} = {v:.2f}")

    # figures
    apply_style()
    e = np.array(h["spectrum_edges"]); c = 0.5 * (e[1:] + e[:-1])
    fig, ax = plt.subplots(figsize=(10, 5.6))
    for sel, col in (("all", OKABE_ITO[0]), ("barrel", OKABE_ITO[1])):
        y = np.array(h["hist"][sel]["spectrum"]) / np.diff(e)
        ax.step(c, y, where="mid", lw=1.2, color=col, label={"all": "opposite-sign pairs",
                                                             "barrel": r"both $|\eta| < 1.2$"}[sel])
    labels = {r"$\eta$": 0.548, r"$\rho,\omega$": 0.78, r"$\phi$": 1.019, r"J/$\psi$": 3.097, r"$\psi$(2S)": 3.686,
              r"$\Upsilon$(1S,2S,3S)": 9.8, "Z": 91.19}
    ytop = np.array(h["hist"]["all"]["spectrum"]) / np.diff(e)
    for t, m in labels.items():
        k = np.searchsorted(c, m)
        if 0 < k < len(c) and ytop[k] > 0:
            ax.annotate(t, (m, ytop[max(0, k - 3):k + 4].max() * 1.6), ha="center", fontsize=10)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"$m_{\mu\mu}$ (GeV)")
    ax.set_ylabel("pairs / GeV"); ax.legend(frameon=False)
    ax.set_title(f"CMS Open Data 2012 (record 12341), {h['n_events_read'] / 1e6:.1f} M events", fontsize=11)
    save_figure(fig, FIG / f"m08_spectrum{args.tag}"); plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(14, 6.6), gridspec_kw={"height_ratios": [3, 1]}, sharex="col")
    for j, w in enumerate(WINDOWS):
        sel = "barrel" if w == "upsilon" else "all"
        r = res["fits"].get(sel, {}).get(w)
        if r is None:
            continue
        ed = np.array(h["windows"][w]["edges"]); cc = 0.5 * (ed[1:] + ed[:-1])
        n = np.array(h["hist"][sel][w], float); mu = np.array(r["expected"])
        ax, axp = axes[0, j], axes[1, j]
        ax.errorbar(cc, n, yerr=np.sqrt(n), fmt="o", ms=2, color="k", lw=0.8, label=f"data ({sel})")
        ax.plot(cc, mu, color=OKABE_ITO[1], lw=1.6, label="fit")
        ax.set_title({"jpsi": r"J/$\psi$, $\psi$(2S)", "upsilon": r"$\Upsilon$(1S, 2S, 3S), barrel", "z": "Z"}[w])
        ax.legend(fontsize=8, frameon=False)
        axp.plot(cc, (n - mu) / np.sqrt(np.clip(mu, 1, None)), ".", color="k", ms=2)
        axp.axhline(0, color="0.5", lw=0.8); axp.set_ylim(-5, 5); axp.set_xlabel(r"$m_{\mu\mu}$ (GeV)")
    axes[0, 0].set_ylabel(r"pairs / bin"); axes[1, 0].set_ylabel("pull")
    fig.tight_layout(); save_figure(fig, FIG / f"m08_fits{args.tag}"); plt.close(fig)
    print(f"[out] {out}, {FIG / f'm08_spectrum{args.tag}.pdf'}, {FIG / f'm08_fits{args.tag}.pdf'}")
    return 0


# ----------------------------------------------------------------------------- toy
TOY_COMPOSITION = {"J/psi": 0.08, "psi(2S)": 0.003, "Upsilon(1S)": 0.02, "Upsilon(2S)": 0.006,
                   "Upsilon(3S)": 0.003, "Z": 0.01}


def toy(args) -> int:
    import awkward as ak
    import uproot
    rng = np.random.default_rng(args.seed)
    t = toy_events(rng, args.n_events, TOY_COMPOSITION, sigma_rel=args.sigma_rel)
    counts = t["nMuon"]
    def jag(x, dtype):                                             # noqa: E306
        return ak.unflatten(np.asarray(x, dtype), counts)
    TOY.parent.mkdir(parents=True, exist_ok=True)
    muon = ak.zip({"pt": jag(t["pt"], np.float32), "eta": jag(t["eta"], np.float32),
                   "phi": jag(t["phi"], np.float32), "mass": jag(t["mass"], np.float32),
                   "charge": jag(t["charge"], np.int32)})
    with uproot.recreate(TOY) as f:
        # NanoAOD naming: Muon_pt, ..., counter nMuon (uproot's defaults differ between versions)
        f.mktree("Events", {"Muon": muon.type.content}, counter_name=lambda c: "n" + c,
                 field_name=lambda outer, inner: inner if outer == "" else outer + "_" + inner)
        f["Events"].extend({"Muon": muon})
    truth = {s: int(round(TOY_COMPOSITION[s] * args.n_events)) for s in TOY_COMPOSITION}
    (TOY.with_suffix(".json")).write_text(json.dumps({"n_events": args.n_events, "seed": args.seed,
                                                      "sigma_rel": args.sigma_rel, "composition": TOY_COMPOSITION,
                                                      "expected_counts": truth}, indent=1))
    print(f"[toy] {TOY} ({args.n_events} events), truth -> {TOY.with_suffix('.json')}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("histogram")
    h.add_argument("--input", default=str(RAW)); h.add_argument("--max-events", type=int, default=None)
    h.add_argument("--step", default="100 MB"); h.add_argument("--tag", default="")
    f = sub.add_parser("fit"); f.add_argument("--tag", default="")
    t = sub.add_parser("toy"); t.add_argument("--n-events", type=int, default=400000)
    t.add_argument("--seed", type=int, default=8); t.add_argument("--sigma-rel", type=float, default=0.01)
    a = ap.parse_args()
    return {"histogram": histogram, "fit": fit, "toy": toy}[a.cmd](a)


if __name__ == "__main__":
    raise SystemExit(main())
