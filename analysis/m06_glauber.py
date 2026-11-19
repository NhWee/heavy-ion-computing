#!/usr/bin/env python3
"""Module 6: Monte Carlo Glauber for Pb-Pb and p-Pb at sqrt(s_NN) = 5.02 TeV.

    python analysis/m06_glauber.py [--n-events 40000] [--n-syst 10000] [--workers 4]

Outputs (results/tables/m06/m06_glauber_summary.json, results/figures/m06/):
  - hadronic cross sections sigma_PbPb and sigma_pPb
  - <N_part>, <N_coll>, <b>, T_AA = <N_coll> / sigma_NN in centrality classes
  - validation against ALICE-PUBLIC-2018-011 (Pb-Pb 5.02 TeV, Table 1):
    0-1 %: <N_part> 401.9, <N_coll> 1949;  4-5 %: 362.9, 1593
  - systematic variations as ALICE does them: R +- 0.06 fm, a +- 0.010 fm,
    sigma_NN +- 0.6 mb, d_min 0 and 0.8 fm
  - figures: sampled nuclear density vs Woods-Saxon, N_part and N_coll vs b,
    one event display

Why it matters for what came before: N_coll converts a pp cross section into
the expectation for Pb-Pb (R_AA = yield_AA / (<N_coll> x yield_pp)); N_part
is the scale of bulk particle production; Module 4's background rho grows
with N_part.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from hicpy.glauber import PB208, PROTON, run, centrality_table, sample_nucleus, collide   # noqa: E402
from hicpy.style import apply_style, save_figure, OKABE_ITO, plt                         # noqa: E402

SIGMA_NN = 67.6           # mb, 5.02 TeV
ALICE_REF = {"0-1%": (401.9, 1949.0), "4-5%": (362.9, 1593.0)}
EDGES = (0, 1, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)


def _run_chunk(args):
    n, a, b, sig, seed = args
    return run(n, a, b, sig, seed=seed)


def run_parallel(n, a, b, sig, workers, seed0):
    chunks = [(n // workers + (i < n % workers), a, b, sig, seed0 + i) for i in range(workers)]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        parts = list(ex.map(_run_chunk, chunks))
    out = {k: np.concatenate([p[k] for p in parts]) for k in ("b", "npart", "ncoll")}
    hit = sum(len(p["b"]) for p in parts)
    frac = hit / n
    out["sigma_inel_b"] = float(np.pi * 20.0 ** 2 * frac / 100)
    out["sigma_inel_err_b"] = float(np.pi * 20.0 ** 2 * np.sqrt(frac * (1 - frac) / n) / 100)
    out["n_events_thrown"] = n
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-events", type=int, default=40000)
    ap.add_argument("--n-syst", type=int, default=10000)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--figdir", default=str(REPO / "results" / "figures" / "m06"))
    ap.add_argument("--tabdir", default=str(REPO / "results" / "tables" / "m06"))
    args = ap.parse_args()
    figdir, tabdir = pathlib.Path(args.figdir), pathlib.Path(args.tabdir)
    figdir.mkdir(parents=True, exist_ok=True); tabdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------ nominal
    pbpb = run_parallel(args.n_events, PB208, PB208, SIGMA_NN, args.workers, 100)
    table = centrality_table(pbpb, EDGES, order="b")
    table_np = centrality_table(pbpb, EDGES, order="npart")
    for r in table:
        r["T_AA_per_mb"] = r["ncoll"] / SIGMA_NN
    print(f"[Pb-Pb] {args.n_events} thrown, sigma_PbPb = {pbpb['sigma_inel_b']:.3f} ± "
          f"{pbpb['sigma_inel_err_b']:.3f} b (stat)")
    print(f"  {'class':>7s} {'<Npart>':>8s} {'<Ncoll>':>8s} {'<b> fm':>7s} {'T_AA 1/mb':>9s}   (ordered by b)")
    for r in table:
        print(f"  {r['class']:>7s} {r['npart']:8.1f} {r['ncoll']:8.0f} {r['b']:7.2f} {r['T_AA_per_mb']:9.2f}")
    ppb = run_parallel(max(args.n_events // 2, 1000), PROTON, PB208, SIGMA_NN, args.workers, 200)
    print(f"[p-Pb]  sigma_pPb = {ppb['sigma_inel_b']:.3f} ± {ppb['sigma_inel_err_b']:.3f} b; "
          f"<N_coll> = {ppb['ncoll'].mean():.2f}, <N_part> = {ppb['npart'].mean():.2f}")

    # ------------------------------------------------------------ validation
    val = {}
    for cls, (np_ref, nc_ref) in ALICE_REF.items():
        rb = next(r for r in table if r["class"] == cls)
        rn = next(r for r in table_np if r["class"] == cls)
        val[cls] = {"alice_npart": np_ref, "alice_ncoll": nc_ref,
                    "npart_b_order": rb["npart"], "ncoll_b_order": rb["ncoll"],
                    "npart_npart_order": rn["npart"], "ncoll_npart_order": rn["ncoll"]}
        print(f"[val] {cls}: <Npart> {rb['npart']:.1f} (b-order) / {rn['npart']:.1f} (Npart-order) vs ALICE {np_ref}; "
              f"<Ncoll> {rb['ncoll']:.0f} / {rn['ncoll']:.0f} vs {nc_ref:.0f}")

    # ------------------------------------------------------------ systematics
    variations = {
        "R+": (replace(PB208, R=6.68), SIGMA_NN), "R-": (replace(PB208, R=6.56), SIGMA_NN),
        "a+": (replace(PB208, a=0.556), SIGMA_NN), "a-": (replace(PB208, a=0.536), SIGMA_NN),
        "sigma+": (PB208, SIGMA_NN + 0.6), "sigma-": (PB208, SIGMA_NN - 0.6),
        "dmin0": (replace(PB208, d_min=0.0), SIGMA_NN), "dmin0.8": (replace(PB208, d_min=0.8), SIGMA_NN),
        "stat (independent seeds)": (PB208, SIGMA_NN),
    }
    # COMMON RANDOM NUMBERS: the reference and every variation use the SAME seeds,
    # so they share the impact parameters (drawn first) and nearly the same nucleon
    # draws.  With independent seeds the b-percentile class boundaries of a 10k
    # sample fluctuate by ~0.1 fm, i.e. 3-13 % in <Npart> of peripheral classes --
    # larger than the physics variations (troubleshooting #23).  "stat" below is an
    # independent-seed nominal run: the noise floor the comparison would have had.
    ref = run_parallel(args.n_syst, PB208, PB208, SIGMA_NN, args.workers, 300)
    ref_tab = {r["class"]: r for r in centrality_table(ref, EDGES)}
    syst = {}
    for k, (nuc, sig) in variations.items():
        v = run_parallel(args.n_syst, nuc, nuc, sig, args.workers, 900 if k.startswith("stat") else 300)
        vt = {r["class"]: r for r in centrality_table(v, EDGES)}
        syst[k] = {"sigma_inel_b": v["sigma_inel_b"],
                   "rel": {c: {"npart": vt[c]["npart"] / ref_tab[c]["npart"] - 1,
                               "ncoll": vt[c]["ncoll"] / ref_tab[c]["ncoll"] - 1}
                           for c in ("0-1%", "5-10%", "30-40%", "70-80%")}}
    phys = {k: v for k, v in syst.items() if not k.startswith("stat")}
    for k, v in syst.items():
        print(f"[syst] {k:26s} sigma {v['sigma_inel_b']:.3f} b  " + "  ".join(
            f"{c}: {100 * r['npart']:+.1f}/{100 * r['ncoll']:+.1f} %" for c, r in v["rel"].items()))
    sig_spread = max(abs(s["sigma_inel_b"] - ref["sigma_inel_b"]) for s in phys.values())
    print(f"[syst] sigma_PbPb spread over variations: ±{sig_spread:.3f} b (reference {ref['sigma_inel_b']:.3f} b, "
          f"stat ±{ref['sigma_inel_err_b']:.3f})")
    for c in ("0-1%", "5-10%", "30-40%", "70-80%"):
        dn = max(abs(s["rel"][c]["npart"]) for s in phys.values())
        dc = max(abs(s["rel"][c]["ncoll"]) for s in phys.values())
        print(f"[syst] {c}: max |d<Npart>| {100 * dn:.1f} %, max |d<Ncoll>| {100 * dc:.1f} %")

    # ------------------------------------------------------------ figures
    apply_style()
    rng = np.random.default_rng(7)
    rs = np.concatenate([np.linalg.norm(sample_nucleus(rng, PB208, recentre=False), axis=1) for _ in range(300)])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    h, e = np.histogram(rs, bins=60, range=(0, 12))
    c = 0.5 * (e[1:] + e[:-1]); shell = 4 * np.pi * c ** 2 * np.diff(e)
    dens = h / shell / 300
    ax.errorbar(c, dens, yerr=np.sqrt(h) / shell / 300, fmt="o", ms=4, color=OKABE_ITO[0], label="sampled nucleons")
    ws = 1 / (1 + np.exp((c - PB208.R) / PB208.a)); rho0 = 208 / np.sum(4 * np.pi * c ** 2 * ws * np.diff(e))
    ax.plot(c, rho0 * ws, "-", color=OKABE_ITO[1], lw=2, label=rf"Woods–Saxon, R = {PB208.R} fm, a = {PB208.a} fm")
    ax.set_xlabel("r (fm)"); ax.set_ylabel(r"$\rho(r)$ (nucleons / fm$^3$)"); ax.legend(fontsize=10)
    save_figure(fig, figdir / "m06_density"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    sel = rng.choice(len(pbpb["b"]), size=min(4000, len(pbpb["b"])), replace=False)
    ax.plot(pbpb["b"][sel], pbpb["npart"][sel], ".", ms=2, color=OKABE_ITO[0], alpha=0.5, label=r"$N_{\rm part}$")
    ax.plot(pbpb["b"][sel], pbpb["ncoll"][sel] / 5, ".", ms=2, color=OKABE_ITO[1], alpha=0.5, label=r"$N_{\rm coll}/5$")
    ax.set_xlabel("impact parameter b (fm)"); ax.set_ylabel("count"); ax.legend(markerscale=6)
    ax.set_title(rf"Pb–Pb, $\sigma_{{NN}}$ = {SIGMA_NN} mb", fontsize=11)
    save_figure(fig, figdir / "m06_npart_ncoll_vs_b"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    pa, pb_ = sample_nucleus(rng, PB208), sample_nucleus(rng, PB208)
    bb = 7.0
    _, _, wa, wb = collide(pa, pb_, bb, SIGMA_NN)
    for pos, w, col, sh in ((pa, wa, OKABE_ITO[0], +bb / 2), (pb_, wb, OKABE_ITO[1], -bb / 2)):
        x, y = pos[:, 0] + sh, pos[:, 1]
        ax.scatter(x[~w], y[~w], s=40, facecolors="none", edgecolors=col, lw=0.8)
        ax.scatter(x[w], y[w], s=40, color=col, alpha=0.8)
    ax.set_aspect("equal"); ax.set_xlim(-14, 14); ax.set_ylim(-14, 14)
    ax.set_xlabel("x (fm)"); ax.set_ylabel("y (fm)")
    ax.set_title(f"b = {bb} fm: N_part = {int(wa.sum() + wb.sum())} (filled)", fontsize=11)
    save_figure(fig, figdir / "m06_event_display"); plt.close(fig)

    summary = {"inputs": {"R_fm": PB208.R, "a_fm": PB208.a, "d_min_fm": PB208.d_min, "sigma_NN_mb": SIGMA_NN,
                          "sources": ["arXiv:1301.4361 (R, a, d_min)", "arXiv:1710.07098 (sigma_NN at 5.02 TeV)"]},
               "pbpb": {"sigma_inel_b": pbpb["sigma_inel_b"], "sigma_inel_err_b": pbpb["sigma_inel_err_b"],
                        "n_events_thrown": args.n_events, "centrality_b_order": table,
                        "centrality_npart_order": table_np},
               "ppb": {"sigma_inel_b": ppb["sigma_inel_b"], "sigma_inel_err_b": ppb["sigma_inel_err_b"],
                       "ncoll_mean": float(ppb["ncoll"].mean()), "npart_mean": float(ppb["npart"].mean())},
               "validation_alice_public_2018_011": val,
               "systematics": {"reference_sigma_b": ref["sigma_inel_b"], "variations": syst}}
    (tabdir / "m06_glauber_summary.json").write_text(json.dumps(summary, indent=1))
    print(f"[out] {figdir}/m06_*.pdf|png, {tabdir / 'm06_glauber_summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
