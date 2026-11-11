#!/usr/bin/env python3
"""Module 1, step 4: do the three arms give the SAME numbers?

    python analysis/m01_compare_arms.py [--allow-missing]

Three independent implementations read the same file:

    uproot arm   analysis/m01_analysis_uproot.py  -> results/tables/m01_fit_results.json
    PyROOT arm   analysis/m01_analysis_pyroot.py  -> results/tables/m01_pyroot_summary.json
    C++ macro    analysis/m01_analysis_root.C     -> results/tables/m01_root_summary.json

The uproot arm is the reference; BOTH ROOT arms must be present (a missing
summary is a failure, not a pass -- use --allow-missing to compare only what
exists, and the verdict then says so).  Before any number is compared, two
preconditions are checked, because a comparison of the wrong inputs proves
nothing:

  0a. same input    every arm records the provenance of the file it read
                    (produced_utc, git commit, generator config); they must be
                    identical -- a stale summary from an earlier run fails here
  0b. same cuts     the track and muon cut values each arm actually applied

Then three levels are checked, because agreeing on a few scalar summaries is
NOT evidence that two implementations agree:

  1. scalars        <N_ch>, dN/deta, dimuon count and mean mass
  2. histograms     bin edges, then every bin's content and error, for the
                    four histograms (P(N_ch), dN/deta, invariant pT yield,
                    m(mumu))
  3. Tsallis fit    T and n, compared in units of their statistical error.
                    This is the least-squares fit at bin centres that all
                    three arms implement identically (TF1 chi2 / curve_fit) --
                    a check of the implementations.  The physics result of the
                    uproot arm is a Poisson-likelihood fit ("tsallis_poisson"),
                    which the ROOT arms do not run; its inputs (the selected
                    particles and the binning) are covered by levels 1-2.

When the arms disagree, in this order of likelihood the cause is:

  1. bin edges          (linear vs log, or a different number of bins; the C++
                         macro duplicates the config numbers -- check them)
  2. cut set            (one arm applies pt < pt_max, the other does not)
  3. normalisation      (ROOT's Scale(1/N,"width") vs dividing by hand;
                         per-particle 1/pT weight vs the bin-centre 1/pT)
  4. cut ordering       (cutting on a *derived* quantity computed from clipped
                         values)
  5. float32 vs float64 arithmetic at a cut boundary: a particle with
                         |eta| = 0.7999999 lands on different sides of the cut
                         in the two arms.  Shows up as a ONE-COUNT difference
                         in a bin, never as a percent-level shift.

Only the fifth is allowed to survive.  The per-bin tolerance below is designed
to pass exactly that case and nothing else: a bin agrees if the relative
difference is < 1e-6, OR if the absolute difference is smaller than half the
statistical error of the bin (one migrating particle in a bin of N entries is
1/sqrt(N) sigma).  The number of bins that needed the second clause is printed,
so a systematic effect hiding behind it is still visible.

The fit tolerances have NOT yet been tuned against a real ROOT run (see
docs/next_steps.md); if MINUIT and scipy converge to slightly different minima,
loosen FIT_TOL_SIGMA and record why in docs/troubleshooting.md.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
TAB = REPO / "results" / "tables"

REFERENCE_FILE = "m01_fit_results.json"
ARM_FILES = [("pyroot", "m01_pyroot_summary.json"),
             ("root_macro", "m01_root_summary.json")]

# (label, key in uproot summary, key in ROOT summaries, relative tolerance)
SCALARS = [
    ("n_events",            "n_events",                "n_events",             0.0),
    ("<N_ch> (|eta|<0.8)",  "mean_Nch_in_acceptance",  "mean_Nch",             1e-6),
    ("dN/deta (|eta|<0.5)", "dNdeta_central",          "dNdeta_central",       1e-6),
    ("dimuon candidates",   "n_dimuon_candidates",     "n_dimuon_candidates",  0.0),
    ("<m(mumu)> (GeV)",     "mean_dimuon_mass_GeV",    "mean_dimuon_mass_GeV", 1e-5),
]
HISTS = ["mult", "eta", "pt", "mass"]
BIN_REL_TOL = 1e-6        # clause 1: float round-off
BIN_SIGMA_TOL = 0.5       # clause 2: one migrating particle at a cut boundary
EDGE_REL_TOL = 1e-9
FIT_TOL_SIGMA = 0.05      # |T_a - T_b| < 0.05 sigma_T  (MINUIT vs scipy noise)
FIT_REL_TOL = 1e-3


def rel_diff(a, b):
    return 0.0 if a == b else abs(a - b) / max(abs(a), abs(b), 1e-300)


def compare_scalars(up, rt):
    ok = True
    print(f"  {'quantity':22s} {'uproot':>16s} {'ROOT':>16s} {'rel. diff':>12s}")
    for label, k_up, k_rt, tol in SCALARS:
        a, b = up.get(k_up), rt.get(k_rt)
        if a is None or b is None:
            print(f"  {label:22s} {'--':>16s} {'--':>16s} {'MISSING KEY':>12s}  FAIL")
            ok = False
            continue
        r = rel_diff(a, b)
        good = r <= tol
        ok &= good
        print(f"  {label:22s} {a:16.6f} {b:16.6f} {r:12.2e}  "
              f"{'ok' if good else 'MISMATCH'}")
    return ok


def compare_hist(name, ha, hb):
    """Return (ok, message) for one histogram."""
    if ha is None or hb is None:
        return False, "missing in one arm"
    ea, eb = ha["edges"], hb["edges"]
    if len(ea) != len(eb):
        return False, f"different number of bins ({len(ea)-1} vs {len(eb)-1})"
    worst_edge = max(rel_diff(x, y) for x, y in zip(ea, eb))
    if worst_edge > EDGE_REL_TOL:
        return False, f"bin edges differ (max rel. diff {worst_edge:.1e})"

    n_bad, n_clause2, worst = 0, 0, 0.0
    for i, (ca, cb, sa, sb) in enumerate(zip(ha["contents"], hb["contents"],
                                             ha["errors"], hb["errors"])):
        r = rel_diff(ca, cb)
        worst = max(worst, r)
        if r <= BIN_REL_TOL and rel_diff(sa, sb) <= max(BIN_REL_TOL, 1e-5):
            continue
        sigma = max(sa, sb)
        if sigma > 0 and abs(ca - cb) <= BIN_SIGMA_TOL * sigma:
            n_clause2 += 1
            continue
        n_bad += 1
    nb = len(ha["contents"])
    msg = (f"{nb} bins, max rel. diff {worst:.1e}"
           + (f", {n_clause2} bin(s) within {BIN_SIGMA_TOL} sigma "
              f"(boundary migration)" if n_clause2 else ""))
    if n_bad:
        msg += f", {n_bad} bin(s) MISMATCH"
    if n_clause2 > max(3, nb // 10):
        msg += "  <- too many boundary bins for round-off alone; investigate"
        return False, msg
    return n_bad == 0, msg


def compare_fit(fa, fb):
    if fa is None or fb is None:
        return False, "Tsallis fit missing in one arm"
    ok = True
    parts = []
    for p in ("T", "n"):
        a, b = fa[p], fb[p]
        sig = max(fa[f"{p}_err"], fb[f"{p}_err"], 1e-300)
        d_sig = abs(a - b) / sig
        r = rel_diff(a, b)
        good = d_sig <= FIT_TOL_SIGMA or r <= FIT_REL_TOL
        ok &= good
        parts.append(f"{p}: {a:.6g} vs {b:.6g} ({d_sig:.3f} sigma, "
                     f"{'ok' if good else 'MISMATCH'})")
    chi = f"chi2/ndf {fa['chi2']/max(fa['ndf'],1):.3f} vs {fb['chi2']/max(fb['ndf'],1):.3f}"
    return ok, "; ".join(parts) + "; " + chi


PROVENANCE_KEYS = ("produced_utc", "git_commit", "config")
CUT_KEYS = {"track": ("pt_min", "pt_max", "eta_abs_max"), "muon": ("pt_min", "eta_abs_max")}


def compare_inputs(up, rt):
    """0a: both arms read the same file; 0b: both applied the same cut values."""
    ok = True
    pu, pr = up.get("generator_provenance"), rt.get("input_provenance")
    if not pu or not pr:
        print("  [0a] input provenance MISSING in one summary -- re-run that arm "
              "(summaries written before 2026-11-07 do not record it)")
        ok = False
    else:
        diff = [k for k in PROVENANCE_KEYS if pu.get(k) != pr.get(k)]
        if diff:
            print(f"  [0a] DIFFERENT INPUT FILES: {', '.join(diff)} differ "
                  f"({pu.get('produced_utc')} vs {pr.get('produced_utc')}) -- a stale summary?")
            ok = False
        else:
            print(f"  [0a] same input file (produced {pu['produced_utc']}, "
                  f"git {pu['git_commit']}, seed {pu['config'].get('seed')})")
    cu, cr = up.get("cuts"), rt.get("cuts")
    if not cu or not cr:
        print("  [0b] cut values MISSING in one summary -- re-run that arm")
        return False
    bad = [f"{blk}.{k}: {cu[blk][k]} vs {cr[blk].get(k)}"
           for blk, keys in CUT_KEYS.items() for k in keys
           if cr.get(blk, {}).get(k) is None or rel_diff(cu[blk][k], cr[blk][k]) > 1e-12]
    if bad:
        print("  [0b] DIFFERENT CUTS: " + "; ".join(bad))
        ok = False
    else:
        print("  [0b] same cut values")
    return ok


def compare_arm(name, up, rt):
    print(f"\n=== uproot arm vs {name} ({rt.get('root_version', '?')}) ===")
    print("[0] inputs")
    ok = compare_inputs(up, rt)

    print("[1] scalars")
    ok &= compare_scalars(up, rt)

    print("[2] histograms (bin by bin)")
    hu, hr = up.get("histograms", {}), rt.get("histograms", {})
    for h in HISTS:
        good, msg = compare_hist(h, hu.get(h), hr.get(h))
        ok &= good
        print(f"  {h:6s} {'ok      ' if good else 'MISMATCH'} {msg}")

    print("[3] Tsallis fit")
    good, msg = compare_fit(up.get("tsallis"), rt.get("tsallis"))
    ok &= good
    print(f"  {'ok      ' if good else 'MISMATCH'} {msg}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--allow-missing", action="store_true",
                    help="compare only the ROOT arms whose summary exists")
    ap.add_argument("--tabdir", default=str(TAB), help="directory with the summaries")
    args = ap.parse_args()
    tab = pathlib.Path(args.tabdir)
    REFERENCE = ("uproot", tab / REFERENCE_FILE)
    ARMS = [(n, tab / f) for n, f in ARM_FILES]
    if not REFERENCE[1].exists():
        print(f"missing {REFERENCE[1].name}; run analysis/m01_analysis_uproot.py first")
        return 1
    up = json.loads(REFERENCE[1].read_text())
    if "histograms" not in up:
        print(f"{REFERENCE[1].name} has no 'histograms' block: re-run "
              "analysis/m01_analysis_uproot.py (schema changed)")
        return 1

    found = [(n, p) for n, p in ARMS if p.exists()]
    if not found:
        print("no ROOT-arm summary found; expected one of:")
        for n, p in ARMS:
            print(f"  {p}   ({n})")
        print("(those arms need ROOT -- they are skipped automatically without it)")
        return 1

    ok = True
    for name, path in found:
        ok &= compare_arm(name, up, json.loads(path.read_text()))
    missing = [n for n, p in ARMS if not p.exists()]
    if missing:
        print(f"\nNOT CHECKED: {', '.join(missing)} (no summary file)")
        if not args.allow_missing:
            ok = False

    if ok and missing:
        print("\nARMS AGREE for the arms present -- NOT a full three-arm check")
    else:
        print("\nARMS " + ("AGREE" if ok else "DISAGREE -- see the checklist in this "
                                             "file's docstring"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
