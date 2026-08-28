#!/usr/bin/env python3
"""Module 1, step 4: do the uproot arm and the ROOT arm give the SAME numbers?

    python analysis/m01_compare_arms.py

Two independent implementations reading the same file must agree.  When they do
not, in this order of likelihood the cause is:

  1. bin edges          (linear vs log, or a different number of bins)
  2. normalisation      (ROOT's Scale(1/N,"width") vs dividing by hand;
                         per-particle 1/pT weight vs the bin-centre 1/pT)
  3. cut ordering       (cutting on eta before or after a pT cut changes nothing;
                         cutting on a *derived* quantity computed from clipped
                         values does)
  4. float32 vs float64 arithmetic at a cut boundary

Only the fourth is allowed to survive.  Everything else is a bug in one arm.

Inputs
  results/tables/m01_fit_results.json     written by m01_analysis_uproot.py
  results/tables/m01_pyroot_summary.json  written by m01_analysis_pyroot.py
"""
from __future__ import annotations

import json
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
TAB = REPO / "results" / "tables"

# (label, key in uproot summary, key in pyroot summary, relative tolerance)
CHECKS = [
    ("n_events",            "n_events",                "n_events",             0.0),
    ("<N_ch> (|eta|<0.8)",  "mean_Nch_in_acceptance",  "mean_Nch",             1e-4),
    ("dN/deta (|eta|<0.5)", "dNdeta_central",          "dNdeta_central",       1e-4),
    ("dimuon candidates",   "n_dimuon_candidates",     "n_dimuon_candidates",  0.0),
    ("<m(mumu)> (GeV)",     "mean_dimuon_mass_GeV",    "mean_dimuon_mass_GeV", 1e-4),
]


def main() -> int:
    f_up = TAB / "m01_fit_results.json"
    f_rt = TAB / "m01_pyroot_summary.json"
    if not f_up.exists():
        print(f"missing {f_up.name}; run analysis/m01_analysis_uproot.py first")
        return 1
    if not f_rt.exists():
        print(f"missing {f_rt.name}; run analysis/m01_analysis_pyroot.py first\n"
              "(that arm needs ROOT -- it is skipped automatically without it)")
        return 1

    up = json.loads(f_up.read_text())
    rt = json.loads(f_rt.read_text())

    print(f"uproot arm vs ROOT arm ({rt.get('root_version', '?')})\n")
    print(f"{'quantity':22s} {'uproot':>16s} {'ROOT':>16s} {'rel. diff':>12s}  ")
    ok = True
    for label, k_up, k_rt, tol in CHECKS:
        a, b = up.get(k_up), rt.get(k_rt)
        if a is None or b is None:
            print(f"{label:22s} {'--':>16s} {'--':>16s} {'MISSING KEY':>12s}  FAIL")
            ok = False
            continue
        rel = 0.0 if a == b else abs(a - b) / max(abs(a), 1e-30)
        good = rel <= tol
        ok &= good
        print(f"{label:22s} {a:16.6f} {b:16.6f} {rel:12.2e}  "
              f"{'ok' if good else 'MISMATCH'}")

    print("\nARMS " + ("AGREE" if ok else "DISAGREE -- see the checklist in this "
                                         "file's docstring"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
