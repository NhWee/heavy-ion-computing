#!/usr/bin/env python3
"""Module 3, step 2: read the HepMC3 file and check what an event record promises.

    python analysis/m03_event_record.py [--in data/generated/m03_pythia_mb_13TeV.hepmc]
                                        [--figdir results/figures/m03] [--tabdir results/tables/m03]

Checks (each one says something about the physics bookkeeping, not the code):

  1. global      sum of final-state 4-momenta = sum of the two beams, and total
                 final charge = +2 (two protons), event by event
  2. counting    HepMC status 1 particles = PYTHIA isFinal particles; no
                 non-final particle is left without an end vertex
  3. per vertex  4-momentum balance at every vertex, broken down by the stage of
                 the history that produced the outgoing particles (decays,
                 hadronization, showers, MPI, ...)
  4. ancestry    one final-state particle followed back to the beam

Expected outcome of check 3 -- and the reason it is the interesting one: decays
and string fragmentation conserve 4-momentum vertex by vertex, but parton
showers do not, because PYTHIA's dipole showers take the recoil of a branching
from a SECOND parton that is recorded as a new copy with its own history.  The
energy flows between vertices; only the whole event balances.  "Every vertex
conserves momentum" is therefore NOT a valid check of a PYTHIA record.
"""
from __future__ import annotations

import argparse
import collections
import json
import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from hicpy.eventrecord import STAGES, vertex_balance, global_checks, ancestry, _int_attr  # noqa: E402
from hicpy.style import apply_style, save_figure, OKABE_ITO, plt                        # noqa: E402

PDG_NAMES = {2212: "p", -2212: "p̄", 211: "π+", -211: "π−", 111: "π0", 321: "K+", -321: "K−",
             113: "ρ0", 213: "ρ+", -213: "ρ−", 223: "ω", 221: "η", 331: "η′", 313: "K*0",
             323: "K*+", -323: "K*−", 2224: "Δ++", 2114: "Δ0", 21: "g", 1: "d", 2: "u",
             3: "s", -1: "d̄", -2: "ū", -3: "s̄", 22: "γ", 2101: "(ud)0", 2103: "(ud)1",
             2203: "(uu)1", 1103: "(dd)1", 92: "string", 93: "cluster", 9902210: "p_diffr"}


def name(pid):
    return PDG_NAMES.get(pid, str(pid))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="infile",
                    default=str(REPO / "data" / "generated" / "m03_pythia_mb_13TeV.hepmc"))
    ap.add_argument("--figdir", default=str(REPO / "results" / "figures" / "m03"))
    ap.add_argument("--tabdir", default=str(REPO / "results" / "tables" / "m03"))
    ap.add_argument("--rel-tol", type=float, default=1e-6)
    args = ap.parse_args()
    import pyhepmc

    figdir, tabdir = pathlib.Path(args.figdir), pathlib.Path(args.tabdir)
    figdir.mkdir(parents=True, exist_ok=True); tabdir.mkdir(parents=True, exist_ok=True)

    agg = collections.defaultdict(lambda: {"n": 0, "conserved": 0, "worst": 0.0})
    glob = {"n_events": 0, "max_p4_rel_diff": 0.0, "charge_ok": 0, "count_ok": 0,
            "dangling": 0, "mother_conflicts": 0, "n_vertices": 0, "n_particles": 0}
    example = None
    with pyhepmc.open(args.infile) as f:
        for ev in f:                       # the event object stays alive in the loop body
            glob["n_events"] += 1
            glob["n_vertices"] += len(ev.vertices)
            glob["n_particles"] += len(ev.particles)
            g = global_checks(ev)
            glob["max_p4_rel_diff"] = max(glob["max_p4_rel_diff"], g["p4_rel_diff"])
            glob["charge_ok"] += int(abs(g["charge_final"] - 2.0) < 1e-9)
            n_final_pythia = sum(1 for p in ev.particles if _int_attr(p, "pythia_status") > 0)
            glob["count_ok"] += int(n_final_pythia == g["n_final"])
            glob["dangling"] += g["n_no_end_vertex_nonbeam"]
            if "n_mother_conflicts" in ev.attributes:
                glob["mother_conflicts"] += _int_attr(ev, "n_mother_conflicts")
            for k, r in vertex_balance(ev, args.rel_tol).items():
                a = agg[k]; a["n"] += r["n"]; a["conserved"] += r["conserved"]
                a["worst"] = max(a["worst"], r["worst"])
            if example is None:
                # a final charged pion from a resonance decay: a short, readable chain
                for p in ev.particles:
                    if p.status == 1 and abs(p.pid) == 211 and _int_attr(p, "pythia_status") == 91:
                        example = {"event": ev.event_number,
                                   "chain": [(name(pid), st) for pid, st in ancestry(p)]}
                        break

    print(f"[read] {glob['n_events']} events, {glob['n_particles']} particles, "
          f"{glob['n_vertices']} vertices")
    print(f"[1] global: max |Σp(final) − Σp(beams)| / E = {glob['max_p4_rel_diff']:.1e}; "
          f"charge = +2 in {glob['charge_ok']}/{glob['n_events']} events")
    print(f"[2] counting: status-1 = PYTHIA final in {glob['count_ok']}/{glob['n_events']}; "
          f"non-final without end vertex: {glob['dangling']}; mother conflicts: "
          f"{glob['mother_conflicts']}")
    print(f"[3] per-vertex 4-momentum balance (|Δp| < {args.rel_tol:g} E_in), "
          f"by stage of the outgoing particles:")
    rows = []
    for k in sorted(agg, key=lambda s: [int(x) for x in s.split("+")]):
        a = agg[k]
        label = " + ".join(STAGES.get(int(x), x) for x in k.split("+"))
        frac = a["conserved"] / a["n"]
        rows.append({"stage": k, "label": label, **a, "fraction": frac})
        print(f"    {label:32s} {a['n']:7d} vertices  {100 * frac:6.1f} % conserved  "
              f"worst {a['worst']:.1e}")
    if example:
        print("[4] ancestry of a final π± from a decay (event "
              f"{example['event']}): " + " ← ".join(f"{n}({s})" for n, s in example["chain"]))

    # figure: fraction of conserving vertices for the single-stage rows
    single = [r for r in rows if "+" not in r["stage"]]
    apply_style()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    y = np.arange(len(single))
    cols = [OKABE_ITO[0] if r["fraction"] > 0.99 else OKABE_ITO[1] for r in single]
    ax.barh(y, [100 * r["fraction"] for r in single], color=cols)
    for yi, r in zip(y, single):
        ax.text(101, yi, f"{r['n']} vertices", va="center", fontsize=10, color="0.35")
    ax.set_yticks(y); ax.set_yticklabels([r["label"] for r in single]); ax.invert_yaxis()
    ax.set_xlim(0, 125); ax.set_xlabel("vertices conserving 4-momentum (%)")
    ax.tick_params(axis="y", which="both", right=False, left=False)
    save_figure(fig, figdir / "m03_vertex_balance")
    plt.close(fig)

    summary = {"input": args.infile, "rel_tol": args.rel_tol, "global": glob,
               "per_stage": rows, "example_ancestry": example}
    (tabdir / "m03_event_record_summary.json").write_text(json.dumps(summary, indent=2,
                                                                     ensure_ascii=False))
    print(f"[out] {figdir / 'm03_vertex_balance.pdf'}, {tabdir / 'm03_event_record_summary.json'}")
    ok = (glob["n_events"] > 0 and glob["max_p4_rel_diff"] < 1e-6
          and glob["charge_ok"] == glob["n_events"]
          and glob["count_ok"] == glob["n_events"])
    print("GLOBAL CHECKS " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
