#!/usr/bin/env python3
"""Do the regenerated results match the committed ones?

    python scripts/check_results_drift.py [--ref HEAD]

Run after scripts/run_all.sh.  Same seed + same code must give the same
numbers, so every results/tables/*.json that is committed at --ref is compared
with the file on disk:

* provenance blocks are skipped -- they record the machine, the time and the
  commit, which legitimately differ (keys in SKIP_KEYS);
* floats must agree to 1e-6 relative (1e-9 absolute near zero): enough slack
  for round-off between machines, far tighter than any physics change;
* everything else (counts, strings, structure) must be identical.

Figures are not compared byte by byte here: PNG rasterisation can change
between matplotlib releases.  `git diff --stat -- results/figures` shows them.
Exit code 1 lists every path that drifted.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess

REPO = pathlib.Path(__file__).resolve().parent.parent
SKIP_KEYS = {"generator_provenance", "input_provenance"}
REL_TOL, ABS_TOL = 1e-6, 1e-9


def walk(a, b, path, problems):
    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            problems.append(f"{path}: keys differ {sorted(a.keys() ^ b.keys())}")
            return
        for k in a:
            if k not in SKIP_KEYS:
                walk(a[k], b[k], f"{path}/{k}", problems)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            problems.append(f"{path}: length {len(a)} vs {len(b)}")
            return
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f"{path}[{i}]", problems)
    elif isinstance(a, float) or isinstance(b, float):
        if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
            problems.append(f"{path}: {a!r} vs {b!r}")
        elif abs(a - b) > max(REL_TOL * max(abs(a), abs(b)), ABS_TOL):
            problems.append(f"{path}: {a!r} vs {b!r}")
    elif a != b:
        problems.append(f"{path}: {a!r} vs {b!r}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default="HEAD", help="git revision holding the committed results")
    args = ap.parse_args()

    listed = subprocess.check_output(["git", "-C", str(REPO), "ls-tree", "--name-only",
                                      args.ref, "results/tables/"], text=True).split()
    committed = [p for p in listed if p.endswith(".json")]
    if not committed:
        print(f"no committed results/tables/*.json at {args.ref}")
        return 1
    bad = 0
    for rel in committed:
        new_path = REPO / rel
        if not new_path.exists():
            print(f"MISSING  {rel} (committed but not regenerated)")
            bad += 1
            continue
        old = json.loads(subprocess.check_output(["git", "-C", str(REPO), "show",
                                                  f"{args.ref}:{rel}"], text=True))
        new = json.loads(new_path.read_text())
        problems: list[str] = []
        walk(new, old, "", problems)
        if problems:
            bad += 1
            print(f"DRIFT    {rel}: {len(problems)} value(s)")
            for p in problems[:15]:
                print(f"           {p}")
        else:
            print(f"same     {rel}")
    print("\ncommitted results " + ("REPRODUCED" if not bad else f"NOT reproduced ({bad} file(s))"))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
