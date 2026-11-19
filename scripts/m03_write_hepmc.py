#!/usr/bin/env python3
"""Module 3, step 1: generate PYTHIA events and write the FULL history as HepMC3.

    python scripts/m03_write_hepmc.py [--n-events 100] [--card generators/pythia_mb_13TeV.cmnd]
                                      [--out data/generated/m03_pythia_mb_13TeV.hepmc]

Module 2 kept only final-state particles (a flat list per event).  Here every
entry of the PYTHIA event record is kept -- beams, hard process, MPI, showers,
remnants, string fragmentation, decays -- and written as a HepMC3 graph with
hicpy.eventrecord.to_hepmc (one vertex per distinct set of mothers).

The same card and seed as Module 2 are used, so event N here is event N there.
The .hepmc file is large (~1 MB per 10 events) and gitignored; the provenance
JSON next to it is what gets committed.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

from hicpy.eventrecord import records_from_pythia, to_hepmc   # noqa: E402
from m02_run_pythia import card_value                        # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--card", default=str(REPO / "generators" / "pythia_mb_13TeV.cmnd"))
    ap.add_argument("--n-events", type=int, default=100)
    ap.add_argument("--out", default=str(REPO / "data" / "generated" / "m03_pythia_mb_13TeV.hepmc"))
    args = ap.parse_args()
    try:
        import pythia8
    except ImportError:
        import pythia8mc as pythia8
    import pyhepmc

    card = pathlib.Path(args.card); card_text = card.read_text()
    out = pathlib.Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    py = pythia8.Pythia("", False)
    if not py.readFile(str(card)):
        print(f"PYTHIA rejected {card}"); return 1
    py.readString("Print:quiet = on")            # silence; Next:* irrelevant here
    if not py.init():
        print("PYTHIA init() failed"); return 1

    t0, n_ok = time.time(), 0
    with pyhepmc.open(str(out), "w") as f:
        while n_ok < args.n_events:
            if not py.next():
                continue
            ev = to_hepmc(records_from_pythia(py.event), event_number=n_ok)
            f.write(ev)
            n_ok += 1
    try:
        commit = subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                                         stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        commit = "unknown"
    prov = {"generator": "pythia8", "pythia_version": py.parm("Pythia:versionNumber"),
            "module": pythia8.__name__, "card_file": str(card.relative_to(REPO)) if card.is_relative_to(REPO) else str(card),
            "seed": int(card_value(card_text, "Random:seed", 0)), "n_events": n_ok,
            "converter": "hicpy.eventrecord.to_hepmc (one vertex per mother set)",
            "pyhepmc": pyhepmc.__version__, "git_commit": commit,
            "produced_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    out.with_suffix(".json").write_text(json.dumps(prov, indent=2))
    print(f"[out] {out} ({out.stat().st_size / 1e6:.1f} MB), {n_ok} events in "
          f"{time.time() - t0:.1f} s; provenance -> {out.with_suffix('.json').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
