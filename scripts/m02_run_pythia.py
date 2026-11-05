#!/usr/bin/env python3
"""Module 2, step 1: run PYTHIA 8 from a card and WRITE THE SAME ROOT LAYOUT AS MODULE 1.

    python scripts/m02_run_pythia.py --card generators/pythia_mb_13TeV.cmnd \
        [--n-events 1000] [--out data/generated/m02_pythia_mb_13TeV.root] [--hepmc]

STATUS: written for the pythia8 Python bindings of conda-forge `pythia8`
(PYTHIA 8.3 API) but NOT executed in the bootstrap container (no PYTHIA
there).  Run it inside the WSL2 `hic` environment; record API drift in
docs/troubleshooting.md.

The point of this script
------------------------
Module 1's analysis chain (analysis/m01_analysis_uproot.py, the ROOT arms,
the tests) must run on this file WITHOUT MODIFICATION.  So the TTree has the
same branches as scripts/m01_make_toy_events.py writes:

    nparticle                      counter (written by uproot's mktree)
    particle_px/py/pz/E   float32  per final-state particle, GeV
    particle_pdg          int32    PDG Monte Carlo ID
    has_z                 bool     a Z boson (id 23) appears in the event record

plus branches Module 1 does not read but Module 2 needs:

    particle_charge       int8     3 * electric charge is PYTHIA's convention;
                                   here the charge itself (-2..+2)
    process_code          int32    PYTHIA process code: 101 ND, 103/104 SD,
                                   105 DD, 106 CD -- lets the analysis select
                                   NSD, INEL>0, ... afterwards
    n_ch_eta1             int32    charged particles with |eta| < 1 (INEL>0 flag)
    weight                float64  event weight (1 for unweighted MB)

"Final state" = PYTHIA isFinal() after decays limited by
ParticleDecays:tau0Max in the card, i.e. the experimental primary-particle
definition.  Neutrinos are kept (they are final-state particles; the analysis
selects on pdg/charge), beam remnants are not special-cased.

Provenance written into the file (TObjString "provenance") and next to it as
JSON: the FULL card text, the PYTHIA version, the seed, the cross section and
error PYTHIA reports at the end (sigmaGen), the number of events tried /
accepted, git commit, platform and package versions.  The cross section is
what turns "events" into "mb" when comparing to data.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import platform
import re
import subprocess
import sys
import time

import awkward as ak
import numpy as np
import uproot

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent


def git_hash(repo: pathlib.Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def card_value(card_text: str, key: str, default=None):
    """Read `key = value` from the card text (first match, comments stripped)."""
    for line in card_text.splitlines():
        line = line.split("!", 1)[0].strip()
        m = re.match(rf"^{re.escape(key)}\s*=\s*(\S+)", line)
        if m:
            return m.group(1)
    return default


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--card", default=str(REPO / "generators" / "pythia_mb_13TeV.cmnd"))
    ap.add_argument("--n-events", type=int, default=None,
                    help="override Main:numberOfEvents from the card")
    ap.add_argument("--seed", type=int, default=None,
                    help="override Random:seed from the card (for independent samples)")
    ap.add_argument("--out", default=None,
                    help="output ROOT file (default data/generated/<card stem>.root)")
    ap.add_argument("--hepmc", action="store_true",
                    help="also write a HepMC3 file next to the ROOT file (Module 3)")
    args = ap.parse_args()

    try:
        import pythia8
    except ImportError:
        print("pythia8 Python bindings not found -- run inside the WSL2 `hic` "
              "environment (environment/setup_wsl2.sh)")
        return 2

    card = pathlib.Path(args.card)
    card_text = card.read_text()
    out = (pathlib.Path(args.out) if args.out
           else REPO / "data" / "generated" / f"m02_{card.stem}.root")
    out.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------ configure
    pythia = pythia8.Pythia()
    if not pythia.readFile(str(card)):
        print(f"PYTHIA rejected the card {card}")
        return 1
    n_events = args.n_events or int(card_value(card_text, "Main:numberOfEvents", 1000))
    seed = args.seed or int(card_value(card_text, "Random:seed", 0))
    pythia.readString(f"Main:numberOfEvents = {n_events}")
    pythia.readString("Random:setSeed = on")
    pythia.readString(f"Random:seed = {seed}")
    if not pythia.init():
        print("PYTHIA init() failed -- see the log above")
        return 1

    hepmc_writer = None
    if args.hepmc:
        import pyhepmc
        hepmc_path = out.with_suffix(".hepmc")
        hepmc_writer = pyhepmc.io.WriterAscii(str(hepmc_path))
        to_hepmc = pythia8.Pythia8ToHepMC3() if hasattr(pythia8, "Pythia8ToHepMC3") else None
        if to_hepmc is None:
            print("this pythia8 build has no HepMC3 interface; --hepmc ignored")
            hepmc_writer = None

    # ------------------------------------------------------------ generate
    px_l, py_l, pz_l, e_l, id_l, q_l = [], [], [], [], [], []
    has_z, proc, nch1, wgt = [], [], [], []
    n_tried = 0
    t0 = time.time()
    print(f"[gen] {n_events} events, seed = {seed}, card = {card.name}")
    max_tries = 2 * n_events + 100        # never spin forever on a broken setup
    while len(proc) < n_events:
        n_tried += 1
        if n_tried > max_tries:
            print(f"[gen] giving up: {n_tried} tries for {len(proc)} events")
            break
        if not pythia.next():
            continue                      # failed event: PYTHIA counts errors itself
        ev = pythia.event
        px, py, pz, e, pid, q = [], [], [], [], [], []
        z_seen = False
        n_ch_eta1 = 0
        for p in ev:
            if abs(p.id()) == 23:
                z_seen = True
            if not p.isFinal():
                continue
            px.append(p.px()); py.append(p.py()); pz.append(p.pz()); e.append(p.e())
            pid.append(p.id()); q.append(int(round(p.charge())))
            if p.isCharged() and abs(p.eta()) < 1.0:
                n_ch_eta1 += 1
        px_l.append(np.asarray(px, np.float32)); py_l.append(np.asarray(py, np.float32))
        pz_l.append(np.asarray(pz, np.float32)); e_l.append(np.asarray(e, np.float32))
        id_l.append(np.asarray(pid, np.int32)); q_l.append(np.asarray(q, np.int8))
        has_z.append(z_seen)
        info = pythia.infoPython() if hasattr(pythia, "infoPython") else pythia.info
        proc.append(int(info.code()))
        nch1.append(n_ch_eta1)
        wgt.append(float(info.weight()))
        if hepmc_writer is not None:
            hev = pyhepmc.GenEvent()
            to_hepmc.fill_next_event(pythia, hev)
            hepmc_writer.write_event(hev)
    dt = time.time() - t0
    pythia.stat()                         # cross-section table into the log
    if hepmc_writer is not None:
        hepmc_writer.close()

    info = pythia.infoPython() if hasattr(pythia, "infoPython") else pythia.info
    n_acc = len(proc)
    print(f"[gen] {n_acc} events accepted / {n_tried} tried in {dt:.1f} s "
          f"({n_acc / max(dt, 1e-9):.1f} ev/s)")
    print(f"[gen] sigma_gen = {info.sigmaGen():.4f} +- {info.sigmaErr():.4f} mb")

    # ------------------------------------------------------------ write
    particle = ak.zip({
        "px": ak.Array(px_l), "py": ak.Array(py_l), "pz": ak.Array(pz_l),
        "E": ak.Array(e_l), "pdg": ak.Array(id_l), "charge": ak.Array(q_l),
    })
    arrays = {
        "particle": particle,
        "has_z": np.asarray(has_z, bool),
        "process_code": np.asarray(proc, np.int32),
        "n_ch_eta1": np.asarray(nch1, np.int32),
        "weight": np.asarray(wgt, np.float64),
    }
    provenance = {
        "generator": "pythia8",
        "pythia_version": str(info.pythiaVersion() if hasattr(info, "pythiaVersion")
                              else pythia.settings.parm("Pythia:versionNumber")),
        "card_file": str(card.relative_to(REPO)) if card.is_relative_to(REPO) else str(card),
        "card_text": card_text,
        "config": {"n_events": n_acc, "seed": seed,
                   "eCM_GeV": float(card_value(card_text, "Beams:eCM", "nan")),
                   "tune_pp": card_value(card_text, "Tune:pp"),
                   "pdf_pSet": card_value(card_text, "PDF:pSet"),
                   "process": "SoftQCD:inelastic",
                   "tau0Max_mm": float(card_value(card_text, "ParticleDecays:tau0Max", "nan"))},
        "cross_section_mb": {"sigma_gen": info.sigmaGen(), "sigma_err": info.sigmaErr(),
                             "n_tried": n_tried, "n_accepted": n_acc},
        "produced_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_commit": git_hash(REPO),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "versions": {"numpy": np.__version__, "awkward": ak.__version__,
                     "uproot": uproot.__version__},
    }
    # same TTree layout as Module 1 (see docs/troubleshooting.md #1 and #2)
    with uproot.recreate(out) as f:
        branch_types = {k: (v.type if hasattr(v, "type") else v.dtype)
                        for k, v in arrays.items()}
        tree = f.mktree("events", branch_types, title="PYTHIA 8 final-state particles",
                        counter_name=lambda name: "n" + name)
        tree.extend(arrays)
        f["provenance"] = json.dumps(provenance, indent=2)
    out.with_suffix(".json").write_text(json.dumps(provenance, indent=2))

    n_part = int(sum(len(a) for a in px_l))
    print(f"[out] {out}  ({out.stat().st_size / 1e6:.2f} MB)")
    n_inel0 = int(np.count_nonzero(np.asarray(nch1) > 0))
    print(f"[out] {n_acc} events, {n_part} final-state particles, "
          f"{n_inel0} INEL>0 ({100.0 * n_inel0 / max(n_acc, 1):.1f} %), "
          f"{int(np.sum(has_z))} with a Z")
    print(f"[out] provenance -> {out.with_suffix('.json').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
