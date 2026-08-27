#!/usr/bin/env python3
"""Module 1, step 1: generate a toy particle-level event sample and WRITE A ROOT FILE.

    python scripts/m01_make_toy_events.py [--config config/analysis.yaml]

What this teaches
-----------------
A collider event sample is *ragged*: event 1 has 137 particles, event 2 has 61.
A flat table cannot hold that, which is precisely why HEP uses ROOT TTrees with
variable-length branches instead of CSV.  Here we build such a TTree with
`uproot` (pure Python, no ROOT installation required) and store, per event:

    px, py, pz, E   -- variable-length float32 arrays, one entry per particle
    pdg             -- variable-length int32 array (PDG Monte Carlo ID)
    n_particles     -- scalar int32
    has_z           -- scalar bool

The file also carries its own provenance: the full generator configuration is
written next to it as JSON *and* embedded in the ROOT file as a TObjString, so
that a file found on disk in six months can still be traced back to the code
and the random seed that produced it.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import platform
import subprocess
import sys
import time

import awkward as ak
import numpy as np
import uproot
import yaml

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "src"))

from hicpy.toy_generator import ToyConfig, generate  # noqa: E402


def git_hash(repo: pathlib.Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=str(REPO / "config" / "analysis.yaml"))
    ap.add_argument("--n-events", type=int, default=None,
                    help="override generator.n_events (useful for quick tests)")
    ap.add_argument("--out", default=None, help="override output ROOT file path")
    ap.add_argument("--format", choices=["ttree", "rntuple"], default="ttree",
                    help="ROOT container: classic TTree (default, universally "
                         "readable) or the new RNTuple (needs ROOT >= 6.34)")
    args = ap.parse_args()

    cfg_all = yaml.safe_load(open(args.config))
    gen_cfg = dict(cfg_all["generator"])
    gen_cfg.pop("name", None)
    if args.n_events is not None:
        gen_cfg["n_events"] = args.n_events
    cfg = ToyConfig(**gen_cfg)

    out = pathlib.Path(args.out) if args.out else REPO / cfg_all["paths"]["generated_file"]
    out.parent.mkdir(parents=True, exist_ok=True)

    print(f"[gen] {cfg.n_events} events, seed = {cfg.seed}")
    t0 = time.time()
    ev = generate(cfg)
    dt = time.time() - t0
    print(f"[gen] done in {dt:.1f} s "
          f"({cfg.n_events / max(dt, 1e-9):.0f} events/s)")

    # NanoAOD-style layout: one *record* branch "particle" holding the
    # per-particle fields, plus a shared counter branch "nparticle" that ROOT
    # uses to know how long each event's arrays are.  This is what real CMS
    # NanoAOD files look like (nMuon / Muon_pt / Muon_eta / ...).
    particle = ak.zip({
        "px":  ak.values_astype(ak.Array(ev["px"]),  np.float32),
        "py":  ak.values_astype(ak.Array(ev["py"]),  np.float32),
        "pz":  ak.values_astype(ak.Array(ev["pz"]),  np.float32),
        "E":   ak.values_astype(ak.Array(ev["E"]),   np.float32),
        "pdg": ak.values_astype(ak.Array(ev["pdg"]), np.int32),
    })
    arrays = {"particle": particle, "has_z": ev["has_z"]}

    provenance = {
        "generator": "hicpy.toy_generator",
        "config_file": str(pathlib.Path(args.config).name),
        "config": cfg.to_dict(),
        "produced_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_commit": git_hash(REPO),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "versions": {"numpy": np.__version__, "awkward": ak.__version__,
                     "uproot": uproot.__version__},
    }

    # NOTE (see docs/troubleshooting.md #1): since uproot 5.7, assigning a dict
    # of arrays with `f["events"] = arrays` writes a ROOT::RNTuple, which only
    # ROOT >= 6.34 can read.  We want a classic TTree so that the ROOT C++ and
    # PyROOT versions of this analysis work with any modern ROOT.  `mktree`
    # forces a TTree; `--format rntuple` opts into the new container instead.
    with uproot.recreate(out) as f:
        if args.format == "ttree":
            branch_types = {k: (v.type if hasattr(v, "type") else v.dtype)
                            for k, v in arrays.items()}
            tree = f.mktree("events", branch_types,
                            title="toy particle-level events",
                            counter_name=lambda name: "n" + name)
            tree.extend(arrays)
        else:
            f["events"] = arrays
        f["provenance"] = json.dumps(provenance, indent=2)

    (out.with_suffix(".json")).write_text(json.dumps(provenance, indent=2))

    size_mb = out.stat().st_size / 1e6
    n_part = int(np.sum(ev["n_particles"]))
    print(f"[out] {out}  ({size_mb:.2f} MB)")
    print(f"[out] {cfg.n_events} events, {n_part} particles, "
          f"{int(ev['has_z'].sum())} events with a Z")
    print(f"[out] provenance -> {out.with_suffix('.json').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
