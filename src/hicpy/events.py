"""Reading a Module-1-layout event file and the standard track selection.

Every Python analysis script reads the same branches the same way: float32
storage promoted to float64 BEFORE any arithmetic (docs/troubleshooting.md
#12), pT and eta from hicpy.kinematics (eta = asinh(pz/pT), the definition all
three arms share), and the provenance JSON the generator stored in the file.
"""
from __future__ import annotations

import json
from types import SimpleNamespace

import awkward as ak
import numpy as np
import uproot

from . import kinematics as kin

PARTICLE_BRANCHES = ("particle_px", "particle_py", "particle_pz", "particle_E",
                     "particle_pdg", "particle_charge")


def read_particles(path, extra_branches=()):
    """Read the particle branches of `path`.

    Returns a namespace with px, py, pz, E (float64, jagged), pdg, charge,
    pt, eta, n_events, provenance (dict) and any requested extra branches.
    """
    with uproot.open(path) as f:
        tree = f["events"]
        arr = tree.arrays(list(PARTICLE_BRANCHES) + list(extra_branches))
        prov = json.loads(str(f["provenance"]))
        n_events = tree.num_entries
    f64 = lambda a: ak.values_astype(a, np.float64)          # noqa: E731
    ev = SimpleNamespace(px=f64(arr.particle_px), py=f64(arr.particle_py),
                         pz=f64(arr.particle_pz), E=f64(arr.particle_E),
                         pdg=arr.particle_pdg, charge=arr.particle_charge,
                         n_events=int(n_events), provenance=prov)
    for b in extra_branches:
        setattr(ev, b, arr[b])
    ev.pt = kin.pt(ev.px, ev.py)
    ev.eta = kin.eta(ev.px, ev.py, ev.pz)
    return ev


def track_selection(ev, track_cuts):
    """The N_ch / pT-spectrum selection: charged (charge != 0), |eta| <
    eta_abs_max, pt_min < pT < pt_max.  Same in all three Module 1 arms."""
    return ((ev.charge != 0)
            & (abs(ev.eta) < track_cuts["eta_abs_max"])
            & (ev.pt > track_cuts["pt_min"]) & (ev.pt < track_cuts["pt_max"]))


__all__ = ["PARTICLE_BRANCHES", "read_particles", "track_selection"]
