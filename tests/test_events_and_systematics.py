"""hicpy.events (the shared reader) and the fit-systematics script."""
import json
import subprocess
import sys

import awkward as ak
import numpy as np

from conftest import REPO
from hicpy import kinematics as kin
from hicpy.events import read_particles, track_selection


def test_read_particles_is_float64_with_shared_eta(small_sample):
    ev = read_particles(small_sample)
    for a in (ev.px, ev.py, ev.pz, ev.E, ev.pt, ev.eta):
        assert ak.flatten(a).layout.dtype == np.float64
    np.testing.assert_array_equal(ak.flatten(ev.eta), ak.flatten(kin.eta(ev.px, ev.py, ev.pz)))
    assert ev.provenance["config"]["seed"] == 20260827
    assert len(ev.px) == ev.n_events


def test_track_selection_matches_its_definition(small_sample):
    ev = read_particles(small_sample)
    cuts = {"pt_min": 0.15, "pt_max": 50.0, "eta_abs_max": 0.8}
    sel = ak.flatten(track_selection(ev, cuts))
    pt, eta, q = (ak.to_numpy(ak.flatten(x)) for x in (ev.pt, ev.eta, ev.charge))
    manual = (q != 0) & (np.abs(eta) < 0.8) & (pt > 0.15) & (pt < 50.0)
    np.testing.assert_array_equal(ak.to_numpy(sel), manual)


def test_fit_systematics_runs_and_separates_stat_and_syst(small_sample, tmp_path):
    r = subprocess.run([sys.executable, str(REPO / "analysis" / "m01_fit_systematics.py"),
                        "--infile", str(small_sample), "--n-bootstrap", "5",
                        "--tabdir", str(tmp_path), "--figdir", str(tmp_path)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    out = json.loads((tmp_path / "m01_tsallis_systematics.json").read_text())
    for key in ("inclusive", "pion_control"):
        s = out[key]["summary"]
        for par in ("T", "n"):
            assert s[par]["stat"] > 0 and s[par]["syst"] >= 0
        rows = out[key]["variations"]
        assert rows[0]["group"] == "nominal"
        assert {r["group"] for r in rows[1:]} == {"range_low", "range_high", "binning"}
        for row in rows[1:]:
            assert row["T_n_bootstrap"] >= 3
    assert (tmp_path / "m01_tsallis_systematics.png").exists()
