"""End-to-end tests of the Module 1 chain on a small sample (no ROOT needed):

    generate -> ROOT file -> uproot analysis -> JSON summary -> closure test

These run in CI and guard the analysis against regressions once the toy
generator is replaced by PYTHIA (Module 2): the analysis code must keep giving
the same answers for the same input.
"""
import json
import pathlib
import subprocess
import sys

import numpy as np
import pytest
import uproot

from conftest import REPO, N_EVENTS_SMALL


def run(script, *args):
    return subprocess.run([sys.executable, str(REPO / script), *map(str, args)],
                          capture_output=True, text=True)


def test_root_file_layout_and_provenance(small_sample):
    with uproot.open(small_sample) as f:
        assert "events" in f and "provenance" in f
        t = f["events"]
        assert t.num_entries == N_EVENTS_SMALL
        # NanoAOD-style: counter branch + jagged particle branches, classic TTree
        assert t.classname == "TTree"
        for b in ("nparticle", "particle_px", "particle_py", "particle_pz",
                  "particle_E", "particle_pdg", "has_z"):
            assert b in t.keys(), b
        prov = json.loads(str(f["provenance"]))
    assert prov["config"]["n_events"] == N_EVENTS_SMALL
    assert prov["config"]["seed"] == 20260827
    assert prov["generator"] == "hicpy.toy_generator"
    assert json.loads(small_sample.with_suffix(".json").read_text()) == prov


@pytest.fixture(scope="module")
def analysed(small_sample, tmp_path_factory):
    out = tmp_path_factory.mktemp("analysis")
    r = run("analysis/m01_analysis_uproot.py", "--infile", small_sample,
            "--figdir", out / "fig", "--tabdir", out / "tab")
    assert r.returncode == 0, r.stderr
    return out


def test_analysis_outputs_exist(analysed):
    for stem in ("m01_multiplicity", "m01_dndeta", "m01_pt_spectrum", "m01_dimuon_mass"):
        assert (analysed / "fig" / f"{stem}.pdf").exists()
        assert (analysed / "fig" / f"{stem}.png").exists()
    assert (analysed / "tab" / "m01_fit_results.json").exists()


def test_summary_schema(analysed):
    s = json.loads((analysed / "tab" / "m01_fit_results.json").read_text())
    assert s["arm"] == "uproot"
    assert s["n_events"] == N_EVENTS_SMALL
    for k in ("mean_Nch_in_acceptance", "dNdeta_central", "n_dimuon_candidates",
              "mean_dimuon_mass_GeV", "tsallis", "histograms", "generator_provenance"):
        assert k in s, k
    for h in ("mult", "eta", "pt", "mass"):
        rec = s["histograms"][h]
        assert len(rec["edges"]) == len(rec["contents"]) + 1 == len(rec["errors"]) + 1
    assert all(f["converged"] for f in s["fits"])


def test_histograms_are_properly_normalised(analysed):
    """P(N_ch) integrates to 1; dN/deta integrates to the mean charged
    multiplicity over the full eta range; the pT yield is positive."""
    s = json.loads((analysed / "tab" / "m01_fit_results.json").read_text())
    h = s["histograms"]["mult"]
    w = np.diff(h["edges"])
    assert np.sum(np.array(h["contents"]) * w) == pytest.approx(1.0, abs=1e-9)
    h = s["histograms"]["eta"]
    w = np.diff(h["edges"])
    total = np.sum(np.array(h["contents"]) * w)
    nbar = s["generator_provenance"]["config"]["nbd_nbar"]
    assert total == pytest.approx(nbar, rel=0.05)      # all charged hadrons, pT > pt_min
    assert np.all(np.array(s["histograms"]["pt"]["contents"]) >= 0)


def test_figures_are_byte_reproducible(small_sample, tmp_path):
    """Same input, same code -> identical PDF and PNG bytes (no timestamps)."""
    for d in ("a", "b"):
        r = run("analysis/m01_analysis_uproot.py", "--infile", small_sample,
                "--figdir", tmp_path / d, "--tabdir", tmp_path / d)
        assert r.returncode == 0, r.stderr
    for f in sorted((tmp_path / "a").glob("*.p*")):
        assert f.read_bytes() == (tmp_path / "b" / f.name).read_bytes(), f.name


def test_closure_passes_on_small_sample(small_sample):
    r = run("analysis/m01_closure_test.py", "--infile", small_sample)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "CLOSURE PASSED" in r.stdout
