"""Shared fixtures: a small, seed-fixed toy sample written once per session.

Tests must never touch data/ or results/ of the repository; everything goes to
pytest's tmp_path_factory directory.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

N_EVENTS_SMALL = 3000     # ~36k pions in |eta|<0.8: enough for every fit to converge


@pytest.fixture(scope="session")
def small_sample(tmp_path_factory):
    """Generate a 3000-event ROOT file with the repository script (so the
    script itself -- argument parsing, TTree layout, provenance -- is tested)."""
    out = tmp_path_factory.mktemp("sample") / "toy_small.root"
    subprocess.run([sys.executable, str(REPO / "scripts" / "m01_make_toy_events.py"),
                    "--n-events", str(N_EVENTS_SMALL), "--out", str(out)],
                   check=True, capture_output=True)
    assert out.exists() and out.with_suffix(".json").exists()
    return out


@pytest.fixture(scope="session")
def analysed(small_sample, tmp_path_factory):
    """The uproot arm run once on the small sample; returns the output dir
    (figures in ./fig, tables in ./tab)."""
    out = tmp_path_factory.mktemp("analysis")
    r = subprocess.run([sys.executable, str(REPO / "analysis" / "m01_analysis_uproot.py"),
                        "--infile", str(small_sample),
                        "--figdir", str(out / "fig"), "--tabdir", str(out / "tab")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return out
