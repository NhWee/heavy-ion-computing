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
