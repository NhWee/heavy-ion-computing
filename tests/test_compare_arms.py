"""The three-arm cross-check has never been run against real ROOT output (no
ROOT in CI).  These tests feed it ROOT-arm summaries synthesised from the
uproot summary and check its VERDICTS: identical inputs agree; a stale summary,
different cuts, a shifted bin or a missing arm must fail."""
import copy
import json
import subprocess
import sys

import pytest

from conftest import REPO


def run_compare(tabdir, *extra):
    return subprocess.run([sys.executable, str(REPO / "analysis" / "m01_compare_arms.py"),
                           "--tabdir", str(tabdir), *extra], capture_output=True, text=True)


def root_arm_from_uproot(up, arm):
    """What a ROOT arm writes, built from the uproot summary (same numbers)."""
    return {"arm": arm, "n_events": up["n_events"],
            "mean_Nch": up["mean_Nch_in_acceptance"],
            "dNdeta_central": up["dNdeta_central"],
            "n_dimuon_candidates": up["n_dimuon_candidates"],
            "mean_dimuon_mass_GeV": up["mean_dimuon_mass_GeV"],
            "cuts": {"track": {k: up["cuts"]["track"][k] for k in ("pt_min", "pt_max", "eta_abs_max")},
                     "muon": dict(up["cuts"]["muon"])},
            "tsallis": {k: up["tsallis"][k] for k in ("T", "T_err", "n", "n_err", "chi2", "ndf")},
            "histograms": copy.deepcopy(up["histograms"]),
            "input_provenance": copy.deepcopy(up["generator_provenance"]),
            "root_version": "synthetic"}


@pytest.fixture()
def tables(analysed, tmp_path):
    up = json.loads((analysed / "tab" / "m01_fit_results.json").read_text())
    (tmp_path / "m01_fit_results.json").write_text(json.dumps(up))
    for arm, fname in (("pyroot", "m01_pyroot_summary.json"),
                       ("root_macro", "m01_root_summary.json")):
        (tmp_path / fname).write_text(json.dumps(root_arm_from_uproot(up, arm)))
    return tmp_path


def edit(path, fn):
    d = json.loads(path.read_text())
    fn(d)
    path.write_text(json.dumps(d))


def test_identical_arms_agree(tables):
    r = run_compare(tables)
    assert r.returncode == 0, r.stdout
    assert "ARMS AGREE" in r.stdout and "same input file" in r.stdout


def test_stale_summary_from_another_file_fails(tables):
    edit(tables / "m01_root_summary.json",
         lambda d: d["input_provenance"].update(produced_utc="2026-01-01T00:00:00Z"))
    r = run_compare(tables)
    assert r.returncode == 1 and "DIFFERENT INPUT FILES" in r.stdout


def test_summary_without_provenance_fails(tables):
    edit(tables / "m01_pyroot_summary.json", lambda d: d.pop("input_provenance"))
    r = run_compare(tables)
    assert r.returncode == 1 and "input provenance MISSING" in r.stdout


def test_different_cut_fails(tables):
    edit(tables / "m01_pyroot_summary.json", lambda d: d["cuts"]["track"].update(pt_max=60.0))
    r = run_compare(tables)
    assert r.returncode == 1 and "DIFFERENT CUTS" in r.stdout


def test_shifted_bin_fails(tables):
    def shift(d):
        h = d["histograms"]["pt"]
        h["contents"][5] *= 1.05
    edit(tables / "m01_root_summary.json", shift)
    r = run_compare(tables)
    assert r.returncode == 1 and "MISMATCH" in r.stdout


def test_missing_arm_fails_unless_allowed(tables):
    (tables / "m01_root_summary.json").unlink()
    r = run_compare(tables)
    assert r.returncode == 1 and "NOT CHECKED: root_macro" in r.stdout
    r = run_compare(tables, "--allow-missing")
    assert r.returncode == 0 and "NOT a full three-arm check" in r.stdout
