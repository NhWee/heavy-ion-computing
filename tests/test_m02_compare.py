"""analysis/m02_compare_data.py on the toy sample (no PYTHIA, no network).

The "measurement" is synthesised from the MC itself, through the HEPData YAML
format and the provenance JSON, so the loader, the event-class selection,
the binning, the ratio and the verdict are all exercised and the answer is
known: MC / data = 1 exactly.
"""
import json
import subprocess
import sys

import numpy as np
import pytest
import yaml

from conftest import REPO

sys.path.insert(0, str(REPO / "analysis"))
import m02_compare_data as cmp   # noqa: E402

SCRIPT = REPO / "analysis" / "m02_compare_data.py"
CFG = REPO / "config" / "analysis_m02_pythia.yaml"


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                          capture_output=True, text=True)


def test_mc_only_runs_on_toy(small_sample, tmp_path):
    r = run("--config", CFG, "--mc", small_sample, "--no-data",
            "--figdir", tmp_path, "--tabdir", tmp_path)
    assert r.returncode == 0, r.stderr
    s = json.loads((tmp_path / "m02_compare_summary.json").read_text())
    assert s["event_class"] == "INEL>0" and s["verdict"] is None
    assert s["n_events_selected"] <= s["n_events_generated"]
    assert (tmp_path / "m02_dndeta_vs_data.pdf").exists()


def test_event_classes(small_sample):
    mc = cmp.read_mc(small_sample)
    inel = cmp.select_events(mc, "INEL")
    inel0 = cmp.select_events(mc, "INEL>0")
    assert inel.all() and inel0.sum() <= inel.sum()
    with pytest.raises(ValueError):           # toy has no process_code branch
        cmp.select_events(mc, "NSD")
    with pytest.raises(ValueError):
        cmp.select_events(mc, "bogus")


def test_dndeta_normalisation(small_sample):
    """Integral of dN/deta over the full range = mean charged multiplicity."""
    mc = cmp.read_mc(small_sample)
    mask = cmp.select_events(mc, "INEL")
    edges = np.linspace(-10, 10, 101)        # wide enough to hold every particle
    y, e = cmp.dndeta(mc, mask, edges)
    import awkward as ak
    nch = ak.to_numpy(ak.sum(mc["charged"], axis=1)).mean()
    assert np.sum(y * np.diff(edges)) == pytest.approx(nch, rel=1e-9)
    # binning-free central density = count in |eta| < 0.5 per event / 1.0
    c, ce = cmp.central_density(mc, mask)
    n05 = ak.to_numpy(ak.sum(mc["charged"] & (abs(mc["eta"]) < 0.5), axis=1))
    assert c == pytest.approx(n05.mean(), rel=1e-12)
    assert ce == pytest.approx(n05.std(ddof=1) / np.sqrt(len(n05)), rel=1e-9)


def _write_fake_dataset(tmp_path, edges, y, stat, sys_):
    """HEPData-style YAML table + filled provenance JSON (synthetic, for the test)."""
    table = {
        "independent_variables": [{"header": {"name": "ETARAP"}, "values": [
            {"low": float(lo), "high": float(hi)} for lo, hi in zip(edges[:-1], edges[1:])]}],
        "dependent_variables": [{"header": {"name": "DN/DETARAP"}, "values": [
            {"value": float(v), "errors": [{"symerror": float(s), "label": "stat"},
                                           {"asymerror": {"plus": float(u), "minus": -float(u)},
                                            "label": "sys"}]}
            for v, s, u in zip(y, stat, sys_)]}],
    }
    (tmp_path / "fake_table.yaml").write_text(yaml.safe_dump(table))
    prov = {"experiment": "SYNTHETIC-TEST", "sqrt_s_TeV": 13.0, "observable": "dN_ch/deta",
            "event_class": "INEL>0", "file": "fake_table.yaml"}
    pj = tmp_path / "fake.json"
    pj.write_text(json.dumps(prov))
    return pj


def test_self_comparison_gives_ratio_one(small_sample, tmp_path):
    mc = cmp.read_mc(small_sample)
    mask = cmp.select_events(mc, "INEL>0")
    edges = np.linspace(-1.5, 1.5, 16)
    y, e = cmp.dndeta(mc, mask, edges)
    pj = _write_fake_dataset(tmp_path, edges, y, e, 0.05 * y)
    r = run("--config", CFG, "--mc", small_sample, "--data", pj,
            "--figdir", tmp_path, "--tabdir", tmp_path)
    assert r.returncode == 0, r.stderr + r.stdout
    s = json.loads((tmp_path / "m02_compare_summary.json").read_text())
    v = s["verdict"]
    assert v["ratio"] == pytest.approx(1.0, abs=1e-9)
    assert v["within_tolerance"] is True
    assert v["data_sys"] == pytest.approx(0.05 * v["data"], rel=1e-9)
    assert s["data"]["provenance"]["experiment"] == "SYNTHETIC-TEST"


def test_csv_loader_equivalent_to_yaml(tmp_path):
    edges = np.array([-1.0, 0.0, 1.0]); y = np.array([5.0, 5.2])
    pj = _write_fake_dataset(tmp_path, edges, y, [0.1, 0.1], [0.3, 0.3])
    a = cmp.load_hepdata(tmp_path / "fake_table.yaml")
    (tmp_path / "t.csv").write_text(
        "#: name: Table 1\nx,xlo,xhi,y,stat+,stat-,sys+,sys-\n"
        "-0.5,-1.0,0.0,5.0,0.1,-0.1,0.3,-0.3\n0.5,0.0,1.0,5.2,0.1,-0.1,0.3,-0.3\n")
    b = cmp.load_hepdata(tmp_path / "t.csv")
    for k in ("x", "xlo", "xhi", "y", "stat", "sys"):
        np.testing.assert_allclose(a[k], b[k])


def test_unverified_provenance_is_refused(tmp_path):
    pj = tmp_path / "p.json"
    pj.write_text(json.dumps({"experiment": "ALICE", "sqrt_s_TeV": 13.0, "verified": False,
                              "observable": "dN_ch/deta", "event_class": "INEL>0", "file": "x.yaml"}))
    with pytest.raises(ValueError, match="verified"):
        cmp.load_dataset(pj)


def test_repository_provenance_json_is_consistent():
    """The committed ALICE JSON must be refused until verified, and its fields
    must be the ones the fetch script and the comparison expect."""
    pj = REPO / "data" / "raw" / "alice_pp13TeV_dndeta.json"
    prov = json.loads(pj.read_text())
    for k in ("experiment", "arxiv", "hepdata", "sqrt_s_TeV", "observable", "event_class", "file"):
        assert k in prov
    assert prov["hepdata"]["record"].startswith("ins")
    if not prov["verified"]:
        assert prov["sha256"] is None and prov["downloaded_utc"] is None
        with pytest.raises(ValueError, match="verified"):
            cmp.load_dataset(pj)
    else:
        assert prov["sha256"] and prov["downloaded_utc"]


def test_incomplete_provenance_is_refused(tmp_path):
    pj = tmp_path / "p.json"
    pj.write_text(json.dumps({"experiment": "ALICE", "sqrt_s_TeV": 13.0,
                              "observable": "<fill>", "event_class": "INEL>0", "file": "x.yaml"}))
    with pytest.raises(ValueError, match="incomplete"):
        cmp.load_dataset(pj)
