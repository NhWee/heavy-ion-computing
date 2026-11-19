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


# ---------------------------------------------------------------- pT spectrum
import m02_compare_pt as cpt   # noqa: E402

SCRIPT_PT = REPO / "analysis" / "m02_compare_pt.py"


def test_pt_self_comparison_gives_ratio_one(small_sample, tmp_path):
    """The MC invariant yield, dressed up as a HEPData table, must give
    MC/data = 1 in every bin and for the integral."""
    mc = cmp.read_mc(small_sample)
    mask = cmp.select_events(mc, "INEL>0")
    edges = np.geomspace(0.15, 5.0, 13)
    y, e = cpt.invariant_yield(mc, mask, edges)
    pj = _write_fake_dataset(tmp_path, edges, y, e, 0.05 * y)
    r = subprocess.run([sys.executable, str(SCRIPT_PT), "--config", str(CFG), "--mc", str(small_sample),
                        "--data", str(pj), "--figdir", str(tmp_path), "--tabdir", str(tmp_path)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr + r.stdout
    s = json.loads((tmp_path / "m02_compare_pt_summary.json").read_text())
    np.testing.assert_allclose(s["ratio"]["per_bin"], 1.0, rtol=1e-12)
    assert s["ratio"]["integral"] == pytest.approx(1.0, rel=1e-12)


def test_pt_yield_integrates_to_dndeta(small_sample):
    """Integral of 2 pi pT x invariant yield over pT = dN/deta in |eta| < 0.8
    (exact per-particle Jacobian: use fine bins so the bin-centre factor is ~1)."""
    import awkward as ak
    mc = cmp.read_mc(small_sample)
    mask = cmp.select_events(mc, "INEL")
    edges = np.geomspace(0.15, 60.0, 4001)
    y, _ = cpt.invariant_yield(mc, mask, edges)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    integral = np.sum(y * 2 * np.pi * ctr * np.diff(edges))
    n = ak.to_numpy(ak.sum(mc["charged"] & (abs(mc["eta"]) < 0.8), axis=1)).mean() / 1.6
    assert integral == pytest.approx(n, rel=2e-3)


# ---------------------------------------------------------------- identified
import m02_compare_identified as cid   # noqa: E402


@pytest.mark.parametrize("invariant,average", [(False, False), (True, True)])
def test_identified_self_comparison(small_sample, tmp_path, invariant, average):
    """MC pion spectrum dressed up as a table in either form must come back
    with ratio 1: checks y (not eta), the per-INEL>0 normalisation and the
    form handling (invariant / average of charges)."""
    mc = cid.read_mc(small_sample)
    edges = np.array([0.2, 0.4, 0.6, 1.0, 2.0])
    form = {"invariant": invariant, "average": average}
    y, e = cid.mc_spectrum(mc, "pi", edges, form)
    header = ("(1/(2*PI*PT)) * D2(N)/DYRAP/DPT" if invariant else "D2(N)/DYRAP/DPT")
    if average:
        header = "(PI+ + PI-)/2 " + header
    table = {"independent_variables": [{"header": {"name": "PT"}, "values": [
                {"low": float(lo), "high": float(hi)} for lo, hi in zip(edges[:-1], edges[1:])]}],
             "dependent_variables": [{"header": {"name": header}, "qualifiers": [
                {"name": "EVENT CLASS", "value": "INEL>0"}], "values": [
                {"value": float(v), "errors": [{"symerror": float(s), "label": "stat"},
                                               {"symerror": float(0.05 * v), "label": "sys"}]}
                for v, s in zip(y, e)]}]}
    (tmp_path / "pi.yaml").write_text(yaml.safe_dump(table))
    (tmp_path / "pi.json").write_text(json.dumps({"experiment": "SYNTHETIC", "sqrt_s_TeV": 13.0,
        "observable": "pi spectrum", "event_class": "INEL>0", "file": "pi.yaml", "verified": True,
        "hepdata": {"record": "test"}}))
    det = cid.detect_form(tmp_path / "pi.json")
    assert det["invariant"] == invariant and det["average"] == average
    r = subprocess.run([sys.executable, str(REPO / "analysis" / "m02_compare_identified.py"), "--mc",
                        str(small_sample), "--data", f"pi={tmp_path / 'pi.json'}", "--figdir", str(tmp_path),
                        "--tabdir", str(tmp_path)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr + r.stdout
    s = json.loads((tmp_path / "m02_identified_summary.json").read_text())
    np.testing.assert_allclose(s["species"]["pi"]["ratio"], 1.0, rtol=1e-12)


def test_identified_inel0_derivation_recovers_known_spectrum(tmp_path):
    """Fake ins1784041 tables: S_i = c_i * S, R_i = c_i with sum f_i c_i = 1.
    Both the ratio method and the weighted sum must give back S exactly."""
    import importlib.util
    import yaml
    spec = importlib.util.spec_from_file_location("derive", REPO / "scripts" / "m02_derive_identified_inel0.py")
    der = importlib.util.module_from_spec(spec); spec.loader.exec_module(der)
    edges = [0.1, 0.2, 0.5, 1.0, 3.0, 10.0, 20.0]
    s_true = np.array([10.0, 5.0, 2.0, 0.5, 0.01, 1e-4])
    c = np.linspace(4.0, 0.3, 10); c = c / np.sum(der.FRAC * c)          # sum f_i c_i = 1
    x = [{"low": lo, "high": hi} for lo, hi in zip(edges[:-1], edges[1:])]

    def table(n, classes, values):
        dvs = [{"header": {"name": "y"}, "qualifiers": [{"name": "Mutliplicity", "value": f"V0M Class {k}"}],
                "values": [{"value": float(v), "errors": [{"label": "stat", "symerror": float(0.01 * v)},
                                                          {"label": "full syst.", "symerror": float(0.05 * v)}]}
                           for v in vals]} for k, vals in zip(classes, values)]
        (tmp_path / f"Table_{n}.yaml").write_text(yaml.safe_dump(
            {"independent_variables": [{"header": {"name": "PT"}, "values": x}], "dependent_variables": dvs}))

    first9, last = der.CLASSES[:9], der.CLASSES[9:]
    table(1, first9, [ci * s_true for ci in c[:9]]); table(2, last, [c[9] * s_true])
    table(7, first9, [np.full(6, ci) for ci in c[:9]]); table(8, last, [np.full(6, c[9])])
    rows, spread, closure = der.derive(tmp_path, (1, 2, 7, 8))
    assert len(rows) == 5                                    # the 10-20 GeV/c bin is cut (PT_MAX)
    np.testing.assert_allclose([r["y"] for r in rows], s_true[:5], rtol=1e-12)
    np.testing.assert_allclose(closure, 1.0, rtol=1e-12)
    assert spread.max() < 1e-12
    np.testing.assert_allclose([r["sys"] / r["y"] for r in rows], 0.05, rtol=1e-12)   # linear sum of 5 %
