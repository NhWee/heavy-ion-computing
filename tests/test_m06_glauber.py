"""hicpy.glauber against limits that follow from the definitions."""
import numpy as np
import pytest

from hicpy.glauber import PB208, PROTON, Nucleus, sample_nucleus, collide, run, centrality_table


def test_sampled_density_follows_woods_saxon():
    rng = np.random.default_rng(1)
    r = np.concatenate([np.linalg.norm(sample_nucleus(rng, PB208, recentre=False), axis=1) for _ in range(50)])
    grid = np.linspace(0, 15, 6001)
    w = grid ** 2 / (1 + np.exp((grid - PB208.R) / PB208.a))
    expected_r2 = np.trapezoid(w * grid ** 2, grid) / np.trapezoid(w, grid)
    # d_min pushes nucleons slightly outwards; allow 1 %
    assert np.mean(r ** 2) == pytest.approx(expected_r2, rel=0.01)


def test_minimum_distance_is_respected():
    pos = sample_nucleus(np.random.default_rng(2), PB208)
    d = np.linalg.norm(pos[:, None] - pos[None, :], axis=2)[np.triu_indices(208, 1)]
    assert d.min() >= PB208.d_min - 1e-12


def test_head_on_and_far_apart_limits():
    rng = np.random.default_rng(3)
    a, b = sample_nucleus(rng, PB208), sample_nucleus(rng, PB208)
    npart0, ncoll0, _, _ = collide(a, b, 0.0, 67.6)
    assert npart0 > 380 and ncoll0 > npart0 / 2
    assert collide(a, b, 25.0, 67.6)[:2] == (0, 0)


def test_point_projectile_npart_is_ncoll_plus_one():
    """p-A: the proton is wounded once, every collision wounds one target
    nucleon, so N_part = N_coll + 1 in every event."""
    res = run(400, PROTON, PB208, 67.6, seed=4)
    np.testing.assert_array_equal(res["npart"], res["ncoll"] + 1)


def test_pp_cross_section_is_sigma_nn():
    """Two point nucleons: sigma_inel must equal sigma_NN (67.6 mb = 0.0676 b)."""
    res = run(200000, PROTON, PROTON, 67.6, b_max=2.0, seed=5)
    assert res["sigma_inel_b"] == pytest.approx(0.0676, abs=3 * res["sigma_inel_err_b"])


def test_centrality_classes_are_ordered():
    res = run(600, PB208, PB208, 67.6, seed=6)
    t = centrality_table(res, (0, 10, 50, 100))
    assert t[0]["npart"] > t[1]["npart"] > t[2]["npart"]
    assert t[0]["b"] < t[1]["b"] < t[2]["b"]


def test_m07_percentile_classes_and_frame():
    """Angantyr analysis helpers: percentile classes and the p-Pb lab frame."""
    import pathlib
    import sys
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "analysis"))
    import m07_angantyr as m7
    rng = np.random.default_rng(0)
    key = np.arange(1000, dtype=float)
    cls = m7.percentile_classes(key, [0, 10, 50, 100], rng)
    assert [len(i) for _, i in cls] == [100, 400, 500]
    assert key[cls[0][1]].min() > key[cls[1][1]].max()          # first class = largest key
    assert m7.Y_SHIFT["pPb"] == -0.465 and m7.FLIP["pPb"]          # proton to -eta, CM at y = -0.465
    assert not m7.FLIP["PbPb"] and m7.Y_SHIFT["PbPb"] == 0.0
