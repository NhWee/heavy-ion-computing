"""Unit tests for hicpy.kinematics -- the checks done by hand in Session 1,
now permanent."""
import numpy as np
import pytest

from hicpy import kinematics as kin

M_MU = 0.1056583755


def test_ptetaphim_round_trip():
    rng = np.random.default_rng(1)
    pt = rng.uniform(0.2, 50, 1000)
    eta = rng.uniform(-5, 5, 1000)
    phi = rng.uniform(-np.pi, np.pi, 1000)
    m = rng.choice([0.13957039, 0.493677, 0.93827208816], 1000)
    px, py, pz, e = kin.from_ptetaphim(pt, eta, phi, m)
    np.testing.assert_allclose(kin.pt(px, py), pt, rtol=1e-12)
    np.testing.assert_allclose(kin.eta(px, py, pz), eta, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(kin.phi(px, py), phi, rtol=1e-12)
    # m^2 = E^2 - p^2 is a catastrophic cancellation: at pT = 50 GeV, eta = 5
    # (p ~ 3.7 TeV) only ~8 significant digits of a 0.14 GeV mass survive in
    # float64 (docs/troubleshooting.md #6).  Assert that bound, not better.
    np.testing.assert_allclose(kin.invariant_mass(e, px, py, pz), m, atol=1e-6)


def test_back_to_back_muons_give_91p2_GeV():
    """Two 45.6 GeV muons back to back at eta = 0: m = 2 E = 91.2 GeV."""
    p = 45.6
    e = np.sqrt(p * p + M_MU * M_MU)
    m = kin.invariant_mass_pair((p, 0.0, 0.0, e), (-p, 0.0, 0.0, e))
    assert m == pytest.approx(2 * e, rel=1e-12)
    assert m == pytest.approx(91.2, abs=1e-3)


def test_rapidity_equals_eta_for_massless():
    px, py, pz, e = kin.from_ptetaphim(3.0, 1.7, 0.4, 0.0)
    assert kin.rapidity(e, pz) == pytest.approx(kin.eta(px, py, pz), rel=1e-12)


def test_rapidity_below_eta_for_massive():
    """|y| < |eta| for any massive particle (y -> eta only as m/pT -> 0)."""
    px, py, pz, e = kin.from_ptetaphim(0.5, 2.0, 0.0, 0.93827208816)
    assert abs(kin.rapidity(e, pz)) < abs(kin.eta(px, py, pz))


def test_eta_is_finite_for_collinear_particle():
    assert np.isfinite(kin.eta(0.0, 0.0, 10.0))
    assert np.isfinite(kin.eta(0.0, 0.0, 0.0))


def test_delta_phi_wraps_into_minus_pi_pi():
    d = kin.delta_phi(np.array([3.0, -3.0, 0.1]), np.array([-3.0, 3.0, 0.0]))
    assert np.all(d > -np.pi) and np.all(d <= np.pi)
    assert d[0] == pytest.approx(6.0 - 2 * np.pi)
    assert d[2] == pytest.approx(0.1)


def test_invariant_mass_clips_negative_roundoff_to_zero():
    assert kin.invariant_mass(1.0, 1.0 + 1e-17, 0.0, 0.0) == 0.0
