"""Module 8 closure: toy dimuon events with known content -> mass, fit -> truth back."""
import numpy as np
import pytest

from hicpy.dimuon import PDG, crystal_ball, fit_window, pair_mass, toy_events


def test_pair_mass_reproduces_the_generated_mass():
    t = toy_events(np.random.default_rng(2), 5000, {"J/psi": 0.3, "Z": 0.1})
    pt, eta, phi = (t[k].reshape(-1, 2) for k in ("pt", "eta", "phi"))
    m = pair_mass(pt[:, 0], eta[:, 0], phi[:, 0], pt[:, 1], eta[:, 1], phi[:, 1])
    np.testing.assert_allclose(m, t["true_mass"], rtol=1e-9)


def test_crystal_ball_is_normalised():
    x = np.linspace(5, 12, 200001)
    assert np.trapezoid(crystal_ball(x, 9.46, 0.1, 1.5, 3.0), x) == pytest.approx(1.0, abs=2e-3)


def test_upsilon_fit_recovers_yields_and_mass_scale():
    rng = np.random.default_rng(11)
    comp = {"Upsilon(1S)": 0.10, "Upsilon(2S)": 0.03, "Upsilon(3S)": 0.015}
    t = toy_events(rng, 200000, comp, sigma_rel=0.01, bkg_slope=0.3)
    pt, eta, phi = (t[k].reshape(-1, 2) for k in ("pt", "eta", "phi"))
    m = pair_mass(pt[:, 0], eta[:, 0], phi[:, 0], pt[:, 1], eta[:, 1], phi[:, 1])
    edges = np.arange(8.5, 11.5001, 0.02)
    counts = np.histogram(m, edges)[0]
    r = fit_window(edges, counts, list(comp), "cb")
    for s, f in comp.items():
        truth = f * 200000
        n, e = r["params"][f"N_{s}"], r["errors"][f"N_{s}"]
        assert e > 0
        assert abs(n - truth) < 4 * e, (s, n, e, truth)
    assert r["params"]["scale"] == pytest.approx(1.0, abs=1e-3)
    assert r["params"]["sigma"] == pytest.approx(0.01 * PDG["Upsilon(1S)"], rel=0.1)
    ratio = r["params"]["N_Upsilon(2S)"] / r["params"]["N_Upsilon(1S)"]
    assert ratio == pytest.approx(0.30, abs=0.05)
