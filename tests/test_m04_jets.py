"""hicpy.jets on controlled inputs (no PYTHIA needed): the three properties
Module 4's conclusions rest on."""
import numpy as np
import pytest

pytest.importorskip("fastjet")

from hicpy.jets import (pseudojets, area_definition, cluster, median_rho,  # noqa: E402
                        toy_background, delta_phi)

RHO_DENSITY, MEAN_PT = 30.0, 0.5            # dN/deta dphi, GeV
RHO_TRUE = RHO_DENSITY * MEAN_PT


def test_isolated_antikt_jet_is_a_circle():
    ad = area_definition(3.0, 0.01)
    jets, cs = cluster(pseudojets([100.0], [0.0], [0.0], [100.0]), "antikt", 0.4, 10.0, 2.0, ad)
    assert len(jets) == 1
    assert jets[0].area() == pytest.approx(np.pi * 0.4 ** 2, rel=0.05)


def test_median_rho_recovers_a_uniform_background():
    rng = np.random.default_rng(3)
    ad = area_definition(3.0, 0.01)
    rhos = [median_rho(pseudojets(*toy_background(rng, RHO_DENSITY, MEAN_PT, 3.0)), ad)[0]
            for _ in range(10)]
    assert np.mean(rhos) == pytest.approx(RHO_TRUE, rel=0.08)


def test_embedded_jet_shifts_by_rho_times_area_and_subtraction_removes_it():
    rng = np.random.default_rng(4)
    ad = area_definition(3.0, 0.01)
    raw, sub, areas = [], [], []
    for _ in range(30):
        sig = pseudojets([100.0], [0.0], [0.0], [100.0])
        full = sig + pseudojets(*toy_background(rng, RHO_DENSITY, MEAN_PT, 3.0))
        jets, cs = cluster(full, "antikt", 0.4, 10.0, 2.4, ad)
        j = jets[0]
        rho, _ = median_rho(full, ad)
        raw.append(j.pt() - 100.0); sub.append(j.pt() - rho * j.area() - 100.0); areas.append(j.area())
    fluct = np.sqrt(np.pi * 0.4 ** 2 * RHO_DENSITY * 2 * MEAN_PT ** 2)     # ~2.7 GeV per jet
    tol = 3 * fluct / np.sqrt(len(raw))
    assert np.mean(raw) == pytest.approx(RHO_TRUE * np.mean(areas), abs=tol)
    assert np.mean(sub) == pytest.approx(0.0, abs=tol)


def test_kt_jets_are_not_circles():
    rng = np.random.default_rng(5)
    ad = area_definition(3.0, 0.01)
    full = pseudojets([100.0], [0.0], [0.0], [100.0]) + pseudojets(*toy_background(rng, 5.0, 0.5, 3.0))
    # keep each ClusterSequence alive while its jets are used (troubleshooting #20)
    jets_akt, cs_akt = cluster(full, "antikt", 0.4, 10.0, 2.0, ad)
    jets_kt, cs_kt = cluster(full, "kt", 0.4, 10.0, 2.0, ad)
    a_akt, a_kt = jets_akt[0].area(), jets_kt[0].area()
    assert abs(a_akt / (np.pi * 0.16) - 1) < 0.1
    assert abs(a_kt / (np.pi * 0.16) - 1) > abs(a_akt / (np.pi * 0.16) - 1)


def test_delta_phi_folds_into_0_pi():
    d = delta_phi(np.array([0.1, 6.2, 3.0]), np.array([6.2, 0.1, -3.0]))
    assert np.all((d >= 0) & (d <= np.pi))
    assert d[0] == pytest.approx(0.1 + 2 * np.pi - 6.2)
