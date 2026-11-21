"""Module 9 geometry: path integrals and the participant plane."""
import numpy as np
import pytest

from hicpy.glauber import participant_plane, path_integrals


def _numeric(x0, phi, src, w, zeta):
    t = np.linspace(0, 60, 300001)
    pts = x0 + t[:, None] * np.array([np.cos(phi), np.sin(phi)])
    rho = np.sum(np.exp(-np.sum((pts[:, None, :] - src[None]) ** 2, axis=2) / (2 * w * w)), axis=1) / (2 * np.pi * w * w)
    return np.trapezoid(rho * t ** zeta, t)


@pytest.mark.parametrize("zeta", [0, 1])
def test_path_integral_closed_form(zeta):
    rng = np.random.default_rng(4)
    src = rng.normal(0, 2.5, (60, 2)); x0 = np.array([0.3, 0.8]); phi = 2.1
    got = path_integrals(x0[None], np.array([phi]), src, 0.5, zeta)[0]
    assert got == pytest.approx(_numeric(x0, phi, src, 0.5, zeta), rel=1e-6)


def test_path_length_scaling_in_a_uniform_slab():
    """Uniform density strip of length L ahead of the start: I_0 ~ L, I_1 ~ L^2 / 2."""
    xs = np.arange(0.05, 8.0, 0.1); ys = np.arange(-6, 6.01, 0.1)
    src = np.array([(x, y) for x in xs for y in ys])                   # 100 sources / fm^2
    i0 = path_integrals(np.array([[0.0, 0.0]]), np.array([0.0]), src, 0.3, 0)[0]
    i1 = path_integrals(np.array([[0.0, 0.0]]), np.array([0.0]), src, 0.3, 1)[0]
    assert i0 == pytest.approx(100 * 8.0, rel=0.03)
    assert i1 == pytest.approx(100 * 8.0 ** 2 / 2, rel=0.03)


def test_participant_plane_points_along_the_short_axis():
    rng = np.random.default_rng(1)
    pts = np.column_stack([rng.normal(0, 2.0, 5000), rng.normal(0, 4.0, 5000)])   # narrow in x
    rot = 0.4
    R = np.array([[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]])
    psi, eps = participant_plane(pts @ R.T)
    assert eps == pytest.approx((16 - 4) / (16 + 4), abs=0.03)
    d = (psi - rot + np.pi / 2) % np.pi - np.pi / 2                       # Psi_2 is defined modulo pi
    assert abs(d) < 0.05
