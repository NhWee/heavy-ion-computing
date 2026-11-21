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


def test_dijet_selection_and_aj():
    import pathlib
    import sys
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "analysis"))
    import m09_dijet as dj
    jets = np.zeros((4, dj.N_JETS, 3))
    jets[0, :2] = [[150, 0, 0], [100, 0, np.pi]]          # balanced enough: A_J = 0.2
    jets[1, :2] = [[100, 0, 0], [150, 0, np.pi]]          # order swapped: re-sorted, same A_J
    jets[2, :2] = [[150, 0, 0], [40, 0, np.pi]]           # subleading below 50: rejected
    jets[3, :2] = [[150, 0, 0], [100, 0, 1.0]]            # not back-to-back: rejected
    aj, ok = dj.select_aj(jets)
    assert list(ok) == [True, True, False, False]
    assert aj == pytest.approx([0.2, 0.2])


def test_background_sigma_formula_and_ordering():
    """Poisson background width: the Module 4 case (dN/detadphi 30, <pT> 0.5, R 0.4) gives 2.75 GeV,
    and the Angantyr-based widths fall from central to peripheral."""
    import pathlib
    import sys
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "analysis"))
    import m09_dijet as dj
    assert np.sqrt(np.pi * 0.4 ** 2 * 30 * 2 * 0.5 ** 2) == pytest.approx(2.746, abs=0.01)
    s = dj.background_sigma()
    assert len(s) == 9 and np.all(np.diff(s) < 0) and 8 < s[0] < 15
