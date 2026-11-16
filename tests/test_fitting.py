"""Tests for hicpy.fitting and hicpy.histograms: the integration rule is exact,
the Poisson fit recovers known parameters, and the samplers the closure test
relies on reproduce their densities inside the analysis windows."""
import numpy as np
import pytest

from hicpy import fitting as fitm
from hicpy import histograms as hh
from hicpy import toy_generator as tg

T0, N0 = 0.135, 6.8


def tsallis_cdf(pt, C, T, n, m0):
    """Closed form of the integral of 2 pi pT f_inv dpT (pT dpT = mT dmT)."""
    u = 1.0 + (np.sqrt(pt * pt + m0 * m0) - m0) / (n * T)
    return 2 * np.pi * C * n * T * ((m0 - n * T) * u ** (1 - n) / (1 - n)
                                    + n * T * u ** (2 - n) / (2 - n))


@pytest.mark.parametrize("m0", [tg.M_PI, tg.M_K, tg.M_P])
def test_bin_integrals_match_closed_form(m0):
    edges = hh.make_edges({"type": "log", "n_bins": 40, "lo": 0.15, "hi": 50.0})
    got = fitm.bin_integrals(fitm.tsallis_dndpt, edges, 1.0, T0, N0, m0)
    exact = tsallis_cdf(edges[1:], 1.0, T0, N0, m0) - tsallis_cdf(edges[:-1], 1.0, T0, N0, m0)
    np.testing.assert_allclose(got, exact, rtol=1e-10)


def test_bin_centre_evaluation_is_biased_but_bin_average_is_not():
    """The size of the approximation the chi2 cross-check fit makes: on the
    40 log bins the bin-centre value under-represents the bin average by up to
    ~4.6 % at high pT (troubleshooting #16)."""
    edges = hh.make_edges({"type": "log", "n_bins": 40, "lo": 0.15, "hi": 50.0})
    ctr, w = hh.centers_widths(edges)
    avg = fitm.bin_integrals(fitm.tsallis_invariant, edges, 1.0, T0, N0) / w
    ratio = avg / fitm.tsallis_invariant(ctr, 1.0, T0, N0)
    assert ratio[0] == pytest.approx(1.0, abs=1e-3)
    assert 1.03 < ratio[-1] < 1.06


def test_sub_edges_requires_a_contiguous_selection():
    e = np.arange(6.0)
    np.testing.assert_array_equal(fitm.sub_edges(e, np.array([0, 1, 1, 1, 0], bool)),
                                  [1.0, 2.0, 3.0, 4.0])
    with pytest.raises(ValueError):
        fitm.sub_edges(e, np.array([1, 0, 1, 0, 0], bool))


def test_poisson_deviance_is_zero_for_a_perfect_model_and_handles_empty_bins():
    n = np.array([0.0, 3.0, 10.0])
    assert fitm.poisson_deviance(n + (n == 0) * 1e-12, n) == pytest.approx(0.0, abs=1e-9)
    assert fitm.poisson_deviance([2.0], [0.0]) == pytest.approx(4.0)   # 2 mu


def test_poisson_fit_recovers_tsallis_parameters():
    cfg = tg.ToyConfig()
    rng = np.random.default_rng(11)
    pts = tg._sample_pt(rng, 200_000, tg.M_PI, cfg)
    edges = hh.make_edges({"type": "log", "n_bins": 40, "lo": 0.15, "hi": 50.0})
    k, _ = hh.hist_counts(pts, edges)
    fit = fitm.fit_poisson(fitm.tsallis_expected_counts(edges, 1, 1.0, tg.M_PI), k,
                           p0={"C": 5e4, "T": 0.13, "n": 7.0},
                           limits={"C": (0, None), "T": (0.01, 1.0), "n": (1.0, 50.0)})
    assert fit["converged"] and fit["accurate_covariance"]
    assert not fit["at_limit"]
    T, dT = fit["params"][1], fit["errors"][1]
    n, dn = fit["params"][2], fit["errors"][2]
    assert abs(T - cfg.tsallis_T) < 3 * dT
    assert abs(n - cfg.tsallis_n) < 3 * dn
    # T and n are strongly positively correlated in this form
    assert fit["correlation"][1][2] > 0.8
    # every scattered start converged to the same minimum
    assert fit["multistart"]["max_delta_deviance"] < 1e-3
    assert fit["multistart"]["max_param_shift_sigma"] < 0.05


def test_breit_wigner_sampler_matches_its_density_in_the_analysis_window():
    """_sample_breit_wigner normalises its rejection weights by the maximum of
    each batch instead of a global bound.  The weight f/g peaks at the upper
    truncation edge (~141 GeV), so only the far tail above ~130 GeV can be
    affected; inside 60-120 GeV the accepted shape must be exact."""
    rng = np.random.default_rng(5)
    m = tg._sample_breit_wigner(rng, 400_000, tg.M_Z, tg.GAMMA_Z)
    edges = np.linspace(60.0, 120.0, 61)
    k, _ = hh.hist_counts(m, edges)
    mu = fitm.bw_expected_counts(edges)(k.sum(), tg.M_Z, tg.GAMMA_Z)
    use = mu > 20
    pearson = np.sum((k[use] - mu[use]) ** 2 / mu[use])
    ndf = int(use.sum()) - 1
    assert pearson < ndf + 5 * np.sqrt(2 * ndf)


def test_poisson_fit_recovers_breit_wigner_parameters():
    rng = np.random.default_rng(6)
    m = tg._sample_breit_wigner(rng, 50_000, tg.M_Z, tg.GAMMA_Z)
    edges = np.linspace(60.0, 120.0, 61)
    k, _ = hh.hist_counts(m, edges)
    fit = fitm.fit_poisson(fitm.bw_expected_counts(edges), k,
                           p0={"N": k.sum(), "m0": 91.2, "gamma": 2.5},
                           limits={"N": (0, None), "m0": (80, 100), "gamma": (0.1, 20)})
    assert fit["converged"]
    assert abs(fit["params"][1] - tg.M_Z) < 3 * fit["errors"][1]
    assert abs(fit["params"][2] - tg.GAMMA_Z) < 3 * fit["errors"][2]
    assert fit["params"][0] == pytest.approx(k.sum(), rel=1e-3)   # yield = entries


def test_chi2_fit_reports_dropped_zero_error_points():
    x = np.linspace(1, 10, 10)
    y = 2.0 * x
    err = np.ones_like(x)
    err[[3, 7]] = 0.0
    fit = fitm.fit_chi2(lambda x, a: a * x, x, y, err, p0=[1.0])
    assert fit["converged"] and fit["n_points"] == 8 and fit["n_dropped_zero_error"] == 2
    assert fit["params"][0] == pytest.approx(2.0)
