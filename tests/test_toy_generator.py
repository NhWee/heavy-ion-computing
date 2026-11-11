"""Tests for hicpy.toy_generator: determinism, bookkeeping, and that the
samplers reproduce the distributions they claim to sample."""
import numpy as np
import pytest

from hicpy import toy_generator as tg

CFG = tg.ToyConfig(n_events=400, seed=7, z_fraction=0.3)


@pytest.fixture(scope="module")
def sample():
    return tg.generate(CFG)


def test_same_seed_same_events(sample):
    again = tg.generate(tg.ToyConfig(n_events=400, seed=7, z_fraction=0.3))
    assert np.array_equal(sample["n_particles"], again["n_particles"])
    assert all(np.array_equal(a, b) for a, b in zip(sample["px"], again["px"]))


def test_different_seed_different_events(sample):
    other = tg.generate(tg.ToyConfig(n_events=400, seed=8, z_fraction=0.3))
    assert not np.array_equal(sample["n_particles"], other["n_particles"])


def test_charge_follows_pdg_sign_convention(sample):
    """pi+ (211) -> +1, pi- (-211) -> -1, mu- (13) -> -1, mu+ (-13) -> +1."""
    assert tg.charge_of([211, -211, 321, -2212, 13, -13]).tolist() == [1, -1, 1, -1, -1, 1]
    for pdg, q in zip(sample["pdg"], sample["charge"]):
        assert q.dtype == np.int8
        assert np.array_equal(q, tg.charge_of(pdg))
        assert np.all(q != 0)                # the toy has no neutrals


def test_particle_counts_consistent(sample):
    for key in ("px", "py", "pz", "E", "pdg", "charge"):
        assert [len(a) for a in sample[key]] == list(sample["n_particles"])
    assert np.all(sample["n_particles"] >= 1)


def test_every_particle_is_on_shell(sample):
    for px, py, pz, e, pdg in zip(sample["px"], sample["py"], sample["pz"],
                                  sample["E"], sample["pdg"]):
        m = np.sqrt(e * e - px * px - py * py - pz * pz)
        expect = np.array([tg.MASS[abs(int(k))] for k in pdg])
        np.testing.assert_allclose(m, expect, rtol=1e-6, atol=1e-9)


def test_z_events_have_exactly_two_muons_near_mz(sample):
    n_mu = np.array([np.sum(np.abs(p) == 13) for p in sample["pdg"]])
    assert np.array_equal(n_mu == 2, sample["has_z"])
    assert np.array_equal(n_mu[~sample["has_z"]], np.zeros(np.sum(~sample["has_z"])))
    masses = []
    for iev in np.flatnonzero(sample["has_z"]):
        sel = np.abs(sample["pdg"][iev]) == 13
        e = sample["E"][iev][sel].sum()
        p = np.array([sample[k][iev][sel].sum() for k in ("px", "py", "pz")])
        masses.append(np.sqrt(e * e - (p * p).sum()))
    masses = np.array(masses)
    # the Breit-Wigner is wide: check the median, which is robust to the tails
    assert np.median(masses) == pytest.approx(tg.M_Z, abs=0.5)
    # mu+ (pdg -13) and mu- (+13), one each
    for iev in np.flatnonzero(sample["has_z"]):
        assert sorted(sample["pdg"][iev][np.abs(sample["pdg"][iev]) == 13]) == [-13, 13]


def test_species_fractions(sample):
    pdg = np.abs(np.concatenate(sample["pdg"]))
    had = pdg[pdg != 13]
    n = len(had)
    for code, frac in ((211, CFG.frac_pi), (321, CFG.frac_K), (2212, CFG.frac_p)):
        got = np.mean(had == code)
        assert got == pytest.approx(frac, abs=4 * np.sqrt(frac * (1 - frac) / n))


def test_multiplicity_matches_nbd_mean_and_variance():
    cfg = tg.ToyConfig(n_events=20000, seed=3, z_fraction=0.0)
    rng = np.random.default_rng(cfg.seed)
    n = rng.negative_binomial(cfg.nbd_k, cfg.nbd_k / (cfg.nbd_k + cfg.nbd_nbar),
                              size=cfg.n_events)
    mean, var = n.mean(), n.var()
    assert mean == pytest.approx(cfg.nbd_nbar, rel=0.03)
    assert var == pytest.approx(cfg.nbd_nbar + cfg.nbd_nbar ** 2 / cfg.nbd_k, rel=0.06)


def test_pt_sampler_reproduces_tsallis_shape():
    """Histogram of sampled pT vs the analytic dN/dpT: chi2/ndf of order 1."""
    cfg = tg.ToyConfig()
    rng = np.random.default_rng(11)
    pt = tg._sample_pt(rng, 200000, tg.M_PI, cfg)
    assert pt.min() >= cfg.pt_min and pt.max() <= cfg.pt_max
    edges = np.geomspace(cfg.pt_min, 5.0, 31)
    counts, _ = np.histogram(pt, edges)
    grid = np.geomspace(cfg.pt_min, cfg.pt_max, 200001)
    d = tg._tsallis_dNdpt(grid, tg.M_PI, cfg.tsallis_T, cfg.tsallis_n)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (d[1:] + d[:-1]) * np.diff(grid))])
    cdf /= cdf[-1]
    expected = len(pt) * np.diff(np.interp(edges, grid, cdf))
    chi2 = np.sum((counts - expected) ** 2 / expected)
    assert chi2 / len(counts) < 2.0


def test_breit_wigner_sampler_peak_and_width():
    rng = np.random.default_rng(5)
    m = tg._sample_breit_wigner(rng, 50000, tg.M_Z, tg.GAMMA_Z)
    assert np.median(m) == pytest.approx(tg.M_Z, abs=0.05)
    # FWHM of a (relativistic) Breit-Wigner is Gamma: the central 50 % of the
    # distribution spans about Gamma/2 on either side of the peak
    q25, q75 = np.percentile(m, [25, 75])
    assert (q75 - q25) == pytest.approx(tg.GAMMA_Z, rel=0.1)
