"""Fit statistics and fit models shared by the analysis scripts.

Two statistics, each for its own job (docs/troubleshooting.md #15, #16):

fit_poisson -- binned maximum likelihood on RAW COUNTS.  THE PHYSICS FIT.
    The model is integrated over every bin (Gauss-Legendre, exact for the
    smooth models used here), so there is no bin-centre approximation, and
    empty or nearly empty bins enter with their correct Poisson weight.  The
    minimised quantity is the Baker-Cousins likelihood-ratio deviance

        D = 2 * sum_i [ mu_i - n_i + n_i ln(n_i / mu_i) ]  =  -2 ln(L / L_sat),

    whose minimum is also a goodness-of-fit statistic: asymptotically chi2
    with ndf = bins - free parameters, but only where the expected counts are
    not small.  analysis/m01_fit_validation.py measures its real distribution
    with pseudo-experiments.

fit_chi2 -- least squares at bin centres (scipy.optimize.curve_fit) with the
    uncertainties taken from the data.  This is what ROOT's default
    TH1::Fit does, so it is implemented identically in all three Module 1
    arms and is the right tool for the CROSS-CHECK between implementations.
    It is not used for physics numbers: on the log-binned, steeply falling
    pT spectrum the bin-centre evaluation biases T by about +1 sigma, and on
    the 240-event Z peak the data-derived sigmas bias the width low and
    underestimate its error (both measured in m01_fit_validation.py).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit

M_PI = 0.13957039                     # GeV, default mass in the Tsallis form


# ---------------------------------------------------------------- models
def tsallis_invariant(pt, C, T, n, m0=M_PI):
    """Invariant yield 1/(2 pi pT) d2N/(dpT deta) = C [1 + (mT - m0)/(n T)]^(-n)."""
    mt = np.sqrt(pt * pt + m0 * m0)
    return C * (1.0 + (mt - m0) / (n * T)) ** (-n)


def tsallis_dndpt(pt, C, T, n, m0=M_PI):
    """dN/dpT per unit eta and per event: 2 pi pT times the invariant yield."""
    return 2.0 * np.pi * pt * tsallis_invariant(pt, C, T, n, m0)


def rel_breit_wigner(m, m0, gamma):
    """Relativistic Breit-Wigner shape m^2 G^2 / ((m^2 - m0^2)^2 + m0^2 G^2).

    Unnormalised (= 1 at m = m0).  It is the SAME form the toy generator
    samples (hicpy.toy_generator._sample_breit_wigner), so a fit with it must
    close on the generated masses.
    """
    num = m ** 2 * gamma ** 2
    den = (m ** 2 - m0 ** 2) ** 2 + m0 ** 2 * gamma ** 2
    return num / den


def bw_plus_linear(m, A, m0, gamma, a, b):
    """Breit-Wigner peak on a linear background: the Session 1-3 chi2 model,
    kept only for the method comparison in m01_fit_validation.py."""
    return A * rel_breit_wigner(m, m0, gamma) + a + b * m


# ---------------------------------------------------------------- integration
_GL_RULES: dict[int, tuple[np.ndarray, np.ndarray]] = {}


def bin_integrals(func, edges, *params, n_nodes=8):
    """Integral of func(x, *params) over every bin [e_i, e_{i+1}].

    Gauss-Legendre with n_nodes points per bin is exact for polynomials of
    degree 2 n_nodes - 1; for the smooth models here it is accurate to ~1e-12
    on bins a few to a few tens of percent wide (tests/test_fitting.py checks
    it against the closed-form Tsallis integral).
    """
    if n_nodes not in _GL_RULES:
        _GL_RULES[n_nodes] = np.polynomial.legendre.leggauss(n_nodes)
    xg, wg = _GL_RULES[n_nodes]
    edges = np.asarray(edges, dtype=np.float64)
    half = 0.5 * np.diff(edges)
    mid = 0.5 * (edges[1:] + edges[:-1])
    vals = func(mid[:, None] + half[:, None] * xg[None, :], *params)
    return half * np.sum(wg[None, :] * vals, axis=1)


def tsallis_expected_counts(edges, n_events, delta_eta, m0=M_PI):
    """Expected RAW counts per pT bin for the invariant-yield parameters
    (C, T, n): N_ev * delta_eta * integral over the bin of 2 pi pT f_inv."""
    def expected(C, T, n):
        return n_events * delta_eta * bin_integrals(tsallis_dndpt, edges, C, T, n, m0)
    return expected


def bw_expected_counts(edges):
    """Expected counts per bin for N events in the window [edges[0], edges[-1]]
    distributed as rel_breit_wigner (normalised over that window)."""
    def expected(N, m0, gamma):
        s = bin_integrals(rel_breit_wigner, edges, m0, gamma)
        return N * s / s.sum()
    return expected


def sub_edges(edges, mask):
    """Edges of the contiguous run of bins selected by `mask`."""
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        raise ValueError("no bin selected")
    if np.any(np.diff(idx) != 1):
        raise ValueError("selected bins are not contiguous")
    return np.asarray(edges, dtype=np.float64)[idx[0]: idx[-1] + 2]


# ---------------------------------------------------------------- statistics
def poisson_deviance(mu, n):
    """Baker-Cousins deviance 2 sum[mu - n + n ln(n/mu)]; empty bins give 2 mu."""
    mu = np.maximum(np.asarray(mu, dtype=np.float64), 1e-300)
    n = np.asarray(n, dtype=np.float64)
    pos = n > 0
    return 2.0 * (np.sum(mu - n) + np.sum(n[pos] * np.log(n[pos] / mu[pos])))


def _inside(x, lo, hi):
    """Move a start value strictly inside its limits (Minuit needs that)."""
    if lo is not None and hi is not None:
        span = hi - lo
        return min(max(x, lo + 1e-3 * span), hi - 1e-3 * span)
    if lo is not None:
        return max(x, lo + 1e-3 * max(abs(lo), 1.0))
    if hi is not None:
        return min(x, hi - 1e-3 * max(abs(hi), 1.0))
    return x


def fit_poisson(expected, counts, p0, limits=None, label="", function="",
                n_starts=5, seed=0, tol=1e-6):
    """Binned Poisson maximum-likelihood fit (minimises poisson_deviance).

    expected : callable(*params) -> expected counts per bin
    counts   : observed raw counts per bin (same length)
    p0       : dict name -> start value; its order is the parameter order
    limits   : dict name -> (lo, hi), None for an open side.  Use bounds only
               where they are physical (a yield >= 0, a width > 0); a result
               sitting on a bound is listed in "at_limit".
    n_starts : multi-start check: n_starts - 1 extra fits from start values
               scattered by +-30 % around p0.  If one of them reaches a
               deviance lower by more than 1e-3 (a genuinely different
               minimum, not convergence noise) it becomes the result
               ("multistart.improved").
    tol      : MIGRAD tolerance; the EDM goal is 0.002 * tol, so the default
               1e-6 locates the minimum to ~1e-4 sigma.  The iminuit default
               (0.1) leaves ~1e-2 sigma of start-point dependence, enough to
               make results differ between machines.

    Returns a record with params, errors, covariance, correlation, deviance,
    ndf, the multi-start summary and convergence flags.
    """
    from iminuit import Minuit

    counts = np.asarray(counts, dtype=np.float64)
    names = list(p0)
    limits = dict(limits or {})

    def cost(*p):
        return poisson_deviance(expected(*p), counts)

    def run(start):
        m = Minuit(cost, *start, name=names)
        m.errordef = Minuit.LEAST_SQUARES          # D = -2 ln L + const
        m.tol = tol
        for k, (lo, hi) in limits.items():
            m.limits[k] = (lo, hi)
        m.migrad()
        if m.valid:
            m.hesse()
        return m

    start0 = [_inside(float(p0[k]), *limits.get(k, (None, None))) for k in names]
    best = run(start0)
    rng = np.random.default_rng(seed)
    others = []
    for _ in range(max(n_starts - 1, 0)):
        st = [_inside(v * (1.0 + 0.3 * rng.uniform(-1.0, 1.0)), *limits.get(k, (None, None)))
              for k, v in zip(names, start0)]
        others.append(run(st))
    improved = False
    for m in others:
        if m.valid and (not best.valid or m.fval < best.fval - 1e-3):
            best, improved = m, True
    if not best.valid:
        return {"label": label, "function": function, "converged": False,
                "method": "binned Poisson likelihood", "error": str(best.fmin)}

    vals = np.array(best.values, dtype=float)
    errs = np.array(best.errors, dtype=float)
    cov = np.array(best.covariance, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov / np.outer(errs, errs)
    conv = [m for m in others if m.valid]
    d_dev = [float(m.fval - best.fval) for m in conv]
    shifts = [float(np.max(np.abs(np.array(m.values) - vals) / np.where(errs > 0, errs, 1.0)))
              for m in conv]
    ndf = int(counts.size - len(names))

    def near(v, bound, err):
        return bound is not None and abs(v - bound) <= max(1e-9 * max(1.0, abs(bound)),
                                                             0.01 * err)
    at_limit = [k for k, v, e in zip(names, vals, errs)
                if any(near(v, b, e) for b in limits.get(k, (None, None)))]
    return {
        "label": label,
        "function": function,
        "method": ("binned Poisson likelihood (Baker-Cousins deviance), "
                   "model integrated over each bin"),
        "converged": True,
        "accurate_covariance": bool(best.accurate),
        "at_limit": at_limit,
        "n_bins": int(counts.size),
        "n_entries": float(counts.sum()),
        "names": names,
        "p0": [float(p0[k]) for k in names],
        "params": vals.tolist(),
        "errors": errs.tolist(),
        "covariance": cov.tolist(),
        "correlation": corr.tolist(),
        "deviance": float(best.fval),
        "ndf": ndf,
        "deviance_ndf": float(best.fval) / ndf if ndf > 0 else None,
        "multistart": {"n_starts": int(n_starts), "n_converged": len(conv) + 1,
                       "improved": improved,
                       "max_delta_deviance": max(d_dev, default=0.0),
                       "max_param_shift_sigma": max(shifts, default=0.0)},
    }


def fit_chi2(func, x, y, yerr, p0, bounds=(-np.inf, np.inf), label=""):
    """Least-squares fit at the given x (bin centres) with data-derived sigmas.

    Points with zero or non-finite uncertainty are dropped (and counted, so a
    silently shrinking fit range is visible in "n_points").  Cross-check use
    only -- see the module docstring.
    """
    x, y, yerr = (np.asarray(a, dtype=np.float64) for a in (x, y, yerr))
    mask = (yerr > 0) & np.isfinite(y)
    try:
        popt, pcov = curve_fit(func, x[mask], y[mask], p0=p0,
                               sigma=yerr[mask], absolute_sigma=True,
                               bounds=bounds, maxfev=40000)
    except Exception as exc:                       # noqa: BLE001
        return {"label": label, "converged": False, "error": str(exc)}
    resid = (y[mask] - func(x[mask], *popt)) / yerr[mask]
    chi2 = float(np.sum(resid ** 2))
    ndf = int(mask.sum() - len(popt))
    return {
        "label": label,
        "function": func.__name__,
        "method": "least squares at bin centres (curve_fit), sigma from data",
        "converged": True,
        "fit_range": [float(x[mask].min()), float(x[mask].max())],
        "n_points": int(mask.sum()),
        "n_dropped_zero_error": int((~mask).sum()),
        "p0": [float(v) for v in p0],
        "params": [float(v) for v in popt],
        "errors": [float(v) for v in np.sqrt(np.diag(pcov))],
        "covariance": [[float(v) for v in row] for row in pcov],
        "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf if ndf > 0 else None,
    }


# Backwards-compatible name used by the Session 1-3 scripts.
do_fit = fit_chi2

__all__ = ["tsallis_invariant", "tsallis_dndpt", "rel_breit_wigner",
           "bw_plus_linear", "bin_integrals", "tsallis_expected_counts",
           "bw_expected_counts", "sub_edges", "poisson_deviance", "fit_poisson",
           "fit_chi2", "do_fit", "M_PI"]
