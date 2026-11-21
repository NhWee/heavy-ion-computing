"""Dimuon invariant mass and resonance fits for Module 8.

THE MEASUREMENT IN ONE PARAGRAPH
--------------------------------
Two muons from the decay of a particle of mass M satisfy
    M^2 = (E1 + E2)^2 - |p1 + p2|^2.
Histogramming M for all opposite-sign muon pairs shows the resonances that
decay to mu+ mu- (omega/phi, J/psi, psi(2S), the Upsilon(1S, 2S, 3S) family,
the Z) on top of a smooth continuum (Drell-Yan, heavy-flavour decays, random
pairs).  The peak positions measure the masses (and calibrate the momentum
scale), the peak widths measure the detector's momentum resolution (the
J/psi and Upsilon natural widths are tiny), and the yields per state give
the excited-to-ground state ratios that heavy-ion physics compares between
pp and Pb-Pb (sequential suppression, Module 9).

FIT MODEL
---------
- quarkonia: Crystal Ball function (a Gaussian core with a power-law tail on
  the LOW-mass side from final-state photon radiation), one per state.
  States of one family share alpha and n; their widths scale with mass
  (resolution is a fixed fraction of the momentum): sigma_i = sigma_1 m_i/m_1.
- Z: Breit-Wigner (Gamma_Z fixed to PDG) convolved with a Gaussian
  resolution: the Voigt profile.
- masses: one free MASS SCALE per family, m_i = scale x m_i(PDG): the
  splittings are fixed to the PDG values (as CMS does), which keeps the fit
  from swapping the three Upsilon states, and scale - 1 is the momentum-scale
  offset of the detector.
- background: exponential times (1 + c1 u + c2 u^2), u = position in the
  window from 0 to 1.
Binned Poisson likelihood, minimised with scipy.
"""
from __future__ import annotations

import numpy as np
from scipy import optimize, special

MU_MASS = 0.1056584
# PDG values (GeV), Review of Particle Physics
PDG = {"J/psi": 3.096900, "psi(2S)": 3.686097, "Upsilon(1S)": 9.46040, "Upsilon(2S)": 10.02340,
       "Upsilon(3S)": 10.35520, "Z": 91.1876}
GAMMA_Z = 2.4952


def pair_mass(pt1, eta1, phi1, pt2, eta2, phi2, m=MU_MASS):
    """Invariant mass of two particles given (pT, eta, phi) and a common mass."""
    px1, py1, pz1 = pt1 * np.cos(phi1), pt1 * np.sin(phi1), pt1 * np.sinh(eta1)
    px2, py2, pz2 = pt2 * np.cos(phi2), pt2 * np.sin(phi2), pt2 * np.sinh(eta2)
    e1 = np.sqrt(px1 ** 2 + py1 ** 2 + pz1 ** 2 + m * m)
    e2 = np.sqrt(px2 ** 2 + py2 ** 2 + pz2 ** 2 + m * m)
    m2 = (e1 + e2) ** 2 - (px1 + px2) ** 2 - (py1 + py2) ** 2 - (pz1 + pz2) ** 2
    return np.sqrt(np.clip(m2, 0, None))


# ----------------------------------------------------------------------------- shapes
def crystal_ball(x, mu, sigma, alpha, n):
    """Unit-normalised Crystal Ball density (tail on the low side)."""
    alpha = abs(alpha)
    t = (x - mu) / sigma
    a = (n / alpha) ** n * np.exp(-0.5 * alpha ** 2)
    b = n / alpha - alpha
    core = np.exp(-0.5 * t ** 2)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        tail = a * np.power(np.where(b - t > 0, b - t, np.inf), -n)
    f = np.where(t > -alpha, core, tail)
    c = n / alpha / (n - 1) * np.exp(-0.5 * alpha ** 2)
    d = np.sqrt(np.pi / 2) * (1 + special.erf(alpha / np.sqrt(2)))
    return f / (sigma * (c + d))


def voigt(x, mu, sigma, gamma=GAMMA_Z):
    """Breit-Wigner (full width gamma) convolved with a Gaussian (sigma)."""
    return special.voigt_profile(x - mu, sigma, gamma / 2)


def background(x, x0, x1, lam, c1, c2):
    u = (x - x0) / (x1 - x0)
    return np.exp(-np.clip(lam, -50, 50) * u) * np.clip(1 + c1 * u + c2 * u * u, 1e-12, None)


# ----------------------------------------------------------------------------- fit
def _bin_integral(f, edges, n_sub=5):
    """Integral of a density over each bin (Simpson-like sub-sampling)."""
    sub = np.linspace(0, 1, n_sub + 1)
    x = edges[:-1, None] + np.diff(edges)[:, None] * (0.5 * (sub[1:] + sub[:-1]))[None, :]
    return f(x).mean(axis=1) * np.diff(edges)


def fit_window(edges, counts, states, model="cb", fix=None, start=None):
    """Binned Poisson maximum-likelihood fit of a mass window.

    states: list of PDG names in this window (first = reference for the width
    scaling).  model 'cb' (quarkonia) or 'voigt' (Z).  Returns a dict with the
    fitted masses, width, yields (with errors from the inverse Hessian), the
    background yield under each peak (+-2 sigma) and the fit quality.
    """
    edges = np.asarray(edges, float); counts = np.asarray(counts, float)
    x0, x1 = edges[0], edges[-1]; n_tot = counts.sum()
    ref = PDG[states[0]]
    names = (["scale", "sigma"] + (["alpha", "n"] if model == "cb" else [])
             + [f"N_{s}" for s in states] + ["N_bkg", "lam", "c1", "c2"])
    p0 = ([1.0, 0.01 * ref if model == "cb" else 1.5] + ([1.5, 3.0] if model == "cb" else [])
          + [0.1 * n_tot / len(states)] * len(states) + [0.5 * n_tot, 1.0, 0.0, 0.0])
    if start:
        p0 = [start.get(k, v) for k, v in zip(names, p0)]
    fix = fix or {}
    free = [i for i, k in enumerate(names) if k not in fix]

    def unpack(q):
        p = np.array(p0, float)
        p[free] = q
        for k, v in fix.items():
            p[names.index(k)] = v
        return dict(zip(names, p))

    def expected(par):
        mu = np.zeros(len(counts))
        for s in states:
            m = par["scale"] * PDG[s]; sig = par["sigma"] * PDG[s] / ref
            if model == "cb":
                f = lambda x, m=m, sig=sig: crystal_ball(x, m, sig, par["alpha"], par["n"])  # noqa: E731
            else:
                f = lambda x, m=m, sig=sig: voigt(x, m, sig)                                    # noqa: E731
            mu += par[f"N_{s}"] * _bin_integral(f, edges)
        bk = _bin_integral(lambda x: background(x, x0, x1, par["lam"], par["c1"], par["c2"]), edges)
        mu += par["N_bkg"] * bk / max(bk.sum(), 1e-300)
        return mu

    def nll(q):
        par = unpack(q)
        if par["sigma"] <= 0 or not 0.9 < par["scale"] < 1.1 or (model == "cb" and (par["n"] <= 1.05 or abs(par["alpha"]) < 0.1)):
            return 1e30
        if any(par[f"N_{s}"] < 0 for s in states) or par["N_bkg"] < 0:
            return 1e30
        if max(abs(par["lam"]), abs(par["c1"]), abs(par["c2"])) > 30:      # keep the background shape finite
            return 1e30
        mu = np.clip(expected(par), 1e-12, None)
        return float(np.sum(mu - counts * np.log(mu)))

    q0 = np.array(p0, float)[free]
    res = optimize.minimize(nll, q0, method="Nelder-Mead",
                            options={"maxiter": 40000, "maxfev": 40000, "xatol": 1e-7, "fatol": 1e-6, "adaptive": True})
    res = optimize.minimize(nll, res.x, method="Powell", options={"maxiter": 40000, "xtol": 1e-7, "ftol": 1e-9})
    par = unpack(res.x)
    # Errors: numerical Hessian of the NLL in the free parameters EXCEPT the
    # background shape (lam, c1, c2), which stays at its best fit.  Where the data
    # barely constrain the shape (no continuum under the Z) the full Hessian is
    # singular and its finite differences can step into the parameter walls.
    shape = {"lam", "c1", "c2"}
    err_idx = [k for k, i in enumerate(free) if names[i] not in shape]
    best = np.array(res.x, float)

    def nll_sub(q):
        z = best.copy(); z[err_idx] = q
        return nll(z)

    h = _hessian(nll_sub, best[err_idx])
    try:
        cov = np.linalg.inv(h)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(h)
    err_sub = np.sqrt(np.clip(np.diag(cov), 0, None))
    err = dict.fromkeys(names, 0.0)
    for k, e in zip(err_idx, err_sub):
        err[names[free[k]]] = float(e)
    mu = expected(par)
    ok = counts > 0
    chi2 = float(2 * np.sum(mu - counts + np.where(ok, counts * np.log(np.where(ok, counts, 1) / mu), 0)))
    out = {"model": model, "states": states, "params": {k: float(v) for k, v in par.items()},
           "errors": err, "chi2_poisson": chi2, "ndf": int(len(counts) - len(free)), "converged": bool(res.success),
           "expected": mu.tolist()}
    for s in states:
        sig = par["sigma"] * PDG[s] / ref
        out[f"m_{s}"] = float(par["scale"] * PDG[s]); out[f"m_{s}_err"] = float(err["scale"] * PDG[s])
        lo, hi = out[f"m_{s}"] - 2 * sig, out[f"m_{s}"] + 2 * sig
        bk = _bin_integral(lambda x: background(x, x0, x1, par["lam"], par["c1"], par["c2"]), edges)
        bk = par["N_bkg"] * bk / bk.sum()
        c = 0.5 * (edges[1:] + edges[:-1]); sel = (c > lo) & (c < hi)
        out[f"S/B_{s}"] = float(par[f"N_{s}"] * 0.95 / max(bk[sel].sum(), 1e-12))
    return out


def _hessian(f, x, rel=1e-4):
    x = np.asarray(x, float); n = len(x)
    h = np.maximum(np.abs(x) * rel, 1e-6)
    H = np.zeros((n, n)); f0 = f(x)
    for i in range(n):
        for j in range(i, n):
            if i == j:
                xp, xm = x.copy(), x.copy(); xp[i] += h[i]; xm[i] -= h[i]
                H[i, i] = (f(xp) - 2 * f0 + f(xm)) / h[i] ** 2
            else:
                xpp, xpm, xmp, xmm = (x.copy() for _ in range(4))
                xpp[i] += h[i]; xpp[j] += h[j]; xpm[i] += h[i]; xpm[j] -= h[j]
                xmp[i] -= h[i]; xmp[j] += h[j]; xmm[i] -= h[i]; xmm[j] -= h[j]
                H[i, j] = H[j, i] = (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4 * h[i] * h[j])
    return H


# ----------------------------------------------------------------------------- toy data
def toy_events(rng, n_events, composition, sigma_rel=0.01, bkg_slope=0.4):
    """Toy dimuon events with NanoAOD-like branches, for closure tests.

    composition: {state: fraction}, the rest is an exponential continuum in
    1-120 GeV.  Each resonance is drawn from a Gaussian of relative width
    sigma_rel (Z: Breit-Wigner x Gaussian), produced with a random pT and
    rapidity and decayed isotropically into mu+ mu-.  Returns flat numpy
    arrays per muon and the event offsets, i.e. what the analysis reads.
    """
    states = list(composition)
    probs = np.array([composition[s] for s in states] + [1 - sum(composition.values())])
    which = rng.choice(len(probs), size=n_events, p=probs)
    m = np.empty(n_events)
    for k, s in enumerate(states):
        sel = which == k
        if s == "Z":
            m[sel] = (PDG[s] + 0.5 * GAMMA_Z * np.tan(np.pi * (rng.random(sel.sum()) - 0.5))
                      + rng.normal(0, sigma_rel * PDG[s], sel.sum()))
        else:
            m[sel] = rng.normal(PDG[s], sigma_rel * PDG[s], sel.sum())
    bk = which == len(states)
    m[bk] = 1.0 + rng.exponential(1 / bkg_slope, bk.sum())
    m = np.clip(m, 2.5 * MU_MASS, None)
    pt = rng.exponential(8.0, n_events); y = rng.uniform(-2, 2, n_events); ph = rng.uniform(-np.pi, np.pi, n_events)
    mt = np.sqrt(m ** 2 + pt ** 2)
    P = np.stack([mt * np.cosh(y), pt * np.cos(ph), pt * np.sin(ph), mt * np.sinh(y)], axis=1)
    # decay in the rest frame
    pstar = np.sqrt(np.clip(m ** 2 / 4 - MU_MASS ** 2, 0, None))
    ct = rng.uniform(-1, 1, n_events); st = np.sqrt(1 - ct ** 2); fp = rng.uniform(0, 2 * np.pi, n_events)
    d = np.stack([np.full(n_events, 0.0), pstar * st * np.cos(fp), pstar * st * np.sin(fp), pstar * ct], axis=1)
    d[:, 0] = np.sqrt(pstar ** 2 + MU_MASS ** 2)
    mus = []
    for sign in (+1, -1):
        q = d.copy(); q[:, 1:] *= sign
        mus.append(_boost(q, P / m[:, None]))
    out = {}
    for name, f in (("pt", lambda v: np.hypot(v[:, 1], v[:, 2])),
                    ("eta", lambda v: np.arcsinh(v[:, 3] / np.hypot(v[:, 1], v[:, 2]))),
                    ("phi", lambda v: np.arctan2(v[:, 2], v[:, 1]))):
        out[name] = np.stack([f(mus[0]), f(mus[1])], axis=1).ravel()
    out["charge"] = np.tile([1, -1], n_events)
    out["mass"] = np.full(2 * n_events, MU_MASS)
    out["nMuon"] = np.full(n_events, 2)
    out["true_mass"] = m
    return out


def _boost(p, u):
    """Boost four-vectors p (E, px, py, pz) from the rest frame of four-velocity u."""
    g = u[:, 0]; b = u[:, 1:] / g[:, None]
    bp = np.sum(b * p[:, 1:], axis=1); b2 = np.sum(b * b, axis=1)
    coef = np.where(b2 > 0, (g - 1) * bp / np.where(b2 > 0, b2, 1), 0.0) + g * p[:, 0]
    out = np.empty_like(p)
    out[:, 0] = g * (p[:, 0] + bp)
    out[:, 1:] = p[:, 1:] + coef[:, None] * b
    return out


__all__ = ["PDG", "GAMMA_Z", "MU_MASS", "pair_mass", "crystal_ball", "voigt", "background", "fit_window",
           "toy_events"]
