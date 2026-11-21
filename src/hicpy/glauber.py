"""Monte Carlo Glauber model (Module 6), written from scratch in numpy.

THE MODEL IN ONE PARAGRAPH
--------------------------
A nucleus is a set of A nucleons placed at random according to the measured
nuclear charge density, a Woods-Saxon profile

    rho(r) ~ 1 / (1 + exp((r - R) / a)),

with a hard-core exclusion: no two nucleons closer than d_min.  Two nuclei
are shifted by the impact parameter b in the transverse plane and pass through
each other on straight lines.  Two nucleons collide if their transverse
distance d satisfies  pi d^2 < sigma_NN  (the inelastic nucleon-nucleon cross
section at that energy: the "black disk" criterion).  Then

    N_coll  = number of colliding nucleon pairs   (scales hard processes)
    N_part  = number of nucleons with >= 1 collision ("wounded", scales soft)

and the hadronic nucleus-nucleus cross section is the area pi b_max^2 times
the fraction of events with at least one collision.  Centrality classes are
percentiles of the cross section, ordered here by impact parameter (the
geometric definition); experiments order by a measured multiplicity, which
agrees closely for all but the most peripheral classes.

DEFAULT INPUTS (Pb-208 at sqrt(s_NN) = 5.02 TeV)
-------------------------------------------------
R = 6.62 +- 0.06 fm, a = 0.546 +- 0.010 fm, d_min = 0.4 +- 0.4 fm
(ALICE, arXiv:1301.4361, from low-energy electron scattering), and
sigma_NN = 67.6 +- 0.6 mb (arXiv:1710.07098).  1 mb = 0.1 fm^2.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MB_TO_FM2 = 0.1


@dataclass
class Nucleus:
    A: int
    R: float = 6.62        # fm
    a: float = 0.546       # fm
    d_min: float = 0.4     # fm


PB208 = Nucleus(208, 6.62, 0.546, 0.4)
PROTON = Nucleus(1, 0.0, 0.0, 0.0)


def _radial_cdf(R, a, r_max=None, n=4000):
    r_max = r_max or R + 10 * a
    r = np.linspace(0, r_max, n)
    pdf = r ** 2 / (1 + np.exp((r - R) / a))
    cdf = np.cumsum(pdf); cdf /= cdf[-1]
    return r, cdf


def sample_nucleus(rng, nuc: Nucleus, recentre: bool = True):
    """(A, 3) nucleon positions in fm.  r^2 * Woods-Saxon radial density,
    isotropic angles, and every nucleon at least d_min from all earlier ones
    (a rejected nucleon is redrawn).  Recentring puts the centre of mass at
    the origin, as TGlauberMC does."""
    if nuc.A == 1:
        return np.zeros((1, 3))
    r_grid, cdf = _radial_cdf(nuc.R, nuc.a)
    pos = np.empty((nuc.A, 3))
    k = 0
    while k < nuc.A:
        r = np.interp(rng.random(), cdf, r_grid)
        ct = 2 * rng.random() - 1; st = np.sqrt(1 - ct * ct); ph = 2 * np.pi * rng.random()
        p = r * np.array([st * np.cos(ph), st * np.sin(ph), ct])
        if k and nuc.d_min > 0 and np.min(np.sum((pos[:k] - p) ** 2, axis=1)) < nuc.d_min ** 2:
            continue
        pos[k] = p; k += 1
    if recentre:
        pos -= pos.mean(axis=0)
    return pos


def collide(pa, pb, b, sigma_nn_mb):
    """Nucleus A shifted by +b/2, B by -b/2 along x.  Returns (n_part, n_coll,
    wounded_a mask, wounded_b mask)."""
    xa = pa[:, :2] + np.array([b / 2, 0.0])
    xb = pb[:, :2] - np.array([b / 2, 0.0])
    d2 = np.sum((xa[:, None, :] - xb[None, :, :]) ** 2, axis=2)
    hit = d2 < sigma_nn_mb * MB_TO_FM2 / np.pi
    wa, wb = hit.any(axis=1), hit.any(axis=0)
    return int(wa.sum() + wb.sum()), int(hit.sum()), wa, wb


def run(n_events, nuc_a: Nucleus, nuc_b: Nucleus, sigma_nn_mb, b_max=20.0, seed=1):
    """Minimum-bias Glauber events, b drawn with dP ~ b db up to b_max.
    Returns dict of arrays over events WITH at least one collision, plus the
    hadronic cross section in barn."""
    rng = np.random.default_rng(seed)
    b = b_max * np.sqrt(rng.random(n_events))
    npart = np.zeros(n_events, int); ncoll = np.zeros(n_events, int)
    for i in range(n_events):
        npart[i], ncoll[i], _, _ = collide(sample_nucleus(rng, nuc_a), sample_nucleus(rng, nuc_b),
                                           b[i], sigma_nn_mb)
    hit = ncoll > 0
    frac = hit.mean()
    sigma_b = np.pi * b_max ** 2 * frac / 100.0          # fm^2 -> barn
    sigma_err = np.pi * b_max ** 2 * np.sqrt(frac * (1 - frac) / n_events) / 100.0
    return {"b": b[hit], "npart": npart[hit], "ncoll": ncoll[hit],
            "sigma_inel_b": float(sigma_b), "sigma_inel_err_b": float(sigma_err),
            "n_events_thrown": n_events}


def centrality_table(res, edges_pct=(0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), order="b"):
    """<N_part>, <N_coll>, <b> in centrality classes (percent of the cross
    section).  order='b': geometric (smallest b = most central);
    order='npart': by N_part (largest = most central), the closer analogue of
    a multiplicity-based experimental definition."""
    key = res["b"] if order == "b" else -res["npart"]
    idx = np.argsort(key, kind="stable")
    n = len(idx); rows = []
    for lo, hi in zip(edges_pct[:-1], edges_pct[1:]):
        sel = idx[int(n * lo / 100):int(n * hi / 100)]
        rows.append({"class": f"{lo}-{hi}%", "npart": float(res["npart"][sel].mean()),
                     "npart_rms": float(res["npart"][sel].std()), "ncoll": float(res["ncoll"][sel].mean()),
                     "b": float(res["b"][sel].mean()), "n": int(len(sel))})
    return rows


# ----------------------------------------------------------------------------- geometry for Module 9
def event_geometry(pa, pb, b, sigma_nn_mb):
    """Transverse geometry of one collision: midpoints of the colliding nucleon
    pairs (where hard processes happen, density ~ N_coll) and the positions of
    the wounded nucleons (the source of the medium, density ~ N_part).
    Returns (collision_points (N_coll, 2), wounded (N_part, 2)) in fm."""
    xa = pa[:, :2] + np.array([b / 2, 0.0])
    xb = pb[:, :2] - np.array([b / 2, 0.0])
    d2 = np.sum((xa[:, None, :] - xb[None, :, :]) ** 2, axis=2)
    hit = d2 < sigma_nn_mb * MB_TO_FM2 / np.pi
    ia, ib = np.nonzero(hit)
    coll = 0.5 * (xa[ia] + xb[ib])
    wounded = np.vstack([xa[hit.any(axis=1)], xb[hit.any(axis=0)]])
    return coll, wounded


def participant_plane(wounded, n=2):
    """Participant-plane angle Psi_n and eccentricity eps_n of the wounded
    nucleons (r^n weighting, centred).  Psi_2 points along the SHORT axis of
    the almond (along b for a smooth overlap): the in-plane direction, where
    the medium is thinnest, so that both flow and energy loss give v_2 > 0
    with respect to it."""
    x = wounded - wounded.mean(axis=0)
    r2 = np.sum(x ** 2, axis=1); phi = np.arctan2(x[:, 1], x[:, 0])
    c, s_ = np.sum(r2 ** (n / 2) * np.cos(n * phi)), np.sum(r2 ** (n / 2) * np.sin(n * phi))
    eps = np.hypot(c, s_) / np.sum(r2 ** (n / 2))
    psi = (np.arctan2(s_, c) + np.pi) / n          # direction of the shortest extent
    return float(psi), float(eps)


def path_integrals(x0, phi, sources, width=0.5, zeta=0):
    """Line integrals of the medium density along straight paths.

    The density is a sum of 2D Gaussians (width fm) on the wounded nucleons,
    rho(x) = sum_j G(x - s_j), in nucleons / fm^2.  For a path starting at x0
    in direction phi, returns I_zeta = integral_0^inf rho(x0 + t n) t^zeta dt:
    zeta = 0 -> number of participants crossed per fm^2 times fm (~ L),
    zeta = 1 -> weighted by the distance travelled (~ L^2, the BDMPS-like
    radiative case).  Closed form per Gaussian, vectorised (paths x sources).
    """
    from scipy.special import erf
    x0 = np.atleast_2d(x0); phi = np.atleast_1d(phi)
    nvec = np.stack([np.cos(phi), np.sin(phi)], axis=1)                  # (k, 2)
    d = sources[None, :, :] - x0[:, None, :]                            # (k, m, 2)
    tc = np.einsum("kmi,ki->km", d, nvec)                               # projection on the path
    dperp2 = np.sum(d ** 2, axis=2) - tc ** 2
    w2 = width ** 2
    pref = np.exp(-dperp2 / (2 * w2)) / (2 * np.pi * w2)
    half = 0.5 * (1 + erf(tc / (width * np.sqrt(2))))                   # fraction of the Gaussian ahead
    i0 = width * np.sqrt(2 * np.pi) * half                              # int_0^inf exp(-(t-tc)^2/2w^2) dt
    if zeta == 0:
        return np.sum(pref * i0, axis=1)
    i1 = w2 * np.exp(-tc ** 2 / (2 * w2)) + tc * i0                    # int_0^inf t exp(...) dt
    return np.sum(pref * i1, axis=1)


__all__ = ["Nucleus", "PB208", "PROTON", "sample_nucleus", "collide", "run", "centrality_table",
           "event_geometry", "participant_plane", "path_integrals"]
