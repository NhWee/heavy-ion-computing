"""Jet tools for Module 4 (FastJet's Python bindings, classic interface).

SEQUENTIAL RECOMBINATION IN ONE PARAGRAPH
-----------------------------------------
All three algorithms repeat the same step: compute, for every particle i, the
beam distance d_iB = pT_i^(2p) and, for every pair, d_ij = min(pT_i^(2p),
pT_j^(2p)) * dR_ij^2 / R^2; take the smallest; if it is a d_ij merge i and j,
if it is a d_iB call i a jet and remove it.  Only the exponent differs:

    p = +1   kT              soft particles cluster first -> irregular shapes
    p =  0   Cambridge/Aachen   purely geometric          -> irregular shapes
    p = -1   anti-kT         hard particles cluster first -> circles of radius R

All three are infrared and collinear (IRC) safe: adding a soft particle or
splitting one into two collinear ones does not change the hard jets -- the
property that makes jets calculable in perturbative QCD.

JET AREA AND BACKGROUND
-----------------------
The AREA of a jet is measured by adding many infinitely soft "ghost"
particles uniformly in (y, phi) and counting how many end up in the jet.  For
an anti-kT jet it is ~ pi R^2.  A uniform background of transverse-momentum
density rho (GeV per unit area) adds rho * A to every jet; rho is estimated
event by event as the MEDIAN of pT_j / A_j over kT jets (the median ignores
the few hard jets).  This is the standard (Cacciari-Salam) subtraction used in
heavy-ion collisions, where rho reaches hundreds of GeV per unit area.
"""
from __future__ import annotations

import numpy as np

ALGOS = {"antikt": "antikt_algorithm", "kt": "kt_algorithm", "ca": "cambridge_algorithm"}
NEUTRINOS = {12, 14, 16}


def fj():
    import fastjet
    return fastjet


def pseudojets(px, py, pz, e):
    f = fj()
    return [f.PseudoJet(float(a), float(b), float(c), float(d)) for a, b, c, d in zip(px, py, pz, e)]


def jet_definition(algo: str, R: float):
    f = fj()
    return f.JetDefinition(getattr(f, ALGOS[algo]), R)


def area_definition(ghost_maxrap: float = 3.0, ghost_area: float = 0.01, seed: int | None = None):
    """Active area with explicit ghosts.

    Ghost positions are random.  Seeding them would make areas reproducible
    bit for bit, but the PyPI fastjet bindings cannot pass the std::vector<int>
    that set_random_status needs (TypeError), so the seed is applied only
    where the binding accepts it.  Without it, a jet area changes at the
    ~1e-3 level between runs; every physics number quoted is far above that,
    but figure BYTES are not reproducible for Module 4 (docs/troubleshooting.md #20).
    """
    f = fj()
    spec = f.GhostedAreaSpec(ghost_maxrap, 1, ghost_area)
    if seed is not None:
        try:
            spec.set_random_status([int(seed), int(seed) + 1])
        except TypeError:
            pass
    return f.AreaDefinition(f.active_area_explicit_ghosts, spec)


def cluster(particles, algo: str, R: float, pt_min: float, abs_rap_max: float, area_def=None):
    """Cluster and return (jets sorted by pT, cluster sequence).  The sequence
    must be kept alive while the jets are used (constituents, areas)."""
    f = fj()
    jd = jet_definition(algo, R)
    cs = f.ClusterSequenceArea(particles, jd, area_def) if area_def is not None \
        else f.ClusterSequence(particles, jd)
    jets = [j for j in f.sorted_by_pt(cs.inclusive_jets(pt_min)) if abs(j.rap()) < abs_rap_max]
    return jets, cs


def median_rho(particles, area_def, R: float = 0.4, abs_rap_max: float = 2.0, n_hardest_removed: int = 2):
    """Event-by-event background density rho (and its fluctuation sigma) from
    the median of pT/A over kT jets, excluding the n hardest."""
    f = fj()
    sel = f.SelectorAbsRapMax(abs_rap_max) * (~f.SelectorNHardest(n_hardest_removed))
    bge = f.JetMedianBackgroundEstimator(sel, jet_definition("kt", R), area_def)
    bge.set_particles(particles)
    return bge.rho(), bge.sigma()


def toy_background(rng, dndetadphi: float, mean_pt: float, eta_max: float):
    """Uncorrelated 'thermal' background: Poisson multiplicity with density
    dN/(deta dphi), uniform in eta and phi, exponential pT with mean mean_pt,
    massless.  Its true pT density is rho = dndetadphi * mean_pt."""
    n = rng.poisson(dndetadphi * 2 * eta_max * 2 * np.pi)
    pt = rng.exponential(mean_pt, n)
    eta = rng.uniform(-eta_max, eta_max, n)
    phi = rng.uniform(-np.pi, np.pi, n)
    return pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta), pt * np.cosh(eta)


def delta_phi(a, b):
    d = (a - b + np.pi) % (2 * np.pi) - np.pi
    return np.abs(d)


def delta_r(j1, j2):
    return float(np.hypot(j1.rap() - j2.rap(), delta_phi(j1.phi(), j2.phi())))


__all__ = ["ALGOS", "NEUTRINOS", "pseudojets", "jet_definition", "area_definition", "cluster",
           "median_rho", "toy_background", "delta_phi", "delta_r"]
