"""Collider kinematics helpers.

Conventions used throughout this project
----------------------------------------
* Natural units: hbar = c = 1. Momenta and masses in GeV; when we write
  "p_T (GeV/c)" on an axis we are using the experimentalists' convention of
  keeping the c explicit in the label only.
* Beam axis = z.  The transverse plane is (x, y).
* Pseudorapidity   eta = -ln[ tan(theta/2) ]      (purely geometric)
* Rapidity          y  = 0.5 * ln[(E+pz)/(E-pz)]  (additive under z-boosts)
  eta -> y in the massless limit p >> m.  Experiments that only track a
  particle's direction (no mass hypothesis) report eta; theorists prefer y
  because dN/dy is boost-invariant in shape.
"""
from __future__ import annotations

import numpy as np


def pt(px, py):
    """Transverse momentum sqrt(px^2 + py^2). Invariant under boosts along z."""
    return np.hypot(px, py)


def p_mag(px, py, pz):
    return np.sqrt(px * px + py * py + pz * pz)


def eta(px, py, pz):
    """Pseudorapidity. Uses the numerically stable arctanh(pz/|p|) form."""
    p = p_mag(px, py, pz)
    # clip guards against |pz| == |p| (exactly-collinear particles -> +-inf)
    return np.arctanh(np.clip(pz / np.where(p > 0, p, 1e-300), -1 + 1e-15, 1 - 1e-15))


def rapidity(e, pz):
    """Rapidity y = 0.5 ln[(E+pz)/(E-pz)]."""
    num = np.clip(e + pz, 1e-300, None)
    den = np.clip(e - pz, 1e-300, None)
    return 0.5 * np.log(num / den)


def phi(px, py):
    """Azimuthal angle in (-pi, pi]."""
    return np.arctan2(py, px)


def energy(px, py, pz, mass):
    """On-shell energy E = sqrt(p^2 + m^2)."""
    return np.sqrt(px * px + py * py + pz * pz + mass * mass)


def mt(pt_val, mass):
    """Transverse mass m_T = sqrt(p_T^2 + m^2) (not the 'transverse mass' of W physics)."""
    return np.sqrt(pt_val * pt_val + mass * mass)


def from_ptetaphim(pt_val, eta_val, phi_val, mass):
    """(pT, eta, phi, m) -> (px, py, pz, E). The collider-standard parametrisation."""
    px = pt_val * np.cos(phi_val)
    py = pt_val * np.sin(phi_val)
    pz = pt_val * np.sinh(eta_val)
    e = energy(px, py, pz, mass)
    return px, py, pz, e


def invariant_mass(e, px, py, pz):
    """m = sqrt(E^2 - |p|^2) for a (possibly summed) four-vector.

    Negative values under the sqrt can appear from floating-point round-off for
    (nearly) massless systems; they are clipped to zero rather than silently
    producing NaN.
    """
    m2 = e * e - (px * px + py * py + pz * pz)
    return np.sqrt(np.clip(m2, 0.0, None))


def invariant_mass_pair(p1, p2):
    """Invariant mass of two particles given as (px, py, pz, E) tuples/arrays."""
    px = p1[0] + p2[0]
    py = p1[1] + p2[1]
    pz = p1[2] + p2[2]
    e = p1[3] + p2[3]
    return invariant_mass(e, px, py, pz)


def delta_phi(phi1, phi2):
    """Signed azimuthal separation folded into (-pi, pi]."""
    d = phi1 - phi2
    return (d + np.pi) % (2 * np.pi) - np.pi


__all__ = ["pt", "p_mag", "eta", "rapidity", "phi", "energy", "mt",
           "from_ptetaphim", "invariant_mass", "invariant_mass_pair", "delta_phi"]
