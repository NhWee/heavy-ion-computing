"""A physically-motivated toy event generator (Module 1).

WHY A TOY GENERATOR AT ALL?
---------------------------
Module 1 is about the *data format and the analysis chain*, not about QCD.  We
need a particle-level event sample before PYTHIA 8 is installed (Module 2), and
it must be realistic enough that the histograms we make look like real ones --
otherwise we cannot tell a working analysis from a broken one.  Every
distribution below is therefore taken from a functional form that experiments
actually fit to data, with parameters in the right ballpark for pp collisions at
LHC energies.  It is NOT a QCD calculation: there is no hard scattering, no
parton shower, no hadronisation, no correlations between particles.  Module 2
replaces it with PYTHIA, and the *same* analysis code must keep working.

THE FOUR INGREDIENTS
--------------------
1. Multiplicity  N_ch  ~ Negative Binomial Distribution (NBD).
   P(n) = Gamma(n+k)/(Gamma(k) n!) * (nbar/k)^n / (1+nbar/k)^(n+k)
   The NBD has described charged-multiplicity distributions in hadronic
   collisions since the UA5 measurements; it is wider than a Poisson
   (var = nbar + nbar^2/k), which is what "event-by-event fluctuations" means
   in practice.  k -> infinity recovers Poisson.

2. Pseudorapidity  dN/deta: flat plateau for |eta| < eta0 with Gaussian
   shoulders.  The plateau is the boost-invariant central region of a
   high-energy collision; the shoulders are the fragmentation regions.

3. Transverse momentum: Tsallis / Hagedorn form, the standard parametrisation
   ALICE and CMS fit to identified-particle spectra:

       1/(2 pi p_T) d^2N/(dp_T dy) = C * [ 1 + (m_T - m0)/(n T) ]^(-n)

   with m_T = sqrt(p_T^2 + m0^2).  Small p_T is thermal-like (exponential,
   temperature T); large p_T becomes a power law p_T^(-n) -- the hard-scattering
   tail.  A single exponential would badly undershoot the high-p_T tail, which
   is exactly the physics point.

4. Azimuth phi: uniform.  (No flow, no jets -- azimuthal structure is the
   subject of Modules 8 and 10.)

Optionally a Z -> mu+ mu- "signal" is embedded in a fraction of events, so that
Module 1 also exercises invariant-mass reconstruction: the dimuon mass is drawn
from a relativistic Breit-Wigner at m_Z = 91.1876 GeV, Gamma_Z = 2.4952 GeV
(PDG 2024 values), and the decay is isotropic in the Z rest frame.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict

import numpy as np

# --- PDG particle data actually used here (PDG 2024 central values, GeV) -----
M_PI = 0.13957039     # pi+-
M_K = 0.493677        # K+-
M_P = 0.93827208816   # proton
M_MU = 0.1056583755   # muon
M_Z = 91.1876         # Z boson
GAMMA_Z = 2.4952      # Z total width

PDG = {"pi": 211, "K": 321, "p": 2212, "mu": 13}
MASS = {211: M_PI, 321: M_K, 2212: M_P, 13: M_MU}


@dataclass
class ToyConfig:
    """Every knob of the generator. Written to disk with the events."""
    # --- event sample ---
    n_events: int = 2000
    seed: int = 20260827          # fixed and logged: reproducibility rule
    sqrt_s_TeV: float = 13.0      # label only; the shape parameters encode it

    # --- multiplicity (NBD) ---
    nbd_nbar: float = 60.0        # mean charged multiplicity in |eta| < 3
    nbd_k: float = 2.2            # NBD shape; smaller k = wider fluctuations

    # --- pseudorapidity ---
    eta_plateau: float = 2.0      # |eta| < eta_plateau is flat
    eta_sigma: float = 1.6        # Gaussian shoulder width beyond the plateau
    eta_max: float = 6.0          # generation range

    # --- transverse momentum (Tsallis) ---
    tsallis_T: float = 0.135      # GeV, "temperature"-like slope at low pT
    tsallis_n: float = 6.8        # power-law index of the hard tail
    pt_min: float = 0.15          # GeV/c, generation threshold
    pt_max: float = 60.0          # GeV/c

    # --- charged-hadron composition (approximate LHC pp values) ---
    frac_pi: float = 0.83
    frac_K: float = 0.12
    frac_p: float = 0.05

    # --- embedded Z -> mu mu signal ---
    z_fraction: float = 0.02      # fraction of events containing a Z
    z_pt_scale: float = 8.0       # GeV, exponential scale of the Z pT
    z_y_sigma: float = 1.8        # width of the Z rapidity distribution

    def to_dict(self):
        return asdict(self)


# --------------------------------------------------------------------------
# sampling helpers
# --------------------------------------------------------------------------
def _eta_shape(eta, plateau, sigma):
    """Unnormalised dN/deta: flat core + Gaussian shoulders."""
    x = np.abs(eta)
    out = np.ones_like(x)
    tail = x > plateau
    out[tail] = np.exp(-0.5 * ((x[tail] - plateau) / sigma) ** 2)
    return out


def _sample_eta(rng, n, cfg):
    """Inverse-CDF sampling on a fine grid (fast, exact enough, no rejection loop)."""
    grid = np.linspace(-cfg.eta_max, cfg.eta_max, 4001)
    pdf = _eta_shape(grid, cfg.eta_plateau, cfg.eta_sigma)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    return np.interp(rng.random(n), cdf, grid)


def _tsallis_dNdpt(ptv, mass, T, n):
    """dN/dp_T (note the extra p_T from the Jacobian of the invariant form)."""
    mt = np.sqrt(ptv * ptv + mass * mass)
    return ptv * (1.0 + (mt - mass) / (n * T)) ** (-n)


def _sample_pt(rng, n, mass, cfg):
    """Inverse-CDF sampling of the Tsallis spectrum on a log grid."""
    grid = np.geomspace(cfg.pt_min, cfg.pt_max, 3000)
    pdf = _tsallis_dNdpt(grid, mass, cfg.tsallis_T, cfg.tsallis_n)
    # trapezoidal CDF on a non-uniform grid
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(grid))])
    cdf /= cdf[-1]
    return np.interp(rng.random(n), cdf, grid)


def _sample_breit_wigner(rng, n, m0, gamma):
    """Relativistic Breit-Wigner, sampled by rejection against a Cauchy proposal."""
    out = np.empty(n)
    filled = 0
    while filled < n:
        trial = m0 + 0.5 * gamma * np.tan(np.pi * (rng.random(n) - 0.5))
        trial = trial[(trial > m0 - 20 * gamma) & (trial < m0 + 20 * gamma)]
        # weight = relativistic BW / Cauchy proposal
        num = trial ** 2 * gamma ** 2
        den = (trial ** 2 - m0 ** 2) ** 2 + m0 ** 2 * gamma ** 2
        w = (num / den) / (1.0 / (1.0 + ((trial - m0) / (0.5 * gamma)) ** 2))
        w /= w.max()
        acc = trial[rng.random(trial.size) < w]
        take = min(acc.size, n - filled)
        out[filled:filled + take] = acc[:take]
        filled += take
    return out


def _boost_z_decay(rng, m_zz, pt_z, phi_z, y_z):
    """Two-body isotropic decay Z -> mu+ mu-, then boost to the lab frame.

    Returns two (px, py, pz, E) tuples.
    """
    n = m_zz.size
    # --- rest frame: back-to-back muons, isotropic ---
    p_star = np.sqrt(np.clip(0.25 * m_zz ** 2 - M_MU ** 2, 0.0, None))
    cth = 2.0 * rng.random(n) - 1.0
    sth = np.sqrt(1.0 - cth ** 2)
    ph = 2.0 * np.pi * rng.random(n)
    px1 = p_star * sth * np.cos(ph)
    py1 = p_star * sth * np.sin(ph)
    pz1 = p_star * cth
    e1 = np.sqrt(p_star ** 2 + M_MU ** 2)

    # --- lab-frame Z four-momentum from (pT, phi, y, m) ---
    mt_z = np.sqrt(m_zz ** 2 + pt_z ** 2)
    pzz = mt_z * np.sinh(y_z)
    ez = mt_z * np.cosh(y_z)
    bx, by, bz = pt_z * np.cos(phi_z) / ez, pt_z * np.sin(phi_z) / ez, pzz / ez
    b2 = bx * bx + by * by + bz * bz
    gam = 1.0 / np.sqrt(np.clip(1.0 - b2, 1e-16, None))

    def boost(px, py, pz, e):
        bp = bx * px + by * py + bz * pz
        f = np.where(b2 > 0, (gam - 1.0) * bp / np.where(b2 > 0, b2, 1.0), 0.0)
        return (px + f * bx + gam * bx * e,
                py + f * by + gam * by * e,
                pz + f * bz + gam * bz * e,
                gam * (e + bp))

    mu_plus = boost(px1, py1, pz1, e1)
    mu_minus = boost(-px1, -py1, -pz1, e1)
    return mu_plus, mu_minus


# --------------------------------------------------------------------------
# main entry point
# --------------------------------------------------------------------------
def generate(cfg: ToyConfig):
    """Generate the event sample.

    Returns a dict of per-event lists (ragged): px, py, pz, E, pdg
    plus per-event scalars: n_particles, has_z.
    """
    rng = np.random.default_rng(cfg.seed)

    # NBD via numpy's negative_binomial(n=k, p=k/(k+nbar))
    p_nb = cfg.nbd_k / (cfg.nbd_k + cfg.nbd_nbar)
    mult = rng.negative_binomial(cfg.nbd_k, p_nb, size=cfg.n_events).astype(int)
    mult = np.maximum(mult, 1)          # keep at least one particle per event

    has_z = rng.random(cfg.n_events) < cfg.z_fraction
    n_z = int(has_z.sum())

    # Z kinematics for the events that have one
    if n_z:
        m_zz = _sample_breit_wigner(rng, n_z, M_Z, GAMMA_Z)
        pt_z = rng.exponential(cfg.z_pt_scale, n_z)
        phi_z = rng.uniform(-np.pi, np.pi, n_z)
        y_z = rng.normal(0.0, cfg.z_y_sigma, n_z)
        mup, mum = _boost_z_decay(rng, m_zz, pt_z, phi_z, y_z)

    species = np.array([PDG["pi"], PDG["K"], PDG["p"]])
    weights = np.array([cfg.frac_pi, cfg.frac_K, cfg.frac_p], dtype=float)
    weights /= weights.sum()

    px_l, py_l, pz_l, e_l, id_l = [], [], [], [], []
    iz = 0
    for iev in range(cfg.n_events):
        n = mult[iev]
        kind = rng.choice(species, size=n, p=weights)
        sign = rng.choice([-1, 1], size=n)
        pdg = kind * sign
        masses = np.array([MASS[k] for k in kind])

        ptv = np.empty(n)
        for s in species:                      # sample per species (mass matters)
            m = kind == s
            if m.any():
                ptv[m] = _sample_pt(rng, int(m.sum()), MASS[s], cfg)
        etav = _sample_eta(rng, n, cfg)
        phiv = rng.uniform(-np.pi, np.pi, n)

        px = ptv * np.cos(phiv)
        py = ptv * np.sin(phiv)
        pz = ptv * np.sinh(etav)
        e = np.sqrt(px * px + py * py + pz * pz + masses * masses)

        if has_z[iev]:
            px = np.concatenate([px, [mup[0][iz], mum[0][iz]]])
            py = np.concatenate([py, [mup[1][iz], mum[1][iz]]])
            pz = np.concatenate([pz, [mup[2][iz], mum[2][iz]]])
            e = np.concatenate([e, [mup[3][iz], mum[3][iz]]])
            pdg = np.concatenate([pdg, [-13, 13]])   # mu+ has pdg -13
            iz += 1

        px_l.append(px); py_l.append(py); pz_l.append(pz); e_l.append(e)
        id_l.append(pdg.astype(np.int32))

    return {
        "px": px_l, "py": py_l, "pz": pz_l, "E": e_l, "pdg": id_l,
        "n_particles": np.array([len(a) for a in px_l], dtype=np.int32),
        "has_z": has_z.astype(bool),
    }
