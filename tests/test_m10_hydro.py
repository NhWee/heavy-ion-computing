"""Module 10: initial conditions, iSS output reader, and the flow estimator."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
import m10_hydro as m  # noqa: E402


def test_thickness_normalisation():
    x = m.grid(0.1)
    xy = np.array([[0.0, 0.0], [1.0, -2.0], [-3.0, 0.5]])
    g = np.array([0.5, 1.0, 2.0])
    T = m.thickness(xy, g, x)
    assert abs(T.sum() * 0.1 ** 2 - g.sum()) < 1e-6
    assert T.shape == (len(x), len(x))
    # blob at (1, -2): its maximum is at ix for x=1, iy for y=-2
    T1 = m.thickness(xy[1:2], g[1:2], x)
    ix, iy = np.unravel_index(np.argmax(T1), T1.shape)
    assert abs(x[ix] - 1.0) < 1e-9 and abs(x[iy] + 2.0) < 1e-9


def test_eccentricity_of_elliptic_gaussian():
    x = m.grid(0.1)
    X, Y = np.meshgrid(x, x, indexing="ij")
    sx, sy = 2.0, 3.0
    s = np.exp(-X ** 2 / (2 * sx ** 2) - Y ** 2 / (2 * sy ** 2))
    e2, psi2 = m.eccentricity(s, x, 2)
    assert abs(e2 - (sy ** 2 - sx ** 2) / (sy ** 2 + sx ** 2)) < 1e-4   # grid truncation
    # short axis along x -> Psi_2 = 0 (mod pi)
    assert min(abs(psi2), abs(psi2 - np.pi)) < 1e-6


def test_glauber_event_deterministic_and_central():
    a, b = m.glauber_event(7, b=1.0), m.glauber_event(7, b=1.0)
    assert np.allclose(a["xa"], b["xa"]) and np.allclose(a["ga"], b["ga"])
    assert a["npart"] > 300 and len(a["xa"]) + len(a["xb"]) == a["npart"]


def test_read_samples_roundtrip(tmp_path):
    rng = np.random.default_rng(0)
    f = tmp_path / "particle_samples.bin"
    evs = []
    with open(f, "wb") as fh:
        for n in (3, 0, 5):
            a = np.zeros(n, m.REC)
            a["pid"] = rng.choice([211, -211, 2212], n)
            for k in ("m", "t", "x", "y", "z", "E", "px", "py", "pz"):
                a[k] = rng.random(n)
            np.array([n], "<i4").tofile(fh); a.tofile(fh); evs.append(a)
    got = m.read_samples(f)
    assert [len(g) for g in got] == [3, 0, 5]
    assert np.array_equal(got[2]["pid"], evs[2]["pid"]) and np.allclose(got[2]["pz"], evs[2]["pz"])


def _toy_samples(v2, v3, n_samples, mult, rng, psi2=0.3, psi3=1.1):
    """Charged pions with dN/dphi ~ 1 + 2 v2 cos 2(phi-psi2) + 2 v3 cos 3(phi-psi3),
    flat in eta (|eta| < 1) and pT 0.5 GeV."""
    out = []
    fmax = 1 + 2 * v2 + 2 * v3
    for _ in range(n_samples):
        phi = []
        while len(phi) < mult:
            p = rng.uniform(-np.pi, np.pi, 2 * mult)
            f = 1 + 2 * v2 * np.cos(2 * (p - psi2)) + 2 * v3 * np.cos(3 * (p - psi3))
            phi.extend(p[rng.random(len(p)) * fmax < f])
        phi = np.array(phi[:mult])
        eta = rng.uniform(-1, 1, mult); pt = np.full(mult, 0.5); mpi = 0.13957
        a = np.zeros(mult, m.REC)
        a["pid"] = rng.choice([211, -211], mult)
        a["px"], a["py"], a["pz"] = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
        a["E"] = np.sqrt(mpi ** 2 + pt ** 2 * np.cosh(eta) ** 2); a["m"] = mpi
        out.append(a)
    return out


def test_cross_sample_flow_estimator():
    """Unbiased for v_n^2; its statistical error is dominated by the signal x
    noise term, sigma(v_n) ~ 1/sqrt(N_total) with N_total the reference
    particles summed over all samples.  Check at 4 sigma."""
    rng = np.random.default_rng(3)
    o = m.event_observables(_toy_samples(0.08, 0.05, 80, 2000, rng))
    ntot = np.sqrt(o["ref2_den"])                       # ~ N_total
    tol = 4.0 / np.sqrt(ntot)
    v2 = np.sqrt(o["ref2_num"] / o["ref2_den"]); v3 = np.sqrt(max(o["ref3_num"], 0) / o["ref3_den"])
    assert abs(v2 - 0.08) < tol and abs(v3 - 0.05) < tol
    # differential v2 of pions in the pT bin containing 0.5 GeV
    b = np.digitize(0.5, m.PT_EDGES) - 1
    v2pt = o["poi2_pi_num"][b] / o["poi2_pi_den"][b] / v2
    assert abs(v2pt - 0.08) < 2 * tol


def test_music_input_variants():
    txt = m.music_input("ideal", 100.0)
    assert "Viscosity_Flag_Yes_1_No_0 0" in txt and txt.endswith("EndOfData\n")
    txt = m.music_input("bulk", 100.0)
    assert "Include_Bulk_Visc_Yes_1_No_0 1" in txt and "Shear_to_S_ratio 0.12" in txt
    assert "turn_on_bulk=1" in m.iss_args("bulk", 10, 1)
