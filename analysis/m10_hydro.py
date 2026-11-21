#!/usr/bin/env python
"""Module 10: event-by-event viscous hydrodynamics for Pb-Pb at 5.02 TeV.

THE CHAIN
---------
1. Initial condition (this file): Module 6 Monte Carlo Glauber nuclei.  Each
   wounded nucleon deposits a Gaussian blob (width w = 0.5 fm) with a
   fluctuating weight gamma ~ Gamma(k, 1/k) (mean 1).  The two nuclei give
   thickness functions T_A, T_B and the entropy density at tau_0 is taken as
   their geometric mean,
        s(x, y, tau_0) = s_factor * sqrt(T_A T_B),
   the "reduced thickness" with p = 0 of the TRENTo parametrisation
   (Moreland, Bernhard, Bass, PRC 92 011901 (2015), arXiv:1412.4708).
   s_factor is the one number tuned to data: dN_ch/deta (|eta| < 0.5) in
   0-5 % Pb-Pb = 1943 (ALICE, arXiv:1512.06104).
2. MUSIC (2+1D, boost invariant): second-order viscous hydrodynamics with
   lattice-QCD equation of state (hotQCD + hadron resonance gas, EOS_to_use 9),
   constant eta/s, optional bulk viscosity zeta/s(T).  Starts at tau_0 = 0.6
   fm/c from rest (no pre-equilibrium flow).  The freeze-out surface at
   energy density e_sw = 0.18 GeV/fm^3 is found with Cornelius.
3. iSS: Cooper-Frye sampling of hadrons on that surface (with the viscous
   delta f), then resonance decays.  No hadronic afterburner (UrQMD/SMASH):
   hadrons stop interacting at e_sw.  Each hydro event is sampled M times
   ("oversampling") so that a smooth hydro event gives many hadron events.

OBSERVABLES
-----------
dN_ch/deta, identified pi/K/p spectra and <pT> (radial flow), and flow
coefficients v_n{2}.  Because the M samples of one hydro event are independent
given the fluid, correlating particles from DIFFERENT samples measures the
flow of that event with no non-flow (decays, self-correlations):
    v_n^2(event) = (|sum_k Q_k|^2 - sum_k |Q_k|^2) / ((sum N_k)^2 - sum N_k^2),
Q_k = sum_j exp(i n phi_j) over the particles of sample k.  The class value
v_n{2} = sqrt(<v_n^2>) with events weighted by their number of pairs, like
the experimental two-particle cumulant.

USAGE
-----
  bash scripts/m10_build_hydro.sh                            # MUSIC + iSS into external/
  python analysis/m10_hydro.py scan --n 4000                  # Glauber min-bias -> centrality classes
  python analysis/m10_hydro.py calibrate --variant nominal    # s_factor from 0-5 % dN/deta
  python analysis/m10_hydro.py run --variant nominal --per-class 20 --jobs 4
  python analysis/m10_hydro.py analyze
`run` appends to results/tables/m10/m10_events_<variant>.jsonl and skips events
already there, so it can be interrupted and resumed.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from hicpy.glauber import PB208, collide, sample_nucleus  # noqa: E402

EXT = Path(os.environ.get("HIC_EXTERNAL", ROOT / "external"))
MUSIC_DIR, ISS_DIR = EXT / "MUSIC", EXT / "iSS"
MUSIC_EXE, ISS_EXE = MUSIC_DIR / "build/src/MUSIChydro", ISS_DIR / "build/src/iSS.e"
TAB, FIG = ROOT / "results/tables/m10", ROOT / "results/figures/m10"
WORK = ROOT / "data/generated/m10"

# ---- initial condition ----------------------------------------------------
SIGMA_NN = 67.6          # mb, Pb-Pb 5.02 TeV (Module 6)
B_MAX = 20.0             # fm
W_NUC = 0.5              # fm, nucleon (blob) width
K_GAMMA = 1.6            # shape of the Gamma-distributed blob weights (a choice; TRENTo-type fits give O(1))
X_MAX, DX = 15.0, 0.2    # hydro grid: [-15, 15] fm, 151 x 151 cells
TAU0 = 0.6               # fm/c
SEED0 = 20261121
CLASSES = [(0, 5), (5, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60)]
DNDETA_0_5 = (1943.0, 54.0, "ALICE, arXiv:1512.06104, |eta| < 0.5")

# ---- hydro and particlization ---------------------------------------------
E_SW = 0.18              # GeV/fm^3
VARIANTS = {
    "nominal": {"eta_s": 0.12, "bulk": False, "label": r"$\eta/s = 0.12$"},
    "ideal":   {"eta_s": 0.0,  "bulk": False, "label": r"ideal ($\eta/s = 0$)"},
    "etas020": {"eta_s": 0.20, "bulk": False, "label": r"$\eta/s = 0.20$"},
    "bulk":    {"eta_s": 0.12, "bulk": True,  "label": r"$\eta/s = 0.12$ + $\zeta/s(T)$"},
}
Y_SAMPLE = 1.5           # hadrons sampled in |y| < 1.5 (boost invariant)
ETA_FLOW, PT_FLOW = 0.8, (0.2, 3.0)    # ALICE-like reference particles
PT_EDGES = np.round(np.arange(0.0, 3.01, 0.1), 3)
# charges of the hadrons and leptons that survive iSS's decays (UrQMD table)
CHARGE = {211: 1, 321: 1, 2212: 1, 3222: 1, 3112: -1, 3312: -1, 3334: -1, 11: -1, 13: -1,
          111: 0, 130: 0, 310: 0, 311: 0, 2112: 0, 3122: 0, 3212: 0, 3322: 0, 22: 0}
SPECIES = {"pi": 211, "K": 321, "p": 2212}


def charge(pid):
    a = np.abs(pid)
    q = np.array([CHARGE.get(int(x), 99) for x in a])
    return np.where(q == 99, 99, np.sign(pid) * q)


# ===========================================================================
# initial conditions
# ===========================================================================
def glauber_event(i, b=None):
    """Nucleon configuration of event i (deterministic in i)."""
    rng = np.random.default_rng([SEED0, i])
    if b is None:
        b = B_MAX * np.sqrt(rng.random())
    pa, pb = sample_nucleus(rng, PB208), sample_nucleus(rng, PB208)
    npart, ncoll, wa, wb = collide(pa, pb, b, SIGMA_NN)
    xa = pa[wa, :2] + np.array([b / 2, 0.0])
    xb = pb[wb, :2] - np.array([b / 2, 0.0])
    ga = rng.gamma(K_GAMMA, 1.0 / K_GAMMA, len(xa))
    gb = rng.gamma(K_GAMMA, 1.0 / K_GAMMA, len(xb))
    return {"i": i, "b": float(b), "npart": npart, "ncoll": ncoll, "xa": xa, "xb": xb, "ga": ga, "gb": gb}


def thickness(xy, gam, x):
    """sum_i gam_i G(x - x_i) G(y - y_i) on the grid x (same in y): separable,
    so one matrix product.  Units fm^-2 (each blob integrates to gam_i)."""
    if len(xy) == 0:
        return np.zeros((len(x), len(x)))
    gx = np.exp(-(x[:, None] - xy[None, :, 0]) ** 2 / (2 * W_NUC ** 2))
    gy = np.exp(-(x[:, None] - xy[None, :, 1]) ** 2 / (2 * W_NUC ** 2))
    return (gx * gam) @ gy.T / (2 * np.pi * W_NUC ** 2)        # [ix, iy]


def grid(dx=DX):
    n = int(round(2 * X_MAX / dx)) + 1
    return np.linspace(-X_MAX, X_MAX, n)


def entropy_profile(ev, x):
    return np.sqrt(thickness(ev["xa"], ev["ga"], x) * thickness(ev["xb"], ev["gb"], x))


def eccentricity(s, x, n):
    """eps_n and participant-plane angle of a profile s[ix, iy] (r^n weight,
    about its centroid; Psi_n points along the short axis for n = 2)."""
    X, Y = np.meshgrid(x, x, indexing="ij")
    w = s / s.sum()
    xc, yc = (w * X).sum(), (w * Y).sum()
    r = np.hypot(X - xc, Y - yc); ph = np.arctan2(Y - yc, X - xc)
    c, sn = (w * r ** n * np.cos(n * ph)).sum(), (w * r ** n * np.sin(n * ph)).sum()
    return float(np.hypot(c, sn) / (w * r ** n).sum()), float((np.arctan2(sn, c) + np.pi) / n)


def _scan_one(i):
    ev = glauber_event(i)
    if ev["npart"] == 0:
        return None
    x = grid(0.5)
    s = entropy_profile(ev, x)
    return {"i": i, "b": ev["b"], "npart": ev["npart"], "ncoll": ev["ncoll"], "S": float(s.sum() * 0.25)}


def cmd_scan(a):
    """Min-bias Glauber events -> centrality classes by total initial entropy
    (the model's multiplicity: experiments classify by measured multiplicity)."""
    with Pool(a.jobs) as p:
        rows = [r for r in p.map(_scan_one, range(a.n), chunksize=50) if r]
    S = np.array([r["S"] for r in rows])
    out = {"n_generated": a.n, "n_inelastic": len(rows),
           "sigma_inel_b": float(np.pi * B_MAX ** 2 * len(rows) / a.n / 100.0),
           "settings": {"w_fm": W_NUC, "k_gamma": K_GAMMA, "sigma_nn_mb": SIGMA_NN, "seed0": SEED0},
           "classes": []}
    for lo, hi in CLASSES:
        s_hi = np.percentile(S, 100 - lo) if lo > 0 else np.inf
        s_lo = np.percentile(S, 100 - hi)
        sel = [r for r in rows if s_lo <= r["S"] < s_hi]
        out["classes"].append({"class": f"{lo}-{hi}%", "S_range": [float(s_lo), float(min(s_hi, 1e30))],
                               "n": len(sel), "b": float(np.mean([r["b"] for r in sel])),
                               "npart": float(np.mean([r["npart"] for r in sel])),
                               "ncoll": float(np.mean([r["ncoll"] for r in sel])),
                               "S": float(np.mean([r["S"] for r in sel])),
                               "events": [r["i"] for r in sel]})
    TAB.mkdir(parents=True, exist_ok=True)
    (TAB / "m10_centrality.json").write_text(json.dumps(out, indent=1))
    print(f"sigma_inel = {out['sigma_inel_b']:.2f} b from {len(rows)}/{a.n} events")
    for c in out["classes"]:
        print(f"  {c['class']:>7}: n={c['n']:4d}  <b>={c['b']:.2f}  <Npart>={c['npart']:.1f}  <S>={c['S']:.1f}")


def load_centrality():
    f = TAB / "m10_centrality.json"
    if not f.exists():
        sys.exit("run `python analysis/m10_hydro.py scan` first")
    return json.loads(f.read_text())


# ===========================================================================
# MUSIC + iSS for one event
# ===========================================================================
def music_input(var, s_factor):
    v = VARIANTS[var]
    visc = 1 if (v["eta_s"] > 0 or v["bulk"]) else 0
    p = {
        "echo_level": 1, "mode": 2, "Initial_profile": 8, "initialize_with_entropy": 1,
        "Initial_Distribution_input_filename": "initial.dat", "s_factor": f"{s_factor:.6g}",
        "Eta_plateau_size": 20.0, "Eta_fall_off": 0.7, "boost_invariant": 1,
        "Initial_time_tau_0": TAU0, "Total_evolution_time_tau": 40.0, "Delta_Tau": 0.02,
        "Eta_grid_size": 14.0, "Grid_size_in_eta": 1,
        "X_grid_size_in_fm": 2 * X_MAX, "Y_grid_size_in_fm": 2 * X_MAX,
        "Grid_size_in_x": len(grid()), "Grid_size_in_y": len(grid()),
        "EOS_to_use": 9,
        "Viscosity_Flag_Yes_1_No_0": visc, "Include_Shear_Visc_Yes_1_No_0": 1 if v["eta_s"] > 0 else 0,
        "Shear_to_S_ratio": v["eta_s"], "T_dependent_Shear_to_S_ratio": 0,
        "Include_Bulk_Visc_Yes_1_No_0": 1 if v["bulk"] else 0, "T_dependent_Bulk_to_S_ratio": 8,
        "Include_second_order_terms": 1, "Include_Rhob_Yes_1_No_0": 0, "turn_on_baryon_diffusion": 0,
        "output_evolution_data": 0, "output_movie_flag": 0,
        "Do_FreezeOut_Yes_1_No_0": 1, "Do_FreezeOut_lowtemp": 1, "freeze_out_method": 4,
        "freeze_surface_in_binary": 1, "average_surface_over_this_many_time_steps": 5,
        "freeze_Ncell_x_step": 1, "freeze_Ncell_eta_step": 1, "freeze_eps_flag": 0, "N_freeze_out": 1,
        "eps_freeze_max": E_SW, "eps_freeze_min": E_SW, "use_eps_for_freeze_out": 1,
    }
    return "".join(f"{k} {val}\n" for k, val in p.items()) + "EndOfData\n"


def iss_args(var, n_samples, seed):
    bulk = 1 if VARIANTS[var]["bulk"] else 0
    p = {"hydro_mode": 1, "afterburner_type": 1, "turn_on_bulk": bulk, "include_deltaf_bulk": bulk,
         "include_deltaf_shear": 1 if VARIANTS[var]["eta_s"] > 0 else 0,
         "perform_decays": 1, "MC_sampling": 4, "number_of_repeated_sampling": n_samples,
         "maximum_sampling_events": max(n_samples, 10), "y_LB": -Y_SAMPLE, "y_RB": Y_SAMPLE,
         "calculate_vn": 0, "calculate_dN_deta": 0, "calculate_dN_dxt": 0, "calculate_dN_dphi": 0,
         "store_samples_in_memory": 1, "use_binary_format": 1, "use_OSCAR_format": 0,
         "use_gzip_format": 0, "output_samples_into_files": 0, "randomSeed": seed}
    return [f"{k}={v}" for k, v in p.items()]


def write_initial(path, s, x):
    n = len(x)
    X, Y = np.meshgrid(x, x, indexing="ij")
    cols = np.zeros((n * n, 11))
    cols[:, 1], cols[:, 2], cols[:, 3], cols[:, 4] = X.ravel(), Y.ravel(), s.ravel(), 1.0
    with open(path, "w") as f:
        f.write(f"# dummy 1 etamax= 1 xmax= {n} ymax= {n} deta= 0 dx= {x[1]-x[0]:.6g} dy= {x[1]-x[0]:.6g}\n")
        np.savetxt(f, cols, fmt="%.6g")


REC = np.dtype([("pid", "<i4"), ("m", "<f4"), ("t", "<f4"), ("x", "<f4"), ("y", "<f4"), ("z", "<f4"),
                ("E", "<f4"), ("px", "<f4"), ("py", "<f4"), ("pz", "<f4")])


def read_samples(path):
    """iSS binary output: per sample an int32 count, then count x (int32 pid, 9 float32)."""
    buf = Path(path).read_bytes()
    out, pos = [], 0
    while pos < len(buf):
        n = int(np.frombuffer(buf, "<i4", 1, pos)[0]); pos += 4
        out.append(np.frombuffer(buf, REC, n, pos)); pos += n * REC.itemsize
    return out


def kinematics(a):
    px, py, pz, E = (a[k].astype(np.float64) for k in ("px", "py", "pz", "E"))
    pt = np.hypot(px, py); p = np.sqrt(pt ** 2 + pz ** 2)
    eta = 0.5 * np.log((p + pz + 1e-12) / (p - pz + 1e-12))
    y = 0.5 * np.log((E + pz) / (E - pz))
    return pt, np.arctan2(py, px), eta, y


def event_observables(samples):
    """Per-hydro-event sums (combined over events later with pair weights)."""
    nb = len(PT_EDGES) - 1
    o = {"n_samples": len(samples), "nch_eta05": 0.0, "unknown_pid": 0,
         "spec": {k: np.zeros(nb) for k in SPECIES}, "sumpt": {k: 0.0 for k in SPECIES},
         "npt": {k: 0.0 for k in SPECIES}, "spec_ch": np.zeros(nb)}
    Qs = {n: [] for n in (2, 3)}; Ns = []
    poi = {(n, s): [] for n in (2, 3) for s in list(SPECIES) + ["ch"]}   # per sample: (sum_q per bin, count per bin)
    for a in samples:
        pt, phi, eta, y = kinematics(a)
        q = charge(a["pid"])
        o["unknown_pid"] += int(np.sum(q == 99))
        ch = (q != 0) & (q != 99)
        o["nch_eta05"] += np.sum(ch & (np.abs(eta) < 0.5))
        o["spec_ch"] += np.histogram(pt[ch & (np.abs(eta) < 0.5)], PT_EDGES)[0]
        for k, pid in SPECIES.items():
            m = (np.abs(a["pid"]) == pid) & (np.abs(y) < 0.5)
            o["spec"][k] += np.histogram(pt[m], PT_EDGES)[0]
            o["sumpt"][k] += pt[m].sum(); o["npt"][k] += m.sum()
        ref = ch & (np.abs(eta) < ETA_FLOW) & (pt > PT_FLOW[0]) & (pt < PT_FLOW[1])
        Ns.append(ref.sum())
        for n in (2, 3):
            Qs[n].append(np.sum(np.exp(1j * n * phi[ref])))
            for s in list(SPECIES) + ["ch"]:
                sel = (ch if s == "ch" else (np.abs(a["pid"]) == SPECIES[s])) & (np.abs(eta) < ETA_FLOW)
                idx = np.digitize(pt[sel], PT_EDGES) - 1
                ok = (idx >= 0) & (idx < nb)
                qb = np.bincount(idx[ok], weights=np.cos(n * phi[sel][ok]), minlength=nb) \
                    + 1j * np.bincount(idx[ok], weights=np.sin(n * phi[sel][ok]), minlength=nb)
                poi[(n, s)].append((qb, np.bincount(idx[ok], minlength=nb)))
    M = len(samples)
    o["nch_eta05"] /= M
    Ns = np.array(Ns, float)
    for n in (2, 3):
        Q = np.array(Qs[n]); Qt = Q.sum()
        o[f"ref{n}_num"] = float(abs(Qt) ** 2 - np.sum(abs(Q) ** 2))
        o[f"ref{n}_den"] = float(Ns.sum() ** 2 - np.sum(Ns ** 2))
        for s in list(SPECIES) + ["ch"]:
            num, den = np.zeros(nb), np.zeros(nb)
            for k, (qb, mb) in enumerate(poi[(n, s)]):
                Qo = Qt - Q[k]                     # reference flow vector from the OTHER samples
                num += np.real(qb * np.conj(Qo)); den += mb * (Ns.sum() - Ns[k])
            o[f"poi{n}_{s}_num"], o[f"poi{n}_{s}_den"] = num, den
    for k in SPECIES:
        o["spec"][k] /= M
    o["spec_ch"] /= M
    return o


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, np.ndarray):
        return [float(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    return o


def run_event(job):
    var, i, s_factor, n_samples, threads, keep = job
    ev = glauber_event(i)
    x = grid()
    s = entropy_profile(ev, x)
    e2, psi2 = eccentricity(s, x, 2); e3, psi3 = eccentricity(s, x, 3)
    d = WORK / var / f"ev{i}"
    shutil.rmtree(d, ignore_errors=True); (d / "results").mkdir(parents=True)
    write_initial(d / "initial.dat", s, x)
    (d / "music_input").write_text(music_input(var, s_factor))
    env = dict(os.environ, HYDROPROGRAMPATH=str(MUSIC_DIR), OMP_NUM_THREADS=str(threads))
    t0 = time.time()
    with open(d / "music.log", "w") as log:
        r = subprocess.run([str(MUSIC_EXE), "music_input"], cwd=d, env=env, stdout=log, stderr=subprocess.STDOUT)
    t_music = time.time() - t0
    surf = sorted(d.glob("surface_eps_*.dat"))
    if r.returncode != 0 or not surf:
        tail = (d / "music.log").read_text().splitlines()[-15:]
        raise RuntimeError(f"MUSIC failed for event {i} ({var}), exit {r.returncode}:\n" + "\n".join(tail))
    shutil.move(str(surf[0]), d / "results/surface.dat")
    shutil.copy(d / "music_input", d / "results/music_input")
    (d / "iSS_tables").symlink_to(ISS_DIR / "iSS_tables")
    t0 = time.time()
    with open(d / "iss.log", "w") as log:
        r = subprocess.run([str(ISS_EXE), str(ISS_DIR / "iSS_parameters.dat"), "results", "surface.dat"]
                           + iss_args(var, n_samples, SEED0 % 100000 + i), cwd=d, env=env,
                           stdout=log, stderr=subprocess.STDOUT)
    t_iss = time.time() - t0
    if r.returncode != 0 or not (d / "particle_samples.bin").exists():
        tail = (d / "iss.log").read_text().splitlines()[-15:]
        raise RuntimeError(f"iSS failed for event {i} ({var}), exit {r.returncode}:\n" + "\n".join(tail))
    obs = event_observables(read_samples(d / "particle_samples.bin"))
    if not keep:
        shutil.rmtree(d, ignore_errors=True)
    rec = {"variant": var, "i": i, "b": ev["b"], "npart": ev["npart"], "ncoll": ev["ncoll"],
           "S0": float(s.sum() * (x[1] - x[0]) ** 2), "eps2": e2, "psi2": psi2, "eps3": e3, "psi3": psi3,
           "s_factor": s_factor, "t_music_s": round(t_music, 1), "t_iss_s": round(t_iss, 1)}
    rec.update(obs)
    return _jsonable(rec)


# ===========================================================================
# driver: calibrate, run
# ===========================================================================
def calib_file():
    return TAB / "m10_calibration.json"


def load_calib():
    f = calib_file()
    return json.loads(f.read_text()) if f.exists() else {}


def cmd_calibrate(a):
    """s_factor such that <dN_ch/deta> of the first n 0-5 % events = 1943.
    Multiplicity ~ entropy ~ s_factor (exact for ideal hydro; viscous entropy
    production makes it slightly non-linear), so iterate the linear update."""
    cen = load_centrality()
    ids = cen["classes"][0]["events"][:a.n]
    cal = load_calib()
    sf = a.start or cal.get(a.variant, {}).get("s_factor", 30.0)
    hist = []
    for it in range(a.iter):
        with Pool(min(a.jobs, len(ids))) as p:
            recs = p.map(run_event, [(a.variant, i, sf, a.samples, a.threads, False) for i in ids])
        dn = float(np.mean([r["nch_eta05"] for r in recs]))
        hist.append({"s_factor": sf, "dndeta": dn, "t_music_s": [r["t_music_s"] for r in recs],
                     "t_iss_s": [r["t_iss_s"] for r in recs]})
        print(f"  [{a.variant}] iteration {it}: s_factor = {sf:.4g} -> dN_ch/deta(0-5%) = {dn:.1f}"
              f"  (MUSIC {np.mean(hist[-1]['t_music_s']):.0f} s, iSS {np.mean(hist[-1]['t_iss_s']):.0f} s per event)")
        sf_new = sf * DNDETA_0_5[0] / dn
        if abs(dn / DNDETA_0_5[0] - 1) < 0.01:
            break
        sf = sf_new
    cal = load_calib()
    cal[a.variant] = {"s_factor": sf, "target": DNDETA_0_5, "events": ids, "iterations": hist}
    calib_file().write_text(json.dumps(cal, indent=1))
    print(f"  [{a.variant}] s_factor = {sf:.4g}  -> {calib_file().relative_to(ROOT)}")


def events_file(var):
    return TAB / f"m10_events_{var}.jsonl"


def cmd_run(a):
    cen = load_centrality()
    cal = load_calib()
    if a.variant not in cal:
        sys.exit(f"no s_factor for {a.variant}: run `calibrate --variant {a.variant}` first")
    sf = cal[a.variant]["s_factor"]
    f = events_file(a.variant)
    done = set()
    if f.exists():
        old = [json.loads(l) for l in f.read_text().splitlines() if l.strip()]
        keep = [r for r in old if abs(r["s_factor"] / sf - 1) < 1e-9]
        if len(keep) < len(old):      # calibration changed: those events are no longer valid
            print(f"  dropping {len(old) - len(keep)} events made with another s_factor")
            f.write_text("".join(json.dumps(r) + "\n" for r in keep))
        done = {r["i"] for r in keep}
    jobs = []
    s_central = cen["classes"][0]["S"]
    for c in cen["classes"]:
        if a.classes and c["class"] not in a.classes:
            continue
        # fewer particles per sample in peripheral events: more samples (iSS is cheap)
        ns = int(round(a.samples * min(8.0, s_central / c["S"])))
        for i in c["events"][:a.per_class]:
            if i not in done:
                jobs.append((a.variant, i, sf, ns, a.threads, False))   # central (longest) first
    print(f"[{a.variant}] {len(jobs)} hydro events to run ({len(done)} already done), s_factor = {sf:.4g}")
    t0 = time.time()
    with Pool(a.jobs) as p, open(f, "a") as out:
        for k, rec in enumerate(p.imap_unordered(run_event, jobs)):
            out.write(json.dumps(rec) + "\n"); out.flush()
            if (k + 1) % max(1, len(jobs) // 20) == 0 or k + 1 == len(jobs):
                el = time.time() - t0
                print(f"  {k+1}/{len(jobs)} events, {el/60:.1f} min elapsed, ~{el/(k+1)*(len(jobs)-k-1)/60:.0f} min left",
                      flush=True)


# ===========================================================================
# analysis
# ===========================================================================
def load_events(var):
    f = events_file(var)
    if not f.exists():
        return []
    return [json.loads(l) for l in f.read_text().splitlines() if l.strip()]


def class_of(cen):
    m = {}
    for c in cen["classes"]:
        for i in c["events"]:
            m[i] = c["class"]
    return m


def summarize(var, cen):
    evs = load_events(var)
    cmap = class_of(cen)
    out = []
    for c in cen["classes"]:
        E = [e for e in evs if cmap.get(e["i"]) == c["class"]]
        if not E:
            continue
        r = {"class": c["class"], "n_hydro_events": len(E), "n_samples": int(sum(e["n_samples"] for e in E)),
             "npart": float(np.mean([e["npart"] for e in E])), "b": float(np.mean([e["b"] for e in E])),
             "eps2_rms": float(np.sqrt(np.mean([e["eps2"] ** 2 for e in E]))),
             "eps3_rms": float(np.sqrt(np.mean([e["eps3"] ** 2 for e in E])))}
        dn = np.array([e["nch_eta05"] for e in E])
        r["dndeta"], r["dndeta_err"] = float(dn.mean()), float(dn.std(ddof=1) / np.sqrt(len(dn)) if len(dn) > 1 else 0)
        for k in SPECIES:
            npt = sum(e["npt"][k] for e in E)
            r[f"meanpt_{k}"] = float(sum(e["sumpt"][k] for e in E) / npt)
            r[f"dndy_{k}"] = float(npt / r["n_samples"])
            r[f"spec_{k}"] = list(np.mean([e["spec"][k] for e in E], axis=0) / np.diff(PT_EDGES))  # dN/dy dpT
            # bootstrap over hydro events for <pT>
            rng = np.random.default_rng(1)
            bs = []
            for _ in range(200):
                j = rng.integers(0, len(E), len(E))
                bs.append(sum(E[t]["sumpt"][k] for t in j) / sum(E[t]["npt"][k] for t in j))
            r[f"meanpt_{k}_err"] = float(np.std(bs))
        for n in (2, 3):
            num = np.array([e[f"ref{n}_num"] for e in E]); den = np.array([e[f"ref{n}_den"] for e in E])
            c2 = num.sum() / den.sum()
            rng = np.random.default_rng(2)
            bs = []
            for _ in range(300):
                j = rng.integers(0, len(E), len(E))
                bs.append(np.sqrt(max(num[j].sum() / den[j].sum(), 0)))
            r[f"v{n}"], r[f"v{n}_err"] = float(np.sqrt(max(c2, 0))), float(np.std(bs))
            r[f"v{n}_event"] = [float(np.sqrt(max(a / b, 0))) for a, b in zip(num, den)]
            for s in list(SPECIES) + ["ch"]:
                pn = np.sum([e[f"poi{n}_{s}_num"] for e in E], axis=0)
                pd = np.sum([e[f"poi{n}_{s}_den"] for e in E], axis=0)
                with np.errstate(invalid="ignore", divide="ignore"):
                    r[f"v{n}pt_{s}"] = list(np.where(pd > 0, pn / pd, np.nan) / np.sqrt(c2))
        r["eps2_event"] = [e["eps2"] for e in E]; r["eps3_event"] = [e["eps3"] for e in E]
        r["t_music_s"] = float(np.mean([e["t_music_s"] for e in E]))
        r["unknown_pid"] = int(sum(e["unknown_pid"] for e in E))
        out.append(r)
    return out


def cmd_analyze(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cen = load_centrality()
    cal = load_calib()
    FIG.mkdir(parents=True, exist_ok=True)
    versions = {}
    vf = EXT / "hydro_versions.json"
    if vf.exists():
        versions = json.loads(vf.read_text())
    S = {}
    for var in VARIANTS:
        rows = summarize(var, cen)
        if not rows:
            continue
        S[var] = rows
        res = {"variant": var, "settings": VARIANTS[var] | {"tau0": TAU0, "e_sw": E_SW, "w": W_NUC,
                                                             "k_gamma": K_GAMMA, "grid": [X_MAX, DX]},
               "s_factor": cal.get(var, {}).get("s_factor"), "codes": versions, "classes": rows}
        (TAB / f"m10_{var}.json").write_text(json.dumps(res, indent=1))
        print(f"\n[{var}] s_factor = {res['s_factor']:.4g}")
        print("  class    nev  Npart  dNch/deta   <pT> pi    K      p     eps2   v2{2}          eps3   v3{2}")
        for r in rows:
            print(f"  {r['class']:>7} {r['n_hydro_events']:4d} {r['npart']:6.1f} {r['dndeta']:8.1f}   "
                  f"{r['meanpt_pi']:.3f} {r['meanpt_K']:.3f} {r['meanpt_p']:.3f}  {r['eps2_rms']:.3f} "
                  f"{r['v2']:.4f}+-{r['v2_err']:.4f}  {r['eps3_rms']:.3f} {r['v3']:.4f}+-{r['v3_err']:.4f}")
    if not S:
        sys.exit("no events: run `run` first")
    ang = json.loads((ROOT / "results/tables/m07/m07_PbPb.json").read_text())
    ang = {c["class"]: c for c in ang["classes"]}
    cls_mid = {c["class"]: 0.5 * sum(map(float, c["class"].rstrip("%").split("-"))) for c in cen["classes"]}
    colors = {"nominal": "C0", "ideal": "C3", "etas020": "C2", "bulk": "C1"}

    # --- initial conditions
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.9))
    x = grid()
    for ax, ci in zip(axs, (0, 3, 5)):
        c = cen["classes"][ci]; ev = glauber_event(c["events"][0]); s = entropy_profile(ev, x)
        e2, _ = eccentricity(s, x, 2); e3, _ = eccentricity(s, x, 3)
        ax.imshow(s.T, origin="lower", extent=[-X_MAX, X_MAX, -X_MAX, X_MAX], cmap="inferno")
        ax.set_xlim(-10, 10); ax.set_ylim(-10, 10); ax.set_xlabel("x (fm)")
        ax.set_title(f"{c['class']}: b = {ev['b']:.1f} fm, " + r"$\varepsilon_2$" + f" = {e2:.2f}, "
                     + r"$\varepsilon_3$" + f" = {e3:.2f}", fontsize=9)
    axs[0].set_ylabel("y (fm)")
    fig.suptitle(r"initial entropy density $\propto\sqrt{T_A T_B}$ (Module 6 Glauber, $w$ = 0.5 fm, $\Gamma$ weights)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "m10_initial.png", dpi=130); plt.close(fig)

    # --- multiplicity and <pT>
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axs[0]
    for var, rows in S.items():
        ax.errorbar([r["npart"] for r in rows], [r["dndeta"] / (r["npart"] / 2) for r in rows],
                    yerr=[r["dndeta_err"] / (r["npart"] / 2) for r in rows], fmt="o-", ms=4, color=colors[var],
                    label="MUSIC " + VARIANTS[var]["label"])
    ax.plot([ang[c]["npart"] for c in ang if c in cls_mid], [ang[c]["dndeta"] / (ang[c]["npart"] / 2) for c in ang if c in cls_mid],
            "s--", color="gray", ms=4, label="Angantyr (Module 7)")
    n05 = S["nominal"][0]["npart"] if "nominal" in S else 383
    ax.errorbar([n05], [DNDETA_0_5[0] / (n05 / 2)], yerr=[DNDETA_0_5[1] / (n05 / 2)], fmt="D", color="k",
                label="ALICE 0-5% (calibration point)")
    ax.set_xlabel(r"$\langle N_{\rm part}\rangle$"); ax.set_ylabel(r"$\frac{dN_{\rm ch}/d\eta}{\langle N_{\rm part}\rangle/2}$")
    ax.legend(fontsize=7); ax.set_title(r"multiplicity per participant pair, $|\eta|<0.5$", fontsize=10)
    ax = axs[1]
    for k, mk in zip(SPECIES, "osD"):
        for var in [v for v in ("nominal", "bulk") if v in S]:
            rows = S[var]
            ax.errorbar([cls_mid[r["class"]] for r in rows], [r[f"meanpt_{k}"] for r in rows],
                        yerr=[r[f"meanpt_{k}_err"] for r in rows], fmt=mk + ("-" if var == "nominal" else ":"),
                        color=colors[var], ms=4, label=f"{k}, MUSIC {VARIANTS[var]['label']}")
        ax.plot([cls_mid[c] for c in cls_mid if c in ang], [ang[c]["meanpt"][k] for c in cls_mid if c in ang],
                mk + "--", color="gray", ms=4, mfc="none", label=f"{k}, Angantyr")
    ax.set_xlabel("centrality (%)"); ax.set_ylabel(r"$\langle p_{\rm T}\rangle$ (GeV/$c$), $|y|<0.5$")
    ax.legend(fontsize=6, ncol=2); ax.set_title("radial flow: mass ordering of $\\langle p_T\\rangle$", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "m10_multiplicity_meanpt.png", dpi=130); plt.close(fig)

    # --- spectra 0-5 %
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    ptc = 0.5 * (PT_EDGES[1:] + PT_EDGES[:-1])
    a05 = ang.get("0-5%")
    for k, col in zip(SPECIES, ("C0", "C1", "C2")):
        if "nominal" in S:
            ax.semilogy(ptc, S["nominal"][0][f"spec_{k}"], "-", color=col, label=f"{k}$^+$+{k}$^-$ MUSIC")
        if a05 is not None and "counts" in a05:
            pass
    ax.set_xlabel(r"$p_{\rm T}$ (GeV/$c$)"); ax.set_ylabel(r"$dN/dy\,dp_{\rm T}$ (GeV$^{-1}c$), $|y|<0.5$")
    ax.set_title("0-5% Pb-Pb, identified hadrons (MUSIC + iSS, decays, no afterburner)", fontsize=9)
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(FIG / "m10_spectra.png", dpi=130); plt.close(fig)

    # --- flow vs centrality, and response
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axs[0]
    for var, rows in S.items():
        xc = [cls_mid[r["class"]] for r in rows]
        ax.errorbar(xc, [r["v2"] for r in rows], yerr=[r["v2_err"] for r in rows], fmt="o-", ms=4,
                    color=colors[var], label=r"$v_2\{2\}$ " + VARIANTS[var]["label"])
        ax.errorbar(xc, [r["v3"] for r in rows], yerr=[r["v3_err"] for r in rows], fmt="s:", ms=4,
                    color=colors[var], mfc="none")
    ax.set_xlabel("centrality (%)"); ax.set_ylabel(r"$v_n\{2\}$, $0.2<p_{\rm T}<3$ GeV/$c$, $|\eta|<0.8$")
    ax.set_title(r"elliptic (filled) and triangular (open) flow", fontsize=10); ax.legend(fontsize=7)
    ax = axs[1]
    for var, rows in S.items():
        e = np.concatenate([r["eps2_event"] for r in rows]); v = np.concatenate([r["v2_event"] for r in rows])
        k2 = float(np.sum(e * v) / np.sum(e * e))
        ax.scatter(e, v, s=6, color=colors[var], alpha=0.5, label=VARIANTS[var]["label"] + f": $v_2/\\varepsilon_2$ = {k2:.2f}")
        for r in rows:
            r["kappa2"] = k2
    xx = np.linspace(0, 0.7, 10)
    ax.set_xlabel(r"$\varepsilon_2$ (initial state)"); ax.set_ylabel(r"$v_2$ of the hydro event")
    ax.set_title("linear response: $v_2 \\approx \\kappa_2\\,\\varepsilon_2$", fontsize=10); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "m10_flow.png", dpi=130); plt.close(fig)

    # --- v2(pT) mass ordering
    if "nominal" in S:
        rows = {r["class"]: r for r in S["nominal"]}
        cl = "30-40%" if "30-40%" in rows else list(rows)[-1]
        fig, ax = plt.subplots(figsize=(5.6, 4.4))
        for s, col in zip(list(SPECIES) + ["ch"], ("C0", "C1", "C2", "k")):
            y = np.array(rows[cl][f"v2pt_{s}"])
            ax.plot(ptc, y, "o-" if s != "ch" else "--", ms=3, color=col, label=s if s != "ch" else "charged")
        ax.set_xlim(0, 3); ax.set_ylim(0, None)
        ax.set_xlabel(r"$p_{\rm T}$ (GeV/$c$)"); ax.set_ylabel(r"$v_2\{2\}(p_{\rm T})$, $|\eta|<0.8$")
        ax.set_title(f"{cl}, MUSIC {VARIANTS['nominal']['label']}: mass ordering", fontsize=10); ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(FIG / "m10_v2pt.png", dpi=130); plt.close(fig)
    print(f"\nfigures -> {FIG.relative_to(ROOT)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("scan"); p.add_argument("--n", type=int, default=4000); p.add_argument("--jobs", type=int, default=os.cpu_count())
    for name in ("calibrate", "run"):
        p = sub.add_parser(name)
        p.add_argument("--variant", default="nominal", choices=list(VARIANTS))
        p.add_argument("--samples", type=int, default=100, help="iSS oversamples per hydro event")
        p.add_argument("--jobs", type=int, default=os.cpu_count(), help="hydro events in parallel")
        p.add_argument("--threads", type=int, default=1, help="OpenMP threads per MUSIC run")
    sub.choices["calibrate"].add_argument("--n", type=int, default=4)
    sub.choices["calibrate"].add_argument("--iter", type=int, default=3)
    sub.choices["calibrate"].add_argument("--start", type=float, default=None)
    sub.choices["run"].add_argument("--per-class", type=int, default=20)
    sub.choices["run"].add_argument("--classes", nargs="*", default=None)
    sub.add_parser("analyze")
    a = ap.parse_args()
    {"scan": cmd_scan, "calibrate": cmd_calibrate, "run": cmd_run, "analyze": cmd_analyze}[a.cmd](a)


if __name__ == "__main__":
    main()
