#!/usr/bin/env python3
"""Module 9: jet quenching with Glauber geometry -- charged-hadron R_AA and high-pT v2.

    python analysis/m09_quenching.py pp      [--n-events 400000] [--workers N]   # PYTHIA pp 5.02 TeV hard QCD
    python analysis/m09_quenching.py glauber [--n-events 30000]                   # Pb-Pb geometry pool (Module 6)
    python analysis/m09_quenching.py model                                         # calibrate, predict, figures

QUESTION
--------
In central Pb-Pb collisions at 5.02 TeV, charged hadrons at pT 6-9 GeV/c are
7-8 times fewer than N_coll x pp (CMS, arXiv:1611.01664: R_AA ~ 0.13-0.14 in
0-5 %).  Angantyr (Module 7), which has no medium, gives R_AA ~ 1.  The
simplest explanation: the scattered partons lose energy crossing the dense
matter made by the participants.  This module builds that explanation from
pieces made earlier and asks what ONE number (the R_AA dip) predicts for
everything else.

MODEL (deliberately minimal, every piece visible)
-------------------------------------------------
1. pp reference: PYTHIA 8 hard QCD 2 -> 2 at 5.02 TeV in pT-hat slices with
   slice weights sigma/N (generators/pythia_hardqcd_5TeV.cmnd).  Per event the
   two outgoing partons and every charged hadron with |eta| < 1, pT > 3 GeV/c,
   each hadron assigned to the nearer parton in Delta R.  (A first version
   left hadrons with Delta R > 0.5 untouched: 30-40 % of the 6-9 GeV/c
   hadrons -- shower, ISR, MPI -- then set a floor R_AA ~ 0.5 that no kappa
   could reach below.  In a medium they lose energy too.)
2. Geometry: Module 6 Monte Carlo Glauber Pb-Pb events in centrality classes
   (impact-parameter percentiles).  A hard scattering happens at a random
   nucleon-nucleon COLLISION point (density ~ N_coll); the medium density is
   the sum of Gaussians (width 0.5 fm) on the WOUNDED nucleons (~ N_part).
3. Energy loss of each parton along its direction:
       Delta E = kappa * I_zeta * (E / 10 GeV)^beta,
       I_zeta = integral_0^inf rho(x0 + t n) t^zeta dt
   zeta = 0: loss ~ path length L (collisional-like); zeta = 1: ~ L^2
   (radiative, BDMPS-like, static medium); beta: energy dependence.
4. Hadrons of a jet keep their momentum FRACTION: pT -> pT (1 - Delta E / E).
5. R_AA(pT) = quenched / unquenched spectrum of the SAME events (N_coll
   scaling is built in: hard processes are placed with the collision
   density).  v2(pT): every event is also quenched with its partons turned
   IN the participant plane Psi_2 (Module 6 geometry) and OUT of it;
   R(dphi) = R0 (1 + 2 v2 cos 2 dphi) gives v2 = (R_in - R_out)/(2 (R_in + R_out)).
   (The direct estimator <cos 2(phi - Psi_2)> with random orientations is kept
   for comparison; with 400k weighted events its noise was +-0.02.)
kappa is fixed for each (zeta, beta) by R_AA(0-5 %, 6-9 GeV/c) = 0.133
(the CMS factor 7.5); everything else is a prediction.

NOT IN THE MODEL: expansion and time dependence of the medium, fluctuations of
the energy loss, medium response, nuclear PDFs (Module 5: small at x ~ 0.01),
the change of the jet's hadron composition.  It is a geometry exercise.

OUTPUT
  data/generated/m09_pp.npz, m09_glauber.npz          (not in git)
  the n_eff of the weighted hadron sample is printed (pp step)
  results/tables/m09/m09_model.json; results/figures/m09/m09_raa.png, m09_v2.png, m09_raa_npart.png
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from hicpy.glauber import PB208, sample_nucleus, event_geometry, participant_plane, path_integrals  # noqa: E402

GEN, TAB, FIG = REPO / "data" / "generated", REPO / "results" / "tables" / "m09", REPO / "results" / "figures" / "m09"
CARD = REPO / "generators" / "pythia_hardqcd_5TeV.cmnd"
SLICES = [8, 12, 18, 27, 40, 60, 100, -1]                # pT-hat slice edges (GeV); -1 = no upper limit
PT_SWITCH = 10.0             # GeV/c: minimum bias below, hard slices above (unquenched hadron pT)
SIGMA_NN = 67.6
CLASSES = [(0, 5), (5, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 70), (70, 80)]
PT_EDGES = np.array([5, 6, 7, 8, 9, 10, 12, 14, 17, 20, 25, 30, 40, 50, 70, 100], float)
CALIB = {"class": "0-5%", "pt": (6.0, 9.0), "raa": 0.133,
         "source": "CMS, arXiv:1611.01664: maximal suppression by a factor 7-8 at pT 6-9 GeV, 0-5 %"}
MODELS = {"L, E^0": (0, 0.0), "L^2, E^0": (1, 0.0), "L, E^1/3": (0, 1 / 3), "L^2, E^1/3": (1, 1 / 3)}
WIDTH = 0.5


# ----------------------------------------------------------------------------- pp
def _pp_chunk(job):
    """One generation job.  kind 'mb': SoftQCD:inelastic (unweighted, MPI-aware:
    the right reference below ~12 GeV/c); kind 'hard': HardQCD in one pT-hat
    slice.  Returns events [weight, sample], partons [event, pT, eta, phi] (all
    outgoing partons of the hard process AND of the MPI scatterings, status
    23 and 33, pT > 1 GeV), hadrons [event, pT, eta, phi, parton, Delta R]."""
    kind, n_events, seed, lo, hi = job
    try:
        import pythia8
    except ImportError:
        import pythia8mc as pythia8
    py = pythia8.Pythia("", False)
    py.readFile(str(CARD))
    settings = [f"Random:seed = {seed}"]
    if kind == "mb":
        settings += ["HardQCD:all = off", "SoftQCD:inelastic = on"]
    else:
        settings += [f"PhaseSpace:pTHatMin = {lo}", f"PhaseSpace:pTHatMax = {hi}"]
    for x in settings:
        py.readString(x)
    if not py.init():
        raise RuntimeError("PYTHIA init() failed")
    ev_rows, par_rows, had_rows = [], [], []
    n = 0
    while n < n_events:
        if not py.next():
            continue
        ev = py.event
        partons = [(ev[i].pT(), ev[i].eta(), ev[i].phi()) for i in range(ev.size())
                   if abs(ev[i].status()) in (23, 33) and ev[i].pT() > 1.0]
        first = len(par_rows)
        for q in partons:
            par_rows.append([n, *q])
        for i in range(ev.size()):
            p = ev[i]
            if not (p.isFinal() and p.isCharged()) or p.pT() < 3.0 or abs(p.eta()) > 1.0:
                continue
            if partons:
                dr = [np.hypot(p.eta() - q[1], (p.phi() - q[2] + np.pi) % (2 * np.pi) - np.pi) for q in partons]
                k = int(np.argmin(dr)); had_rows.append([n, p.pT(), p.eta(), p.phi(), first + k, dr[k]])
            else:                                   # (diffractive event without a hard parton: no loss)
                had_rows.append([n, p.pT(), p.eta(), p.phi(), -1, 99.0])
        ev_rows.append([1.0, 0 if kind == "mb" else 1])
        n += 1
    # infoPython() is a partial SNAPSHOT: weight() follows the events, but
    # sigmaGen() and nAccepted() keep the values of the moment it was called
    # (0 before generation) -> ask for it now (troubleshooting #29)
    info = py.infoPython() if hasattr(py, "infoPython") else py.info
    return (np.array(ev_rows), np.array(par_rows).reshape(-1, 4), np.array(had_rows).reshape(-1, 6),
            py.settings.parm("Pythia:versionNumber"), (kind, lo, hi, info.sigmaGen(), n))


def make_pp(args) -> int:
    w = args.workers or os.cpu_count() or 1
    n_sl = len(SLICES) - 1
    n_mb = args.n_events // 2                         # half minimum bias, half in the hard slices
    jobs = [("mb", n_mb // w, 20261122 + i, 0, 0) for i in range(w)]
    per_job = max(1, (args.n_events - n_mb) // (n_sl * w))
    jobs += [("hard", per_job, 20261222 + 17 * (i * w + k), SLICES[i], SLICES[i + 1]) for i in range(n_sl) for k in range(w)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=w) as ex:
        parts = list(ex.map(_pp_chunk, jobs))
    # weight per event = (cross section of its sample or slice) / (events in it)
    sig = {}
    for *_, (kind, lo, hi, s, n) in parts:
        a = sig.setdefault((kind, lo, hi), [0.0, 0]); a[0] += s * n; a[1] += n
    evs, pars, hads, off_e, off_p = [], [], [], 0, 0
    for e, pa, h, _, (kind, lo, hi, s, n) in parts:
        tot_s, tot_n = sig[(kind, lo, hi)]
        if len(e):
            e[:, 0] = (tot_s / tot_n) / tot_n
        pa = pa.copy(); pa[:, 0] += off_e
        h = h.copy(); h[:, 0] += off_e; h[:, 4] = np.where(h[:, 4] >= 0, h[:, 4] + off_p, -1)
        evs.append(e); pars.append(pa); hads.append(h); off_e += len(e); off_p += len(pa)
    ev, par, had = np.vstack(evs), np.vstack(pars), np.vstack(hads)
    GEN.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(GEN / "m09_pp.npz", events=ev, partons=par, hadrons=had, version=parts[0][3])
    for (kind, lo, hi), (s, n) in sorted(sig.items(), key=lambda x: (x[0][0] != "mb", x[0][1])):
        name = "minimum bias" if kind == "mb" else f"pT-hat {lo}-{'inf' if hi < 0 else hi} GeV"
        print(f"[sample] {name:22s}: sigma {s / n:.4g} mb, {n} events")
    use = domain_mask(ev, had)
    for lo, hi in ((6, 9), (10, 30), (40, 70)):
        sel = use & (had[:, 1] >= lo) & (had[:, 1] < hi)
        _, inv = np.unique(had[sel, 0], return_inverse=True)
        we = np.bincount(inv, weights=ev[had[sel, 0].astype(int), 0])
        print(f"[n_eff] hadrons {lo}-{hi} GeV/c: {len(we)} events, effective {we.sum() ** 2 / np.sum(we ** 2):.0f}")
    # continuity check where the two samples meet
    for k, name in ((0, "minimum bias"), (1, "hard slices")):
        sel = (ev[had[:, 0].astype(int), 1] == k) & (had[:, 1] >= PT_SWITCH) & (had[:, 1] < PT_SWITCH + 4)
        print(f"[match] {name:12s}: yield {PT_SWITCH:.0f}-{PT_SWITCH + 4:.0f} GeV/c = {np.sum(ev[had[sel, 0].astype(int), 0]):.4g} mb "
              f"({sel.sum()} hadrons)")
    print(f"[pp] {len(ev)} events, {len(par)} partons, {len(had)} charged hadrons (|eta|<1, pT>3); "
          f"PYTHIA {parts[0][3]:.3f}, {time.time() - t0:.0f} s")
    return 0


def domain_mask(ev, had):
    """Which sample describes a hadron: minimum bias below PT_SWITCH (UNQUENCHED
    pT), the hard slices above.  Each hadron is counted from one sample only."""
    sample = ev[had[:, 0].astype(int), 1]
    return np.where(had[:, 1] < PT_SWITCH, sample == 0, sample == 1)


# ----------------------------------------------------------------------------- glauber
def _glauber_chunk(job):
    n, seed = job
    rng = np.random.default_rng(seed)
    out = []
    bmax = 20.0
    for _ in range(n):
        b = bmax * np.sqrt(rng.random())
        coll, wounded = event_geometry(sample_nucleus(rng, PB208), sample_nucleus(rng, PB208), b, SIGMA_NN)
        if len(coll) == 0:
            continue
        psi, eps = participant_plane(wounded)
        out.append((b, coll.astype(np.float32), wounded.astype(np.float32), psi, eps))
    return out


def _objarr(items):
    """1D object array of variable-length arrays (np.array would try to stack them)."""
    out = np.empty(len(items), dtype=object)
    for i, x in enumerate(items):
        out[i] = x
    return out


def make_glauber(args) -> int:
    w = args.workers or os.cpu_count() or 1
    jobs = [(args.n_events // w + (i < args.n_events % w), 900 + i) for i in range(w)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=w) as ex:
        evs = [e for part in ex.map(_glauber_chunk, jobs) for e in part]
    b = np.array([e[0] for e in evs]); order = np.argsort(b); n = len(b)
    cls = np.full(n, -1)
    for k, (lo, hi) in enumerate(CLASSES):
        cls[order[int(n * lo / 100):int(n * hi / 100)]] = k
    np.savez_compressed(GEN / "m09_glauber.npz", b=b, cls=cls,
                        npart=np.array([len(e[2]) for e in evs]), ncoll=np.array([len(e[1]) for e in evs]),
                        psi2=np.array([e[3] for e in evs]), eps2=np.array([e[4] for e in evs]),
                        coll=_objarr([e[1] for e in evs]), wounded=_objarr([e[2] for e in evs]), allow_pickle=True)
    print(f"[glauber] {n} Pb-Pb events with a collision ({args.n_events} thrown), {time.time() - t0:.0f} s")
    for k, (lo, hi) in enumerate(CLASSES):
        sel = cls == k
        print(f"  {lo:2d}-{hi:2d}%: {sel.sum():5d} events, <Npart> {np.mean([len(evs[i][2]) for i in np.nonzero(sel)[0]]):6.1f}, "
              f"<eps2> {np.mean([evs[i][4] for i in np.nonzero(sel)[0]]):.3f}")
    return 0


# ----------------------------------------------------------------------------- model
def path_lengths(pp, gl, cls_index, rng, n_use):
    """Hadron-level arrays for one centrality class.

    n_use pp events of each sample are placed in Glauber events of the class,
    at a random collision point.  I_zeta is computed for every parton that has
    a hadron of interest (pT > 5 GeV/c, in its sample's domain), in three
    orientations: 0 = as generated (random w.r.t. the Glauber event), 1 = the
    event turned so that its leading parton points along Psi_2 (in plane),
    2 = along Psi_2 + pi/2 (out of plane).  In/out give v2 without azimuthal
    sampling noise."""
    ev, par, had = pp["events"], pp["partons"], pp["hadrons"]
    chosen = np.zeros(len(ev), bool)
    for k in (0, 1):
        pool = np.nonzero(ev[:, 1] == k)[0]
        chosen[rng.choice(pool, size=min(n_use, len(pool)), replace=False)] = True
    hsel = chosen[had[:, 0].astype(int)] & domain_mask(ev, had) & (had[:, 1] >= PT_EDGES[0])
    h = had[hsel]
    e_of = h[:, 0].astype(int)
    events = np.unique(e_of)
    g_pool = np.nonzero(gl["cls"] == cls_index)[0]
    g_of_ev = dict(zip(events, rng.choice(g_pool, size=len(events))))
    pj = h[:, 4].astype(int)
    has_p = pj >= 0
    plist = np.unique(pj[has_p])                                   # partons that matter
    p_ev = par[plist, 0].astype(int)
    # leading parton of each event (for the in/out orientation)
    lead_phi = {}
    order = np.lexsort((-par[:, 1], par[:, 0]))
    first = np.r_[True, par[order[1:], 0] != par[order[:-1], 0]]
    for row in order[first]:
        lead_phi[int(par[row, 0])] = par[row, 3]
    Ip = np.zeros((len(plist), 2, 3))
    psi_ev = {e: gl["psi2"][g] for e, g in g_of_ev.items()}
    x0_ev = {}
    g_of_p = np.array([g_of_ev[e] for e in p_ev])
    for g in np.unique(g_of_p):
        sel = np.nonzero(g_of_p == g)[0]
        coll, wounded = gl["coll"][g].astype(float), gl["wounded"][g].astype(float)
        for e in np.unique(p_ev[sel]):
            x0_ev[e] = coll[rng.integers(0, len(coll))]
        x0 = np.array([x0_ev[e] for e in p_ev[sel]])
        phi = par[plist[sel], 3]
        lead = np.array([lead_phi[e] for e in p_ev[sel]]); psi = gl["psi2"][g]
        for o, rot in ((0, 0.0), (1, psi - lead), (2, psi + np.pi / 2 - lead)):
            for z in (0, 1):
                Ip[sel, z, o] = path_integrals(x0, phi + rot, wounded, WIDTH, z)
    where = np.searchsorted(plist, np.where(has_p, pj, plist[0] if len(plist) else 0))
    I = np.where(has_p[:, None, None], Ip[where] if len(plist) else 0.0, 0.0)
    e_parton = np.where(has_p, par[np.where(has_p, pj, 0), 1], np.inf)
    return {"w": ev[e_of, 0], "pt": h[:, 1], "phi": h[:, 3], "psi": np.array([psi_ev[e] for e in e_of]),
            "e_parton": np.maximum(e_parton, h[:, 1]), "I": I}


def spectra(H, kappa, zeta, beta, edges=PT_EDGES, orient=0):
    """Weighted pT spectra (unquenched, quenched) and v2 numerator."""
    de = kappa * H["I"][:, zeta, orient] * (np.where(np.isfinite(H["e_parton"]), H["e_parton"], 10.0) / 10.0) ** beta
    frac = np.clip(1 - de / H["e_parton"], 0, 1)
    ptq = H["pt"] * frac
    n0 = np.histogram(H["pt"], edges, weights=H["w"])[0]
    nq = np.histogram(ptq, edges, weights=H["w"])[0]
    c2 = np.histogram(ptq, edges, weights=H["w"] * np.cos(2 * (H["phi"] - H["psi"])))[0]
    return n0, nq, c2


def calibrate(H, zeta, beta, target, pt_range):
    edges = np.array(pt_range)
    def raa(kap):                                                          # noqa: E306
        n0, nq, _ = spectra(H, kap, zeta, beta, edges)
        return nq.sum() / n0.sum()
    lo, hi = 0.0, 1.0
    while raa(hi) > target and hi < 1e3:
        hi *= 2
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if raa(mid) > target else (lo, mid)
    return 0.5 * (lo + hi)


def model(args) -> int:
    from hicpy.style import apply_style, save_figure, OKABE_ITO, plt
    z = np.load(GEN / "m09_pp.npz"); pp = {"events": z["events"], "partons": z["partons"], "hadrons": z["hadrons"]}
    gz = np.load(GEN / "m09_glauber.npz", allow_pickle=True); gl = {k: gz[k] for k in gz.files}
    rng = np.random.default_rng(3)
    t0 = time.time()
    geo = {k: path_lengths(pp, gl, k, rng, args.n_use) for k in range(len(CLASSES))}
    print(f"[geom] path integrals for {len(CLASSES)} classes x {args.n_use} events, {time.time() - t0:.0f} s")
    res = {"calibration": CALIB, "width_fm": WIDTH, "n_pp_events_per_class": args.n_use,
           "pt_edges": PT_EDGES.tolist(), "classes": [f"{lo}-{hi}%" for lo, hi in CLASSES],
           "npart": [float(gl["npart"][gl["cls"] == k].mean()) for k in range(len(CLASSES))],
           "eps2": [float(gl["eps2"][gl["cls"] == k].mean()) for k in range(len(CLASSES))], "models": {}}
    for name, (zeta, beta) in MODELS.items():
        kap = calibrate(geo[0], zeta, beta, CALIB["raa"], CALIB["pt"])
        m = {"zeta": zeta, "beta": beta, "kappa": kap, "raa": [], "raa_err": [], "v2": [], "raa_6_9": [], "raa_40_70": [],
             "v2_10_30": []}
        for k in range(len(CLASSES)):
            n0, nq, c2 = spectra(geo[k], kap, zeta, beta)
            m["raa"].append((nq / n0).tolist()); m["v2"].append((c2 / np.where(nq > 0, nq, np.nan)).tolist())
            for key, (a, b) in (("raa_6_9", (6, 9)), ("raa_40_70", (40, 70))):
                s = (PT_EDGES[:-1] >= a) & (PT_EDGES[1:] <= b)
                m[key].append(float(nq[s].sum() / n0[s].sum()))
            # v2 from in-plane / out-of-plane suppression: R(dphi) = R0 (1 + 2 v2 cos 2 dphi)
            _, n_in, _ = spectra(geo[k], kap, zeta, beta, orient=1)
            _, n_out, _ = spectra(geo[k], kap, zeta, beta, orient=2)
            m["v2"][-1] = ((n_in - n_out) / (2 * (n_in + n_out))).tolist()
            s = (PT_EDGES[:-1] >= 10) & (PT_EDGES[1:] <= 30)
            m["v2_10_30"].append(float((n_in[s].sum() - n_out[s].sum()) / (2 * (n_in[s].sum() + n_out[s].sum()))))
            m.setdefault("v2_10_30_direct", []).append(float(c2[s].sum() / nq[s].sum()))
        res["models"][name] = m
        print(f"[{name:10s}] kappa = {kap:.3f}  "
              f"R_AA 6-9: " + " ".join(f"{v:.2f}" for v in m["raa_6_9"]) +
              f" | 40-70 (0-5 %): {m['raa_40_70'][0]:.2f} | v2 10-30: " + " ".join(f"{v:.3f}" for v in m["v2_10_30"][:6]))
    TAB.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    (TAB / "m09_model.json").write_text(json.dumps(res, indent=1))

    apply_style()
    c = 0.5 * (PT_EDGES[1:] + PT_EDGES[:-1])
    cols = dict(zip(MODELS, OKABE_ITO[:4]))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), sharey=True)
    for ax, k in zip(axes, (0, 2, 4)):
        for name, m in res["models"].items():
            ax.plot(c, m["raa"][k], "-o", ms=3, color=cols[name], label=name)
        if k == 0:
            ax.errorbar([7.5], [CALIB["raa"]], xerr=[1.5], yerr=[0.01], fmt="ks", ms=6, label="CMS 6-9 GeV (calibration)")
        ax.set_xscale("log"); ax.set_title(res["classes"][k]); ax.set_xlabel(r"$p_{\rm T}$ (GeV/$c$)")
        ax.axhline(1, color="0.6", lw=0.8); ax.set_ylim(0, 1.1)
    axes[0].set_ylabel(r"$R_{\rm AA}$ (charged hadrons)"); axes[0].legend(fontsize=8, frameon=False)
    fig.tight_layout(); save_figure(fig, FIG / "m09_raa"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    for name, m in res["models"].items():
        ax.plot(res["npart"], m["raa_6_9"], "-o", ms=4, color=cols[name], label=name)
    ax.set_xlabel(r"$\langle N_{\rm part}\rangle$"); ax.set_ylabel(r"$R_{\rm AA}$, $p_{\rm T}$ 6-9 GeV/$c$")
    ax.axhline(1, color="0.6", lw=0.8); ax.legend(fontsize=8, frameon=False)
    fig.tight_layout(); save_figure(fig, FIG / "m09_raa_npart"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    mids = [0.5 * (lo + hi) for lo, hi in CLASSES]
    for name, m in res["models"].items():
        ax.plot(mids, m["v2_10_30"], "-o", ms=4, color=cols[name], label=name)
    ax.set_xlabel("centrality (%)"); ax.set_ylabel(r"$v_2$ of hadrons, $p_{\rm T}$ 10-30 GeV/$c$")
    ax.axhline(0, color="0.6", lw=0.8); ax.legend(fontsize=8, frameon=False)
    fig.tight_layout(); save_figure(fig, FIG / "m09_v2"); plt.close(fig)
    print(f"[out] {TAB / 'm09_model.json'}, {FIG}/m09_*.png")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("pp"); a.add_argument("--n-events", type=int, default=400000); a.add_argument("--workers", type=int)
    g = sub.add_parser("glauber"); g.add_argument("--n-events", type=int, default=30000); g.add_argument("--workers", type=int)
    m = sub.add_parser("model"); m.add_argument("--n-use", type=int, default=150000)
    args = ap.parse_args()
    return {"pp": make_pp, "glauber": make_glauber, "model": model}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
