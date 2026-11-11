#!/usr/bin/env python3
"""Module 2, step 3: confront the PYTHIA dN_ch/deta with a published measurement.

    python analysis/m02_compare_data.py [--config config/analysis_m02_pythia.yaml]
           [--mc data/generated/m02_pythia_mb_13TeV.root]
           [--data data/raw/alice_pp13TeV_dndeta.json] [--event-class INEL>0]
           [--figdir ...] [--tabdir ...] [--no-data]

What is compared, and why it is not the Module 1 histogram
-----------------------------------------------------------
Module 1 plots dN_ch/deta for charged particles with pT > 0.15 GeV/c, because
that is what a tracking detector sees.  Published dN_ch/deta distributions
(ALICE, CMS, ATLAS "INEL>0" or "NSD") are instead *extrapolated to pT = 0* and
refer to *primary* charged particles.  A fair comparison therefore needs

  1. the same EVENT CLASS       applied to the MC the way the experiment
                                defines it (INEL>0: >= 1 charged primary in
                                |eta| < 1; NSD: no single diffraction) --
                                computed from the particles / process code,
                                never assumed;
  2. the same PARTICLE SET      primaries (card: tau0Max = 10 mm), charge != 0,
                                NO pT cut (pt_min = 0 in the comparison block);
  3. the same BINNING           the data's eta bins, applied to the MC
                                particles directly (no rebinning of a
                                finer histogram);
  4. a RATIO panel              MC / data with statistical and systematic
                                uncertainties shown separately.

Verdict rule (docs/roadmap.md Module 2): dN_ch/deta averaged over |eta| < 0.5
within +-15 % of the measurement.  Outside that, suspect the SETUP (event
class, primary definition, charged definition, cross-section normalisation)
before blaming the tune.

Data provenance rule
--------------------
The measurement is read through its provenance JSON (data/raw/<dataset>.json,
template in data/raw/README.md), never from a bare file: the JSON names the
paper, the HEPData table, the event class and the file.  Without it the script
runs MC-only (--no-data or missing JSON) and says so.  Numbers are never typed
into this file.

HEPData formats understood: the YAML "table" export (independent_variables /
dependent_variables; canonical) and the CSV export (header comments starting
with '#', then x, x_low, x_high, y, stat+, stat-, sys+, sys- by position).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import awkward as ak
import numpy as np
import uproot
import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from hicpy.style import apply_style, stamp, save_figure, OKABE_ITO, plt   # noqa: E402
from m01_analysis_uproot import rel, hist_record                          # noqa: E402

# PYTHIA 8 process codes for SoftQCD (see the PYTHIA manual, "QCD processes")
PROCESS = {101: "ND", 102: "EL", 103: "SD(XB)", 104: "SD(AX)", 105: "DD", 106: "CD"}
SD_CODES = {103, 104}
VERDICT_TOLERANCE = 0.15


# --------------------------------------------------------------------- data
def load_hepdata(path: pathlib.Path, dependent_index: int = 0):
    """Return dict(x, xlo, xhi, y, stat, sys) as float arrays.  stat/sys are
    symmetrised: max(|plus|, |minus|) -- conservative, and the published
    asymmetries for dN/deta are small.  `dependent_index` picks the column
    when a HEPData table carries several (e.g. INEL and INEL>0)."""
    path = pathlib.Path(path)
    if path.suffix in (".yaml", ".yml"):
        doc = yaml.safe_load(path.read_text())
        xv = doc["independent_variables"][0]["values"]
        yv = doc["dependent_variables"][dependent_index]["values"]
        xlo = np.array([v["low"] if "low" in v else v["value"] for v in xv], float)
        xhi = np.array([v["high"] if "high" in v else v["value"] for v in xv], float)
        y = np.array([v["value"] for v in yv], float)
        stat = np.zeros_like(y)
        sys_ = np.zeros_like(y)
        for i, v in enumerate(yv):
            for err in v.get("errors", []):
                label = str(err.get("label", "")).lower()
                if "symerror" in err:
                    mag = abs(float(err["symerror"]))
                else:
                    a = err.get("asymerror", {})
                    mag = max(abs(float(a.get("plus", 0))), abs(float(a.get("minus", 0))))
                if "stat" in label:
                    stat[i] = np.hypot(stat[i], mag)
                else:                       # sys, syst, or unlabelled -> systematic
                    sys_[i] = np.hypot(sys_[i], mag)
    else:
        rows = []
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip().strip('"') for p in line.split(",")]
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue                    # header row
        if not rows:
            raise ValueError(f"no numeric rows in {path}")
        arr = np.array(rows, float)
        if arr.shape[1] < 4:
            raise ValueError(f"{path}: need at least x, x_low, x_high, y columns")
        xlo, xhi, y = arr[:, 1], arr[:, 2], arr[:, 3]
        stat = np.maximum(abs(arr[:, 4]), abs(arr[:, 5])) if arr.shape[1] >= 6 else np.zeros_like(y)
        sys_ = np.maximum(abs(arr[:, 6]), abs(arr[:, 7])) if arr.shape[1] >= 8 else np.zeros_like(y)
    order = np.argsort(xlo)
    return {"x": 0.5 * (xlo + xhi)[order], "xlo": xlo[order], "xhi": xhi[order],
            "y": y[order], "stat": stat[order], "sys": sys_[order]}


def load_dataset(prov_json: pathlib.Path):
    """Provenance JSON -> (provenance dict, table dict).  Refuses to load a
    table whose provenance is still the unfilled template."""
    prov = json.loads(prov_json.read_text())
    required = ["experiment", "sqrt_s_TeV", "observable", "event_class", "file"]
    missing = [k for k in required if k not in prov or str(prov[k]).startswith("<")]
    if missing:
        raise ValueError(f"{prov_json.name}: provenance incomplete, fill in {missing} "
                         "(see data/raw/README.md)")
    if not prov.get("verified", True):
        raise ValueError(f"{prov_json.name}: identifiers not yet verified against HEPData "
                         "(verified = false) -- run scripts/m02_fetch_hepdata.py first")
    table = prov_json.parent / prov["file"]
    if not table.exists():
        raise FileNotFoundError(f"{table} named in {prov_json.name} is not on disk "
                                "(data payloads are not committed; download it)")
    return prov, load_hepdata(table, int(prov.get("dependent_index", 0)))


# --------------------------------------------------------------------- MC
def read_mc(path: pathlib.Path):
    with uproot.open(path) as f:
        t = f["events"]
        keys = set(t.keys())
        want = ["particle_px", "particle_py", "particle_pz", "particle_charge"]
        for extra in ("process_code", "weight", "n_ch_eta1"):
            if extra in keys:
                want.append(extra)
        ev = t.arrays(want)
        prov = json.loads(str(f["provenance"])) if "provenance" in f else {}
    f64 = lambda a: ak.values_astype(a, np.float64)          # noqa: E731
    px, py, pz = f64(ev.particle_px), f64(ev.particle_py), f64(ev.particle_pz)
    pt = np.hypot(px, py)
    p = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
    eta = np.arctanh((pz / p) * (1.0 - 1e-12))
    return {"pt": pt, "eta": eta, "charged": ev.particle_charge != 0,
            "process_code": (ak.to_numpy(ev.process_code) if "process_code" in keys else None),
            "weight": (ak.to_numpy(ev.weight) if "weight" in keys
                       else np.ones(len(ev), dtype=float)),
            "n_ch_eta1_stored": (ak.to_numpy(ev.n_ch_eta1) if "n_ch_eta1" in keys else None),
            "provenance": prov}


def select_events(mc, event_class: str):
    """Boolean per-event mask for the experimental event class, derived from
    the particles (INEL>0) and the process code (NSD)."""
    n_ch_eta1 = ak.to_numpy(ak.sum(mc["charged"] & (abs(mc["eta"]) < 1.0), axis=1))
    if mc["n_ch_eta1_stored"] is not None and not np.array_equal(n_ch_eta1, mc["n_ch_eta1_stored"]):
        print("warning: n_ch_eta1 recomputed from particles differs from the stored "
              "branch -- the generator script and this analysis disagree on "
              "'charged' or 'eta'; using the recomputed value")
    cls = event_class.replace(" ", "").upper()
    if cls == "INEL":
        return np.ones(len(n_ch_eta1), bool)
    if cls in ("INEL>0", "INELGT0"):
        return n_ch_eta1 > 0
    if cls == "NSD":
        if mc["process_code"] is None:
            raise ValueError("NSD needs the process_code branch (PYTHIA samples only)")
        return ~np.isin(mc["process_code"], list(SD_CODES))
    raise ValueError(f"unknown event class {event_class!r}; use INEL, INEL>0 or NSD")


def dndeta(mc, mask, edges, pt_min=0.0):
    """Weighted dN_ch/deta per selected event in the given eta bins."""
    sel = mc["charged"] & (mc["pt"] > pt_min)
    eta = mc["eta"][mask][sel[mask]]
    w_ev = mc["weight"][mask]
    w_part = ak.to_numpy(ak.flatten(ak.broadcast_arrays(w_ev, eta)[0]))
    eta = ak.to_numpy(ak.flatten(eta))
    sumw, _ = np.histogram(eta, bins=edges, weights=w_part)
    sumw2, _ = np.histogram(eta, bins=edges, weights=w_part ** 2)
    n_ev = w_ev.sum()
    width = np.diff(edges)
    return sumw / (n_ev * width), np.sqrt(sumw2) / (n_ev * width)


def central_density(mc, mask, pt_min=0.0, eta_half=0.5):
    """<dN/deta> over |eta| < eta_half, binning-free, with the error of the
    per-event mean."""
    sel = mc["charged"] & (mc["pt"] > pt_min) & (abs(mc["eta"]) < eta_half)
    n = ak.to_numpy(ak.sum(sel[mask], axis=1)).astype(float)
    w = mc["weight"][mask]
    mean = np.sum(w * n) / np.sum(w)
    err = np.sqrt(np.sum(w * (n - mean) ** 2) / np.sum(w) / (len(n) - 1))
    return mean / (2 * eta_half), err / (2 * eta_half)


# --------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(REPO / "config" / "analysis_m02_pythia.yaml"))
    ap.add_argument("--mc", default=None, help="PYTHIA ROOT file (default: paths.generated_file)")
    ap.add_argument("--data", default=None,
                    help="provenance JSON of the measurement (default: comparison.data)")
    ap.add_argument("--event-class", default=None, help="INEL, INEL>0 or NSD")
    ap.add_argument("--no-data", action="store_true", help="MC only")
    ap.add_argument("--figdir", default=None)
    ap.add_argument("--tabdir", default=None)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    cmp_cfg = cfg.get("comparison", {})
    event_class = args.event_class or cfg["analysis"].get("event_selection", "INEL>0")
    pt_min = float(cmp_cfg.get("dndeta_pt_min", 0.0))
    mc_path = pathlib.Path(args.mc) if args.mc else REPO / cfg["paths"]["generated_file"]
    figdir = pathlib.Path(args.figdir) if args.figdir else REPO / cfg["paths"]["figure_dir"]
    tabdir = pathlib.Path(args.tabdir) if args.tabdir else REPO / cfg["paths"]["table_dir"]
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)
    sqrt_s = cfg["generator"]["sqrt_s_TeV"]

    # ---- data (optional) ----
    data, dprov = None, None
    if not args.no_data:
        djson = pathlib.Path(args.data) if args.data else (
            REPO / cmp_cfg["data"] if "data" in cmp_cfg else None)
        if djson is None or not djson.exists():
            print(f"[data] no provenance JSON ({djson}); running MC only")
        else:
            try:
                dprov, data = load_dataset(djson)
                print(f"[data] {dprov['experiment']} {dprov['observable']} "
                      f"{dprov['event_class']}, {len(data['y'])} bins, from {djson.name}")
                if dprov["event_class"].replace(" ", "").upper() != event_class.replace(" ", "").upper():
                    print(f"warning: data event class {dprov['event_class']!r} != "
                          f"MC selection {event_class!r}")
                if abs(float(dprov["sqrt_s_TeV"]) - sqrt_s) > 1e-6:
                    print(f"warning: data sqrt(s) {dprov['sqrt_s_TeV']} TeV != MC {sqrt_s} TeV")
            except (ValueError, FileNotFoundError) as exc:
                print(f"[data] {exc}\n[data] running MC only")

    # ---- MC ----
    mc = read_mc(mc_path)
    n_all = len(mc["weight"])
    mask = select_events(mc, event_class)
    n_sel = int(mask.sum())
    print(f"[mc]   {mc_path.name}: {n_all} events, {n_sel} pass {event_class} "
          f"({100.0 * n_sel / max(n_all, 1):.1f} %)")
    if mc["process_code"] is not None:
        codes, counts = np.unique(mc["process_code"], return_counts=True)
        print("[mc]   process mix: " + ", ".join(
            f"{PROCESS.get(int(c), str(c))} {100.0 * n / n_all:.1f} %" for c, n in zip(codes, counts)))

    edges = (np.concatenate([data["xlo"], [data["xhi"][-1]]]) if data is not None
             else np.linspace(-2.0, 2.0, 21))
    if data is not None and not np.allclose(edges[1:-1], data["xhi"][:-1]):
        print("warning: data bins are not contiguous; MC is histogrammed with the "
              "data's low edges + last high edge")
    y_mc, e_mc = dndeta(mc, mask, edges, pt_min)
    c_mc, ce_mc = central_density(mc, mask, pt_min)
    print(f"[mc]   dN_ch/deta (|eta|<0.5, pT > {pt_min:g}) = {c_mc:.3f} +- {ce_mc:.3f}")

    # ---- verdict ----
    verdict = None
    if data is not None:
        central = np.abs(data["x"]) < 0.5
        if central.any():
            w = (data["xhi"] - data["xlo"])[central]
            c_data = np.sum(data["y"][central] * w) / np.sum(w)
            c_data_err = np.sqrt(np.sum(((data["stat"][central] * w) ** 2))) / np.sum(w)
            c_data_sys = np.sum(data["sys"][central] * w) / np.sum(w)   # fully correlated
        else:
            i = np.argmin(np.abs(data["x"]))
            c_data, c_data_err, c_data_sys = data["y"][i], data["stat"][i], data["sys"][i]
        ratio = c_mc / c_data
        verdict = {"mc": c_mc, "mc_err": ce_mc, "data": c_data, "data_stat": c_data_err,
                   "data_sys": c_data_sys, "ratio": ratio,
                   "tolerance": VERDICT_TOLERANCE,
                   "within_tolerance": bool(abs(ratio - 1.0) <= VERDICT_TOLERANCE)}
        print(f"[cmp]  data {c_data:.3f} +- {c_data_err:.3f} (stat) +- {c_data_sys:.3f} (sys); "
              f"MC/data = {ratio:.3f}  -> "
              + ("WITHIN" if verdict["within_tolerance"] else "OUTSIDE")
              + f" +-{100 * VERDICT_TOLERANCE:.0f} %")

    # ---- figure ----
    apply_style()
    ctr = 0.5 * (edges[1:] + edges[:-1])
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(5.6, 5.4), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.05})
    ax.errorbar(ctr, y_mc, yerr=e_mc, fmt="s", ms=4, color=OKABE_ITO[1],
                label=f"PYTHIA 8 ({event_class})")
    if data is not None:
        ax.bar(data["x"], 2 * data["sys"], bottom=data["y"] - data["sys"],
               width=data["xhi"] - data["xlo"], color=OKABE_ITO[0], alpha=0.25,
               linewidth=0, label=f"{dprov['experiment']} syst.")
        ax.errorbar(data["x"], data["y"], yerr=data["stat"], fmt="o", ms=4,
                    color=OKABE_ITO[0], label=f"{dprov['experiment']} ({dprov['event_class']})")
        tot = np.hypot(data["stat"], data["sys"])
        with np.errstate(divide="ignore", invalid="ignore"):
            r = y_mc / data["y"]
            r_err_mc = e_mc / data["y"]
            r_band = tot / data["y"]
        axr.bar(data["x"], 2 * r_band, bottom=1 - r_band, width=data["xhi"] - data["xlo"],
                color=OKABE_ITO[0], alpha=0.25, linewidth=0)
        axr.errorbar(ctr, r, yerr=r_err_mc, fmt="s", ms=4, color=OKABE_ITO[1])
        axr.axhspan(1 - VERDICT_TOLERANCE, 1 + VERDICT_TOLERANCE, color="0.85", zorder=0)
    axr.axhline(1.0, color="0.3", lw=1)
    ax.set_ylabel(r"$\mathrm{d}N_{\rm ch}/\mathrm{d}\eta$")
    ax.set_ylim(0, max(y_mc.max(), data["y"].max() if data is not None else 0) * 1.6)
    axr.set_ylabel("MC / data")
    axr.set_ylim(0.6, 1.4)
    axr.set_xlabel(r"$\eta$")
    stamp(ax, [f"pp, $\\sqrt{{s}}$ = {sqrt_s:.0f} TeV",
               "charged primaries" + (f", $p_{{\\rm T}} > {pt_min:g}$ GeV/$c$" if pt_min > 0 else ""),
               (rf"MC/data$|_{{|\eta|<0.5}}$ = {verdict['ratio']:.3f}" if verdict else "MC only")],
          loc="upper left")
    ax.legend(loc="lower center", ncol=2, fontsize=9)
    save_figure(fig, figdir / "m02_dndeta_vs_data")
    plt.close(fig)

    # ---- summary ----
    summary = {
        "mc_file": rel(mc_path), "event_class": event_class, "pt_min": pt_min,
        "n_events_generated": int(n_all), "n_events_selected": n_sel,
        "process_mix": ({PROCESS.get(int(c), str(c)): int(n) for c, n in
                         zip(*np.unique(mc["process_code"], return_counts=True))}
                        if mc["process_code"] is not None else None),
        "dNdeta_central_mc": c_mc, "dNdeta_central_mc_err": ce_mc,
        "dndeta_mc": hist_record(edges, y_mc, e_mc),
        "data": ({"provenance": dprov, "x": data["x"].tolist(), "y": data["y"].tolist(),
                  "stat": data["stat"].tolist(), "sys": data["sys"].tolist()}
                 if data is not None else None),
        "verdict": verdict,
        "mc_provenance": mc["provenance"],
    }
    out = tabdir / "m02_compare_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"[out]  {rel(figdir / 'm02_dndeta_vs_data.pdf')}, {rel(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
