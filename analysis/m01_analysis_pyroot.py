#!/usr/bin/env python3
"""Module 1, PyROOT / RDataFrame arm.

    python analysis/m01_analysis_pyroot.py

STATUS: written against the ROOT 6.32 API but NOT executed in the bootstrap
container (no ROOT there).  Run it inside the WSL2 `hic` environment; record any
API drift in docs/troubleshooting.md.

Why RDataFrame
--------------
The event loop in m01_analysis_root.C is explicit: you write the `for` loop.
RDataFrame instead lets you *declare* the transformations (Define / Filter /
Histo1D) and hands the loop to ROOT, which can then run it multi-threaded and
lazily -- nothing is read from disk until you ask for a result.  This is the
interface modern ROOT analyses are written in, and it is the closest ROOT
analogue of the uproot+awkward style used in m01_analysis_uproot.py.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import ROOT
import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent

# C++ helpers compiled once by cling and then usable inside Define() strings.
ROOT.gInterpreter.Declare(r"""
#include "ROOT/RVec.hxx"
#include <cmath>

using RVecF = ROOT::VecOps::RVec<float>;
using RVecI = ROOT::VecOps::RVec<int>;

RVecF calcPt(const RVecF &px, const RVecF &py) { return sqrt(px*px + py*py); }

RVecF calcEta(const RVecF &px, const RVecF &py, const RVecF &pz) {
   auto p = sqrt(px*px + py*py + pz*pz);
   RVecF out(px.size());
   for (size_t i = 0; i < px.size(); ++i)
      out[i] = (p[i] > 0.f) ? 0.5f * std::log((p[i] + pz[i]) / (p[i] - pz[i])) : 0.f;
   return out;
}

// invariant mass of exactly-two selected muons (returns -1 if not exactly two)
float dimuonMass(const RVecF &px, const RVecF &py, const RVecF &pz,
                 const RVecF &E, const ROOT::VecOps::RVec<int> &sel) {
   float sx = 0, sy = 0, sz = 0, se = 0; int n = 0;
   for (size_t i = 0; i < sel.size(); ++i)
      if (sel[i]) { sx += px[i]; sy += py[i]; sz += pz[i]; se += E[i]; ++n; }
   if (n != 2) return -1.f;
   float m2 = se*se - (sx*sx + sy*sy + sz*sz);
   return m2 > 0 ? std::sqrt(m2) : 0.f;
}
""")


def main() -> int:
    ROOT.gROOT.SetBatch(True)
    ROOT.EnableImplicitMT()               # RDataFrame goes multi-threaded

    cfg = yaml.safe_load(open(REPO / "config" / "analysis.yaml"))
    acut = cfg["analysis"]["track_cuts"]
    mcut = cfg["analysis"]["muon_cuts"]
    bpt = cfg["analysis"]["binning"]["pt"]
    eta_max, pt_min = acut["eta_abs_max"], acut["pt_min"]
    mu_pt_min, mu_eta_max = mcut["pt_min"], mcut["eta_abs_max"]
    delta_eta = 2.0 * eta_max

    infile = REPO / cfg["paths"]["generated_file"]
    if not infile.exists():
        print(f"missing {infile}; run scripts/m01_make_toy_events.py first")
        return 1

    # identical log bin edges to the uproot arm -- bin edges are the single most
    # common source of disagreement between two implementations
    pt_edges = np.geomspace(bpt["lo"], bpt["hi"], bpt["n_bins"] + 1)

    df = (ROOT.RDataFrame("events", str(infile))
          .Define("pt", "calcPt(particle_px, particle_py)")
          .Define("eta", "calcEta(particle_px, particle_py, particle_pz)")
          .Define("apdg", "abs(particle_pdg)")
          .Define("isCh", "apdg == 211 || apdg == 321 || apdg == 2212")
          .Define("selTrk", f"isCh && pt > {pt_min}f && abs(eta) < {eta_max}f")
          .Define("selMu", f"apdg == 13 && pt > {mu_pt_min}f "
                           f"&& abs(eta) < {mu_eta_max}f")
          .Define("nch", "Sum(selTrk)")
          .Define("etaCh", f"eta[isCh && pt > {pt_min}f]")
          .Define("ptSel", "pt[selTrk]")
          .Define("nEta05", f"Sum(isCh && pt > {pt_min}f && abs(eta) < 0.5f)")
          # per-particle 1/(2 pi pT d_eta) Jacobian, same convention as the
          # uproot arm (NOT the bin centre)
          .Define("ptWgt", f"1.f / (2.f * (float)M_PI * ptSel * {delta_eta}f)")
          .Define("mmumu", "dimuonMass(particle_px, particle_py, particle_pz,"
                           " particle_E, selMu)"))

    n_ev = df.Count().GetValue()
    mean_nch = df.Mean("nch")                 # exact, not the binned TH1 mean
    mean_eta05 = df.Mean("nEta05")
    mean_mass = df.Filter("mmumu > 0").Mean("mmumu")

    h_mult = df.Histo1D(("hMult", ";N_{ch};P(N_{ch})", 60, 0, 240), "nch")
    h_eta = df.Histo1D(("hEta", ";#eta;dN_{ch}/d#eta", 60, -6, 6), "etaCh")
    h_pt = df.Histo1D(
        ("hPt", ";p_{T} (GeV/c);1/(2#pi p_{T}N_{ev}) d^{2}N/dp_{T}d#eta",
         len(pt_edges) - 1, pt_edges), "ptSel", "ptWgt")
    h_mass = df.Filter("mmumu > 0").Histo1D(
        ("hMass", ";m_{#mu#mu} (GeV/c^{2});counts", 60, 60, 120), "mmumu")

    for h in (h_mult, h_eta, h_pt):
        h.GetValue().Sumw2()
        h.GetValue().Scale(1.0 / n_ev, "width")

    summary = {
        "arm": "pyroot",
        "n_events": int(n_ev),
        "mean_Nch": float(mean_nch.GetValue()),
        "dNdeta_central": float(mean_eta05.GetValue() / 1.0),
        "n_dimuon_candidates": int(h_mass.GetEntries()),
        "mean_dimuon_mass_GeV": float(mean_mass.GetValue()),
        "root_version": ROOT.gROOT.GetVersion(),
    }
    print(f"[read] {infile.name}: {n_ev} events")
    print(f"[res]  <Nch> = {summary['mean_Nch']:.3f}")
    print(f"[res]  dN/deta (|eta|<0.5) = {summary['dNdeta_central']:.3f}")
    print(f"[res]  dimuon candidates = {summary['n_dimuon_candidates']}, "
          f"mean m = {summary['mean_dimuon_mass_GeV']:.3f} GeV/c2")
    (REPO / "results" / "tables" / "m01_pyroot_summary.json").write_text(
        json.dumps(summary, indent=2))

    out = REPO / "results" / "tables" / "m01_pyroot_histograms.root"
    out.parent.mkdir(parents=True, exist_ok=True)
    fout = ROOT.TFile(str(out), "RECREATE")
    for h in (h_mult, h_eta, h_pt, h_mass):
        h.GetValue().Write()
    fout.Close()
    print(f"[out]  {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
