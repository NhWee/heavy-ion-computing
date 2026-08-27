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

import pathlib
import sys

import ROOT

REPO = pathlib.Path(__file__).resolve().parent.parent

# C++ helpers compiled once by cling and then usable inside Define() strings.
ROOT.gInterpreter.Declare(r"""
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

    infile = REPO / "data" / "generated" / "m01_toy_events.root"
    if not infile.exists():
        print(f"missing {infile}; run scripts/m01_make_toy_events.py first")
        return 1

    df = (ROOT.RDataFrame("events", str(infile))
          .Define("pt", "calcPt(particle_px, particle_py)")
          .Define("eta", "calcEta(particle_px, particle_py, particle_pz)")
          .Define("apdg", "abs(particle_pdg)")
          .Define("isCh", "apdg == 211 || apdg == 321 || apdg == 2212")
          .Define("selTrk", "isCh && pt > 0.15f && abs(eta) < 0.8f")
          .Define("selMu", "apdg == 13 && pt > 20.f && abs(eta) < 2.4f")
          .Define("nch", "Sum(selTrk)")
          .Define("etaCh", "eta[isCh && pt > 0.15f]")
          .Define("ptSel", "pt[selTrk]")
          .Define("mmumu", "dimuonMass(particle_px, particle_py, particle_pz,"
                           " particle_E, selMu)"))

    n_ev = df.Count().GetValue()

    h_mult = df.Histo1D(("hMult", ";N_{ch};P(N_{ch})", 60, 0, 240), "nch")
    h_eta = df.Histo1D(("hEta", ";#eta;dN_{ch}/d#eta", 60, -6, 6), "etaCh")
    h_pt = df.Histo1D(("hPt", ";p_{T} (GeV/c);counts", 40, 0.15, 50.0), "ptSel")
    h_mass = df.Filter("mmumu > 0").Histo1D(
        ("hMass", ";m_{#mu#mu} (GeV/c^{2});counts", 60, 60, 120), "mmumu")

    for h in (h_mult, h_eta, h_pt):
        h.GetValue().Scale(1.0 / n_ev, "width")

    print(f"[read] {infile.name}: {n_ev} events")
    print(f"[res]  <Nch> = {h_mult.GetMean():.2f}")
    print(f"[res]  dN/deta at eta=0 = "
          f"{h_eta.GetBinContent(h_eta.FindBin(0.0)):.2f}")
    print(f"[res]  dimuon candidates = {int(h_mass.GetEntries())}, "
          f"mean m = {h_mass.GetMean():.2f} GeV/c2")

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
