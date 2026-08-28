// ---------------------------------------------------------------------------
// Module 1, ROOT-native arm:  the same analysis as m01_analysis_uproot.py,
// written the way it would be written inside a collaboration framework.
//
//   root -l -b -q 'analysis/m01_analysis_root.C("data/generated/m01_toy_events.root")'
//
// NOTE: the default arguments below MUST match config/analysis.yaml.  A C++
// macro cannot read the YAML, so this is the one place in the repository where a
// number is duplicated -- if you change the config, change these too.
//
// STATUS: written against the ROOT 6.32+ API but NOT executed here -- the
// container this repository was bootstrapped in has no ROOT installation
// (conda-forge is unreachable from it).  Run it in the WSL2 `hic` environment
// and record the outcome in docs/troubleshooting.md.
//
// What it demonstrates
//   TFile      -- the container; a ROOT file is a small hierarchical filesystem
//   TTree      -- the event table, with variable-length ("jagged") branches
//   TH1F/TH2F  -- histograms that carry their own binning and errors (Sumw2)
//   TLorentzVector / ROOT::Math::PtEtaPhiMVector -- four-momentum algebra
//   TF1 + Fit  -- the ROOT fitting interface
// ---------------------------------------------------------------------------
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TVirtualPad.h>
#include <TStyle.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <iostream>

// Tsallis invariant yield, used as a TF1
Double_t TsallisInvariant(Double_t *x, Double_t *par) {
   const Double_t pt = x[0];
   const Double_t C = par[0], T = par[1], n = par[2], m0 = par[3];
   const Double_t mt = TMath::Sqrt(pt * pt + m0 * m0);
   return C * TMath::Power(1.0 + (mt - m0) / (n * T), -n);
}

void m01_analysis_root(const char *infile = "data/generated/m01_toy_events.root",
                       const char *outfile = "results/tables/m01_root_histograms.root",
                       Double_t etaMax = 0.8, Double_t ptMin = 0.15,
                       Double_t muPtMin = 20.0, Double_t muEtaMax = 2.4) {
   gStyle->SetOptStat(0);

   TFile *fin = TFile::Open(infile, "READ");
   if (!fin || fin->IsZombie()) { std::cerr << "cannot open " << infile << "\n"; return; }
   TTree *tree = dynamic_cast<TTree *>(fin->Get("events"));
   if (!tree) { std::cerr << "no TTree named 'events'\n"; return; }

   // --- branch buffers -----------------------------------------------------
   const Int_t kMax = 20000;
   Int_t nparticle = 0;
   static Float_t px[kMax], py[kMax], pz[kMax], E[kMax];
   static Int_t pdg[kMax];
   tree->SetBranchAddress("nparticle", &nparticle);
   tree->SetBranchAddress("particle_px", px);
   tree->SetBranchAddress("particle_py", py);
   tree->SetBranchAddress("particle_pz", pz);
   tree->SetBranchAddress("particle_E", E);
   tree->SetBranchAddress("particle_pdg", pdg);

   // --- histograms: log bins in pT, linear elsewhere ------------------------
   const Int_t nPt = 40;
   Double_t ptEdges[nPt + 1];
   for (Int_t i = 0; i <= nPt; ++i)
      ptEdges[i] = 0.15 * TMath::Power(50.0 / 0.15, Double_t(i) / nPt);

   TH1F *hPt   = new TH1F("hPt", ";p_{T} (GeV/c);1/(2#pi p_{T} N_{ev}) d^{2}N/dp_{T}d#eta",
                          nPt, ptEdges);
   TH1F *hEta  = new TH1F("hEta", ";#eta;dN_{ch}/d#eta", 60, -6.0, 6.0);
   TH1F *hMult = new TH1F("hMult", ";N_{ch};P(N_{ch})", 60, 0.0, 240.0);
   TH1F *hMass = new TH1F("hMass", ";m_{#mu^{+}#mu^{-}} (GeV/c^{2});counts",
                          60, 60.0, 120.0);
   for (auto *h : {hPt, hEta, hMult, hMass}) h->Sumw2();   // track errors properly

   const Long64_t nev = tree->GetEntries();
   std::cout << "[read] " << infile << ": " << nev << " events" << std::endl;

   for (Long64_t iev = 0; iev < nev; ++iev) {
      tree->GetEntry(iev);
      Int_t nch = 0;
      Int_t nmu = 0;
      ROOT::Math::PxPyPzEVector muSum(0, 0, 0, 0);

      for (Int_t i = 0; i < nparticle; ++i) {
         const Int_t ap = TMath::Abs(pdg[i]);
         const Double_t pt = TMath::Sqrt(px[i] * px[i] + py[i] * py[i]);
         const Double_t p  = TMath::Sqrt(pt * pt + pz[i] * pz[i]);
         if (p <= 0) continue;
         const Double_t eta = 0.5 * TMath::Log((p + pz[i]) / (p - pz[i]));

         const bool chargedHadron = (ap == 211 || ap == 321 || ap == 2212);
         if (chargedHadron && pt > ptMin) {
            hEta->Fill(eta);
            if (TMath::Abs(eta) < etaMax) {
               ++nch;
               // invariant yield: weight 1/(2 pi pT Delta eta)
               hPt->Fill(pt, 1.0 / (2.0 * TMath::Pi() * pt * 2.0 * etaMax));
            }
         }
         if (ap == 13 && pt > muPtMin && TMath::Abs(eta) < muEtaMax) {
            ++nmu;
            muSum += ROOT::Math::PxPyPzEVector(px[i], py[i], pz[i], E[i]);
         }
      }
      hMult->Fill(nch);
      if (nmu == 2) hMass->Fill(muSum.M());
   }

   // --- normalise per event and per bin width ------------------------------
   for (auto *h : {hPt, hEta, hMult}) {
      h->Scale(1.0 / Double_t(nev), "width");
   }

   // --- fit the pT spectrum -------------------------------------------------
   TF1 *fTsallis = new TF1("fTsallis", TsallisInvariant, 0.15, 50.0, 4);
   fTsallis->SetParameters(hPt->GetBinContent(1) * 2.0, 0.13, 7.0, 0.13957039);
   fTsallis->FixParameter(3, 0.13957039);              // pion mass, fixed
   fTsallis->SetParNames("C", "T", "n", "m0");
   hPt->Fit(fTsallis, "R Q S");
   std::cout << "[fit] T = " << fTsallis->GetParameter(1) * 1000.0 << " +- "
             << fTsallis->GetParError(1) * 1000.0 << " MeV, n = "
             << fTsallis->GetParameter(2) << " +- " << fTsallis->GetParError(2)
             << ", chi2/ndf = " << fTsallis->GetChisquare() / fTsallis->GetNDF()
             << std::endl;

   // --- draw & save ---------------------------------------------------------
   TCanvas *c = new TCanvas("c", "Module 1", 1200, 900);
   c->Divide(2, 2);
   c->cd(1); gPad->SetLogy(); hMult->Draw("E");
   c->cd(2); hEta->Draw("E");
   c->cd(3); gPad->SetLogx(); gPad->SetLogy(); hPt->Draw("E");
   c->cd(4); hMass->Draw("E");
   c->SaveAs("results/figures/m01_root_overview.pdf");

   TFile *fout = TFile::Open(outfile, "RECREATE");
   hPt->Write(); hEta->Write(); hMult->Write(); hMass->Write(); fTsallis->Write();
   fout->Close();
   fin->Close();
   std::cout << "[out] histograms -> " << outfile << std::endl;
}
