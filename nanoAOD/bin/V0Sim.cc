#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

#include<iostream>
#include<memory>
#include<vector>

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

#include "boost/any.hpp"

using namespace std;
  const float gPionMass = 0.1396;
  const float gKaonMass = 0.4937;
  const float gJpsiMass = 3.096;
  const float gD0Mass = 1.865;
  const float gDstarMass = 2.010;
  const float gProtonMass = 0.938272;
  const int pdgId_Kp = 321;
  const int pdgId_pip = 211;
  const int pdgId_pi0 = 111;
  const int pdgId_Jpsi = 443;
  const int pdgId_D0 = 421;
  const int pdgId_Dstar = 413;
  const int pdgId_KS = 310;
  const int pdgId_Lambda = 3122;
  const int pdgId_p = 2122;
  const float jpsiMin_ = 2.5;
  const float jpsiMax_ = 3.4;
  const float D0Min_   = 1.7;
  const float D0Max_   = 2.0;
  const float DstarDiffMin_   = 0.14;
  const float DstarDiffMax_   = 0.16;
  const float KSMin_   = 0.4;
  const float KSMax_   = 0.6;
  const float LambdaMin_   = 0.9;
  const float LambdaMax_   = 1.16;

const double feps = 1e-4;
int main(int argc, char *argv[]) {
  if (argc != 2) { cout << "Please give me an AOD to analyse" << endl; exit(111); }
  FWLiteEnabler::enable();
  cout << "Loading file " << argv[1] << endl;
  auto f = new TFile(argv[1]);
  fwlite::Event ev(f);
  fwlite::Handle<vector<SimTrack>> h;

  auto fout = new TFile("v0sim.root", "RECREATE");
  auto tout = new TTree("v0sim", "v0sim");
  int nKs = 0;
  tout->Branch("nKs", &nKs, "nKs/I");

  vector<vector<float>*> vf;
  vector<vector<bool>*> vb;
  vector<vector<int>*> vi;
#define make_branch(typ, name) std::vector<typ> name; tout->Branch(#name, &name);
#define make_branchf(name) make_branch(float, name); vf.push_back(&name);
#define make_branchb(name) make_branch(bool, name); vb.push_back(&name);
#define make_branchi(name) make_branch(int, name); vi.push_back(&name);
  make_branchb(match);
  make_branchb(pi0match);
  make_branchf(pt);  make_branchf(eta);  make_branchf(phi);
  make_branchf(pt1);  make_branchf(eta1);  make_branchf(phi1);  make_branchi(pdg1);
  make_branchf(pt2);  make_branchf(eta2);  make_branchf(phi2);  make_branchi(pdg2);
  
  for (ev.toBegin(); !ev.atEnd(); ++ev) {
    for (auto a : vf) a->clear();
    for (auto b : vb) b->clear();
    for (auto b : vi) b->clear();
    
    h.getByLabel(ev, "g4SimHits");
    auto v = h.product();
    vector<const SimTrack*> ks;
    for (auto& trk : *h) {
      if (trk.type() != pdgId_KS) continue;
      ks.push_back(&trk);
      auto& ksp = trk.momentum();
      bool ismatch = false;
      bool ispi0 = false;
      // Match pions
      for (unsigned int i = 0; i < (v->size()-1); ++i) {
	auto &pi1 = v->at(i);
	auto &pi2 = v->at(i+1);
	if (abs(pi1.type()) == pdgId_pip && abs(pi2.type()) == pdgId_pip) {
	  if (pi1.type() == pi2.type()) continue;
	  auto comb = pi1.momentum() + pi2.momentum();
	  if (fabs(comb.pt() - ksp.pt()) > feps) continue;
	  if (fabs(comb.eta() - ksp.eta()) > feps) continue;
	  if (fabs(comb.phi() - ksp.phi()) > feps) continue;
	  if (fabs(comb.mass() - ksp.mass()) > feps) continue;
	  ismatch = true;

	  pt1.push_back(pi1.momentum().pt());
	  eta1.push_back(pi1.momentum().eta());
	  phi1.push_back(pi1.momentum().phi());
	  pdg1.push_back(pi2.type());

	  pt2.push_back(pi2.momentum().pt());
	  eta2.push_back(pi2.momentum().eta());
	  phi2.push_back(pi2.momentum().phi());
	  pdg2.push_back(pi2.type());
	  
	  break;
	}
	else if (abs(pi1.type()) == pdgId_pi0 && abs(pi2.type()) == pdgId_pi0) {
	  auto comb = pi1.momentum() + pi2.momentum();
	  if (fabs(comb.pt() - ksp.pt()) > feps) continue;
	  if (fabs(comb.eta() - ksp.eta()) > feps) continue;
	  if (fabs(comb.phi() - ksp.phi()) > feps) continue;
	  if (fabs(comb.mass() - ksp.mass()) > feps) continue;
	  ispi0 = true;

	  pt1.push_back(pi1.momentum().pt());
	  eta1.push_back(pi1.momentum().eta());
	  phi1.push_back(pi1.momentum().phi());
	  pdg1.push_back(pi2.type());

	  pt2.push_back(pi2.momentum().pt());
	  eta2.push_back(pi2.momentum().eta());
	  phi2.push_back(pi2.momentum().phi());
	  pdg2.push_back(pi2.type());

	  break;
	}
      }

      if (!ispi0 && !ismatch) {
	  pt1.push_back(-99);
	  eta1.push_back(-99);
	  phi1.push_back(-99);
	  pdg1.push_back(-1);
	  pt2.push_back(-99);
	  eta2.push_back(-99);
	  phi2.push_back(-99);
	  pdg2.push_back(-1);
      }
      
      match.push_back(ismatch);
      pi0match.push_back(ispi0);
      pt.push_back(ksp.pt());
      phi.push_back(ksp.phi());
      eta.push_back(ksp.eta());
    }

    nKs = ks.size();
    tout->Fill();
  }

  tout->Write();
  fout->Close();
}
