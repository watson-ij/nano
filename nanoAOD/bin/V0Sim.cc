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
#include "SimDataFormats/Vertex/interface/SimVertex.h"

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
  fwlite::Handle<vector<SimVertex>> hv;

  auto fout = new TFile("v0sim.root", "RECREATE");
  auto tout = new TTree("v0sim", "v0sim");
  int nKs = 0;
  tout->Branch("nKs", &nKs, "nKs/I");
  float pvx, pvy, pvz;
  tout->Branch("pvx", &pvx, "pvx/F");
  tout->Branch("pvy", &pvy, "pvy/F");
  tout->Branch("pvz", &pvz, "pvz/F");

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
  make_branchf(vx);  make_branchf(vy);  make_branchf(vz);
  make_branchf(vtx);  make_branchf(vty);  make_branchf(vtz);
  // Trackersurfaceposition from Pions give position on final layer it seems, not interesting for us
  // make_branchf(vx1);  make_branchf(vy1);  make_branchf(vz1);
  // make_branchf(vx2);  make_branchf(vy2);  make_branchf(vz2);
  make_branchf(dtxy);  make_branchf(dtxyz);
  make_branchf(tau); make_branchf(bg);
  
  for (ev.toBegin(); !ev.atEnd(); ++ev) {
    for (auto a : vf) a->clear();
    for (auto b : vb) b->clear();
    for (auto b : vi) b->clear();

    h.getByLabel(ev, "g4SimHits");
    hv.getByLabel(ev, "g4SimHits");

    // It appears the Primary Vertex is always the first SimVertex in the list (makes sense)
    pvx = hv.product()->at(0).position().x();
    pvy = hv.product()->at(0).position().y();
    pvz = hv.product()->at(0).position().z();

    auto v = h.product();
    vector<const SimTrack*> ks;
    for (auto& trk : *h) {
      if (trk.type() != pdgId_KS) continue;
      ks.push_back(&trk);
      auto& ksp = trk.momentum();
      auto& ksv = trk.trackerSurfacePosition();
      bool ismatch = false;
      bool ispi0 = false;
      // Match pions
      for (unsigned int i = 0; i < (v->size()-1); ++i) {
	auto &pi1 = v->at(i);
	auto &pi2 = v->at(i+1);

	const SimTrack *p1 = 0, *p2 = 0;
	if (abs(pi1.type()) == pdgId_pip && abs(pi2.type()) == pdgId_pip) {
	  if (pi1.type() == pi2.type()) continue;
	  auto comb = pi1.momentum() + pi2.momentum();
	  if (fabs(comb.pt() - ksp.pt()) > feps) continue;
	  if (fabs(comb.eta() - ksp.eta()) > feps) continue;
	  if (fabs(comb.phi() - ksp.phi()) > feps) continue;
	  if (fabs(comb.mass() - ksp.mass()) > feps) continue;
	  ismatch = true;
	  p1 = &pi1;
	  p2 = &pi2;
	} else if (abs(pi1.type()) == pdgId_pi0 && abs(pi2.type()) == pdgId_pi0) {
	  auto comb = pi1.momentum() + pi2.momentum();
	  if (fabs(comb.pt() - ksp.pt()) > feps) continue;
	  if (fabs(comb.eta() - ksp.eta()) > feps) continue;
	  if (fabs(comb.phi() - ksp.phi()) > feps) continue;
	  if (fabs(comb.mass() - ksp.mass()) > feps) continue;
	  ispi0 = true;
	  p1 = &pi1;
	  p2 = &pi2;
	}
	if (p1) {

	  pt1.push_back(p1->momentum().pt());
	  eta1.push_back(p1->momentum().eta());
	  phi1.push_back(p1->momentum().phi());
	  pdg1.push_back(p2->type());

	  pt2.push_back(p2->momentum().pt());
	  eta2.push_back(p2->momentum().eta());
	  phi2.push_back(p2->momentum().phi());
	  pdg2.push_back(p2->type());

	  // vx1.push_back(p1->trackerSurfacePosition().x());vy1.push_back(p1->trackerSurfacePosition().y());vz1.push_back(p1->trackerSurfacePosition().z());
	  // vx2.push_back(p2->trackerSurfacePosition().x());vy2.push_back(p2->trackerSurfacePosition().y());vz2.push_back(p2->trackerSurfacePosition().z());

	  // As far as I can tell, this is the position of the Kshort decay vertex
	  int vi = p1->vertIndex();
	  if (vi >= 0) {
	    auto vt = hv.product()->at(vi);
	    vtx.push_back(vt.position().x());
	    vty.push_back(vt.position().y());
	    vtz.push_back(vt.position().z());

	    dtxy.push_back(sqrt(pow(vt.position().x() - pvx, 2) + pow(vt.position().y() - pvy, 2)));
	    float dl = sqrt(pow(vt.position().x() - pvx, 2) + pow(vt.position().y() - pvy, 2) + pow(vt.position().z() - pvz, 2));
	    dtxyz.push_back(dl);
	    float bgam = ksp.Beta()*ksp.Gamma();
	    bg.push_back(bgam);
	    tau.push_back(dl * (1.0 / 300.0) / bgam); // tau (in ps) = l / (beta*gamma*c), l in mm (?)
	  }
	  
	  break;
	}
      }

      if (!ispi0 && !ismatch) {
	  pt1.push_back(-99);	  eta1.push_back(-99);	  phi1.push_back(-99);	  pdg1.push_back(-1);
	  pt2.push_back(-99);	  eta2.push_back(-99);	  phi2.push_back(-99);	  pdg2.push_back(-1);
	  // vx1.push_back(-999); vy1.push_back(-999); vz1.push_back(-999);
	  // vx2.push_back(-999); vy2.push_back(-999); vz2.push_back(-999);
	  vtx.push_back(-999); vty.push_back(-999); vtz.push_back(-999);
	  dtxy.push_back(-999); dtxyz.push_back(-999); tau.push_back(-99); bg.push_back(-99);
      }
      
      match.push_back(ismatch);
      pi0match.push_back(ispi0);
      pt.push_back(ksp.pt());
      phi.push_back(ksp.phi());
      eta.push_back(ksp.eta());
      // This seems to be production position, should equal PV?
      vx.push_back(ksv.x());
      vy.push_back(ksv.y());
      vz.push_back(ksv.z());
    }

    nKs = ks.size();
    tout->Fill();
  }

  tout->Write();
  fout->Close();
}
