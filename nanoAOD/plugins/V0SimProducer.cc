/*
  From cmssw RecoVertex/V0Producer/src/V0Fitter.cc
*/

#include<memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

using namespace edm;
using namespace std;
using namespace reco;

class V0SimProducer : public edm::stream::EDProducer<> {
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;
    
public:
  explicit V0SimProducer(const edm::ParameterSet & iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  void produce( edm::Event&, const edm::EventSetup& ) override;

  edm::InputTag SimTkLabel_;
  edm::InputTag SimVtxLabel_;

  const float gPionMass = 0.1396;
  const float gKaonMass = 0.4937;
  const float gJpsiMass = 3.096;
  const float gD0Mass = 1.865;
  const float gDstarMass = 2.010;
  const float gProtonMass = 0.938272;
  const int pdgId_Kp = 321;
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
};

V0SimProducer::V0SimProducer(const edm::ParameterSet & iConfig) :
  SimTkLabel_(iConfig.getParameter<edm::InputTag>("moduleLabelTk")),
  SimVtxLabel_(iConfig.getParameter<edm::InputTag>("moduleLabelVtx"))
{
  produces<nanoaod::FlatTable>("V0Sim");
}

void
V0SimProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::vector<float> v0pt, v0eta, v0phi,
    v0pt1, v0pt2, v0eta1, v0eta2, v0phi1, v0phi2;
  std::vector<int> v0pdgId1, v0pdgId2;
  std::vector<bool> foundD;
  float eps = 1e-4;
    Handle<edm::View<SimTrack>> SimTk;
    //Handle<SimVertexContainer> SimVtx;
    // iEvent.getByLabel(SimTkLabel_, SimTk);
    // iEvent.getByLabel(SimVtxLabel_, SimVtx);
    // iEvent.getByLabel("g4SimHits", SimTk);
    iEvent.get(edm::ProductID(1, 47), SimTk);
    //iEvent.getByLabel(SimVtxLabel_, SimVtx);
    
    for (unsigned int i = 0; i < SimTk->size(); ++i) {
      auto &tk = SimTk->at(i);
      if (tk.type() == 310) {
	v0pt.push_back(tk.momentum().pt());
	v0eta.push_back(tk.momentum().eta());
	v0phi.push_back(tk.momentum().phi());
	
	// find the matching pions
	float pt1 = -1, pt2 = -1, eta1 = -1, eta2 = -1, phi1 = -1, phi2 = -1;
	int pdgId1 = -1, pdgId2 = -1;
	bool found = false;
	for (unsigned int i = 0; i < SimTk->size(); ++i) {
	  auto &tk1 = SimTk->at(i);
	  if (abs(tk1.type()) != 211) { continue; }
	  auto &tk2 = SimTk->at(i+1);
	  if (abs(tk2.type()) != 211) { continue; }
	  if (tk1.type() == tk2.type()) { continue; }
	  auto comb = tk1.momentum() + tk2.momentum();
	  if (fabs(comb.pt() - tk.momentum().pt()) > eps) { continue; }
	  if (fabs(comb.eta() - tk.momentum().eta()) > eps) { continue; }
	  if (fabs(comb.phi() - tk.momentum().phi()) > eps) { continue; }
	  if (fabs(comb.mass() - tk.momentum().mass()) > eps) { continue; }
	  found = true;
	  pt1 = tk1.momentum().pt(); pt2 = tk2.momentum().pt();
	  eta1 = tk1.momentum().eta(); eta2 = tk2.momentum().eta();
	  phi1 = tk1.momentum().phi(); phi2 = tk2.momentum().phi();
	  pdgId1 = tk1.type(); pdgId2 = tk2.type();
	}
	foundD.push_back(found);
	v0pt1.push_back(pt1); v0pt2.push_back(pt2);
	v0eta1.push_back(eta1); v0eta2.push_back(eta2);
	v0phi1.push_back(phi1); v0phi2.push_back(phi2);
	v0pdgId1.push_back(pdgId1); v0pdgId2.push_back(pdgId2);
      }
    }
  
  auto v0simTable = make_unique<nanoaod::FlatTable>(0,"V0Sim",false);
  v0simTable->addColumn<float>("pt", v0pt,"V0 pt",nanoaod::FlatTable::FloatColumn);  
  v0simTable->addColumn<float>("eta", v0eta,"V0 eta",nanoaod::FlatTable::FloatColumn);  
  v0simTable->addColumn<float>("phi", v0phi,"V0 phi",nanoaod::FlatTable::FloatColumn);  
  // v0simTable->addColumn<bool>("hasDaughts", foundD,"true if found daughters and filled daughter info",nanoaod::FlatTable::FloatColumn);  
  v0simTable->addColumn<float>("pt1", v0pt1,"V0 pt of first daughter (-1 if doesn't exist)",nanoaod::FlatTable::FloatColumn);
  v0simTable->addColumn<float>("eta1", v0eta1,"V0 eta of first daughter (-1 if doesn't exist)",nanoaod::FlatTable::FloatColumn);
  v0simTable->addColumn<float>("phi1", v0phi1,"V0 phi of first daughter (-1 if doesn't exist)",nanoaod::FlatTable::FloatColumn);
  v0simTable->addColumn<float>("pdgId1", v0pdgId1,"V0 PDG of first daughter (-1 if doesn't exist)",nanoaod::FlatTable::FloatColumn);
  v0simTable->addColumn<float>("pt2", v0pt2,"V0 pt of second daughter",nanoaod::FlatTable::FloatColumn);  
  v0simTable->addColumn<float>("eta2", v0eta2,"V0 eta of second daughter",nanoaod::FlatTable::FloatColumn);  
  v0simTable->addColumn<float>("phi2", v0phi2,"V0 phi of second daughter",nanoaod::FlatTable::FloatColumn);  
  v0simTable->addColumn<float>("pdgId2", v0pdgId2,"V0 PDG of second daughter",nanoaod::FlatTable::FloatColumn);

  iEvent.put(move(v0simTable),"V0Sim");
}

void
V0SimProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(V0SimProducer);
