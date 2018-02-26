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

using namespace edm;
using namespace std;
using namespace reco;

class MesonProducer : public edm::stream::EDProducer<> {
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;
  struct trackVars {
    float normalizedChi2;
    int nHits;
    float trkPt;
    float ipsigXY;
    float ipsigZ;
  };
    
public:
  explicit MesonProducer(const edm::ParameterSet & iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  void produce( edm::Event&, const edm::EventSetup& ) override;

  reco::VertexCompositeCandidate fit(vector<const pat::PackedCandidate*> cands,
				     reco::Vertex pv, int pdgId,
				     float &dca, float &angleXY, float &angleXYZ);
    
  SVector3 getDistanceVector(int dim, reco::VertexCompositeCandidate vertex,reco::Vertex pv);
  pair<float, float> getDistance(int dim, reco::VertexCompositeCandidate vertex,reco::Vertex pv);
  trackVars getTrackVars(vector<const pat::PackedCandidate*> cands, reco::Vertex pv);

  int findMCmatch(const pat::Jet & aPatJet, vector<const pat::PackedCandidate*> cands, int pdgid);
  
  edm::EDGetTokenT<edm::View<pat::Jet> > jetSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfSrc_;
  edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
  //edm::EDGetTokenT<edm::View<reco::GenParticle> > mcSrc_;
  edm::ESHandle<TransientTrackBuilder> trackBuilder_;

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
  //unsigned int maxNumPFCand_;
  bool doMatch_;
  bool applyCuts_;
  // cuts on initial track selection
  float tkChi2Cut_;
  int tkNHitsCut_;
  float tkPtCut_;
  float tkIPSigXYCut_;
  float tkIPSigZCut_;
  // cuts on the vertex
  float vtxChi2Cut_;
  float vtxDecaySigXYCut_;
  float vtxDecaySigXYZCut_;
  // miscellaneous cuts
  float tkDCACut_;
  float cosThetaXYCut_;
  float cosThetaXYZCut_;  
};

MesonProducer::MesonProducer(const edm::ParameterSet & iConfig) :
  jetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetLabel"))), 
  pfSrc_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfLabel"))),
  vertexLabel_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"))),
  //mcSrc_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("mcLabel"))),
  doMatch_(iConfig.getParameter<bool>("doMCMatching")),
  applyCuts_(iConfig.getParameter<bool>("applySoftLeptonCut"))  
{
  // cuts on initial track selection
  tkChi2Cut_ = iConfig.getParameter<double>("tkChi2Cut");
  tkNHitsCut_ = iConfig.getParameter<int>("tkNHitsCut");
  tkPtCut_ = iConfig.getParameter<double>("tkPtCut");
  tkIPSigXYCut_ = iConfig.getParameter<double>("tkIPSigXYCut");
  tkIPSigZCut_ = iConfig.getParameter<double>("tkIPSigZCut");   
  // cuts on vertex
  vtxChi2Cut_ = iConfig.getParameter<double>("vtxChi2Cut");
  vtxDecaySigXYZCut_ = iConfig.getParameter<double>("vtxDecaySigXYZCut");
  vtxDecaySigXYCut_ = iConfig.getParameter<double>("vtxDecaySigXYCut");
  // miscellaneous cuts
  tkDCACut_ = iConfig.getParameter<double>("tkDCACut");
  cosThetaXYCut_ = iConfig.getParameter<double>("cosThetaXYCut");
  cosThetaXYZCut_ = iConfig.getParameter<double>("cosThetaXYZCut");
  
  produces<nanoaod::FlatTable>("meson");
  produces<reco::VertexCompositeCandidateCollection>();
}

void
MesonProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool runOnMC = !iEvent.isRealData();
  
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);
  reco::Vertex pv = recVtxs->at(0);

  Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByToken(jetSrc_, jetHandle);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder_);

  Handle<edm::View<pat::PackedCandidate> > pfHandle;
  iEvent.getByToken(pfSrc_, pfHandle);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder_);

  auto cmVxCand = make_unique<reco::VertexCompositeCandidateCollection>();
  vector<float> dca, angleXY, angleXYZ, trkChi2, trknHits, trkPt, trkipsigXY, trkipsigZ;
  vector<float> lxy, lxySig, l3D, l3DSig, jetDR, legDR, diffMass;
  vector<int> nJet, mcMatch;

  // edm::Handle<edm::View<reco::GenParticle> > mcHandle;  
  // if ( runOnMC ) {
  //   iEvent.getByToken(mcSrc_, mcHandle);
  // }

  int njet = 0;
  // for (const pat::Jet & aPatJet : *jetHandle){
  //   if (aPatJet.pt() < 30 or abs(aPatJet.eta()) > 3 ) continue;

    vector< const pat::PackedCandidate * > jetDaughters, softlepCands;
    
    // for( unsigned int idx = 0 ; idx < aPatJet.numberOfDaughters() ; ++idx) {

    for (auto& dauCand_ : *pfHandle) {
      const pat::PackedCandidate * dauCand(&dauCand_);
      // const pat::PackedCandidate * dauCand ( dynamic_cast<const pat::PackedCandidate*>(aPatJet.daughter(idx)));
      if ( dauCand->charge() == 0 ) continue;
      
      if ( dauCand->pt() < tkPtCut_ ) continue;

      if (dauCand->bestTrack() == nullptr) continue;
      
      const reco::Track * trk = dauCand->bestTrack();      
      if (trk->normalizedChi2() > tkChi2Cut_) continue;

      if (trk->numberOfValidHits() < tkNHitsCut_) continue;

      math::XYZPoint referencePos = pv.position();
      float ipsigXY = std::abs(trk->dxy(referencePos)/trk->dxyError());
      if (ipsigXY > tkIPSigXYCut_) continue;
      
      float ipsigZ = std::abs(trk->dz(referencePos)/trk->dzError());
      if (ipsigZ > tkIPSigZCut_) continue;
      
      jetDaughters.emplace_back( dauCand );
      if ( abs(dauCand->pdgId()) == 11  || abs(dauCand->pdgId())==13) softlepCands.emplace_back(dauCand);
    }
    unsigned int dau_size = jetDaughters.size();


    auto fill = [&](int mcmatch,
		    reco::VertexCompositeCandidate& cand,
		    vector<const pat::PackedCandidate*>& cands,
		    float angleXY_, float angleXYZ_, float dca_) {
      mcMatch.emplace_back(mcmatch);
      cmVxCand->emplace_back(cand);
      dca.emplace_back(dca_);
      angleXY.emplace_back(angleXY_);
      angleXYZ.emplace_back(angleXYZ_);
      nJet.emplace_back(njet);
      
      trackVars tt = getTrackVars(cands, pv);
      trkChi2.emplace_back(tt.normalizedChi2);
      trknHits.emplace_back(tt.nHits);
      trkPt.emplace_back(tt.trkPt);
      trkipsigXY.emplace_back(tt.ipsigXY);
      trkipsigZ.emplace_back(tt.ipsigZ);
      
      auto d2 = getDistance(2,cand,pv);
      lxy.emplace_back(d2.first);
      lxySig.emplace_back(d2.second);
      auto d3 = getDistance(3,cand,pv);	
      l3D.emplace_back(d3.first);
      l3DSig.emplace_back(d3.second);
      
      // jetDR.emplace_back(reco::deltaR( cand, aPatJet));
      jetDR.emplace_back(0);
      legDR.emplace_back(reco::deltaR( *cands[0], *cands[1]));
      diffMass.emplace_back(0);	
    };

    // dont need to sort since all combinations are done.
    //sort(jetDaughters.begin(), jetDaughters.end(), [](const pat::PackedCandidate * a, const pat::PackedCandidate * b) {return a->pt() > b->pt(); }); 

    // if ( dau_size > maxNumPFCand_ ) dau_size = maxNumPFCand_;
    // jetDaughters.resize( dau_size );
    for ( unsigned int lep1_idx = 0 ; lep1_idx< dau_size-1 ; ++lep1_idx) {
      continue;
      const pat::PackedCandidate* lep1Cand = jetDaughters[lep1_idx];
      if ( abs(lep1Cand->pdgId()) != 13 and abs(lep1Cand->pdgId()) != 11) continue;
      
      for ( unsigned int lep2_idx = lep1_idx+1 ; lep2_idx< dau_size ; ++lep2_idx) {
        const pat::PackedCandidate* lep2Cand = jetDaughters[lep2_idx];
	if ( abs(lep2Cand->pdgId()) != 13 and abs(lep2Cand->pdgId()) != 11) continue;

        int pdgMul = lep1Cand->pdgId() * lep2Cand->pdgId();
        if ( pdgMul != -121 and pdgMul != -169 ) continue; 

	float dca_Jpsi = -2;
	float angleXY_Jpsi = -2;
	float angleXYZ_Jpsi = -2;
	vector<const pat::PackedCandidate*> cands_Jpsi{lep1Cand, lep2Cand};
	reco::VertexCompositeCandidate JpsiCand = this->fit(cands_Jpsi, pv, pdgId_Jpsi,
							    dca_Jpsi, angleXY_Jpsi, angleXYZ_Jpsi);

        if ( JpsiCand.mass() < jpsiMin_ || JpsiCand.mass() > jpsiMax_ ) continue;

	int mc_Jpsi = -5;
	if (runOnMC){
	  // mc_Jpsi = findMCmatch(aPatJet, cands_Jpsi, pdgId_Jpsi);
	}
	fill(mc_Jpsi, JpsiCand, cands_Jpsi, angleXY_Jpsi, angleXYZ_Jpsi, dca_Jpsi);
      }
    }

    // if ( applyCuts_ && softlepCands.size()==0  ) continue;
    
    for ( unsigned int pion_idx = 0 ; pion_idx< dau_size ; ++pion_idx) {
      const pat::PackedCandidate * pionCand = jetDaughters[pion_idx];
      if ( abs(pionCand->pdgId()) == 13 or abs(pionCand->pdgId()) == 11) continue;
      
      for ( unsigned int kaon_idx = 0 ; kaon_idx< dau_size ; ++kaon_idx) {
        if ( pion_idx == kaon_idx ) continue;
	pat::PackedCandidate kaonCand(*jetDaughters[kaon_idx]);	
        if ( abs(kaonCand.pdgId()) == 13 or abs(kaonCand.pdgId()) == 11) continue;
        if ( pionCand->charge() * kaonCand.charge() != -1 ) continue;

	// for KS, only want to pick up one version, since we treat
	// both mesons as pions, i.e. dont change identity
	if (kaon_idx < pion_idx) goto lambda;

	{ // jump over variable definitions, so they must be in a separate block
	  float dca_KS = -2;
	  float angleXY_KS = -2;
	  float angleXYZ_KS = -2;
	  vector<const pat::PackedCandidate*> cands_KS{pionCand, &kaonCand};
	  reco::VertexCompositeCandidate KSCand = fit(cands_KS, pv, pdgId_KS,
						      dca_KS, angleXY_KS, angleXYZ_KS);
	  
	  // if ( applyCuts_ ) KSCand.addDaughter( *softlepCands[0] );
	  
	  if ( KSCand.mass() < KSMin_ || KSCand.mass() > KSMax_ ) goto lambda;
	  int mc_KS = -5;
	  if (runOnMC) {
	    // mc_KS = findMCmatch(aPatJet, cands_KS, pdgId_KS);
	  }
	  fill(mc_KS, KSCand, cands_KS, angleXY_KS, angleXYZ_KS, dca_KS);
	}
      lambda:
	kaonCand.setMass(gProtonMass);
	kaonCand.setPdgId(kaonCand.charge()*pdgId_p);
	float dca_Lambda = -2;
	float angleXY_Lambda = -2;
	float angleXYZ_Lambda = -2;
	vector<const pat::PackedCandidate*> cands_Lambda{pionCand, &kaonCand};	
	reco::VertexCompositeCandidate LambdaCand = fit(cands_Lambda, pv, pdgId_Lambda,
							dca_Lambda, angleXY_Lambda, angleXYZ_Lambda);
	
	// if ( applyCuts_ ) LambdaCand.addDaughter( *softlepCands[0] );

	int mc_Lambda = -5;
	if ( LambdaCand.mass() < LambdaMin_ || LambdaCand.mass() > LambdaMax_ ) goto dmeson;
	if (runOnMC) {
	  // if (!doFullMatch_) mc_Lambda = findMiniMCMatch(aPatJet, LambdaCand, cands_Lambda, packed, pruned, pdgId_Lambda);
	  // else mc_Lambda = findMCmatch(aPatJet, cands_Lambda, pdgId_Lambda);
	}
	fill(mc_Lambda, LambdaCand, cands_Lambda, angleXY_Lambda, angleXYZ_Lambda, dca_Lambda);

      dmeson:
	continue;
        kaonCand.setMass(gKaonMass);
	
	float dca_D0 = -2;
	float angleXY_D0 = -2;
	float angleXYZ_D0 = -2;
	vector<const pat::PackedCandidate*> cands_D0{pionCand, &kaonCand};	
	reco::VertexCompositeCandidate D0Cand = fit(cands_D0, pv, pdgId_D0,
						    dca_D0, angleXY_D0, angleXYZ_D0);

        // if ( applyCuts_ ) D0Cand.addDaughter( *softlepCands[0] );

        if ( D0Cand.mass() < D0Min_ || D0Cand.mass() > D0Max_ ) continue;

	int mc_D0 = -5;
	if (runOnMC){
	  // mc_D0 = findMCmatch(aPatJet, cands_D0, pdgId_D0);
	}
	fill(mc_D0, D0Cand, cands_D0, angleXY_D0, angleXYZ_D0, dca_D0);
	
        if ( dau_size < 3 ) continue;

	for( unsigned int extra_pion_idx = 0 ;  extra_pion_idx < dau_size ; ++extra_pion_idx) {
	  if ( extra_pion_idx== pion_idx || extra_pion_idx == kaon_idx) continue;
	  const pat::PackedCandidate * pion2Cand = jetDaughters[extra_pion_idx];
	  if ( abs(pion2Cand->pdgId()) != 211) continue;
	    
	  float dca_Dstar = -2;
	  float angleXY_Dstar = -2;
	  float angleXYZ_Dstar = -2;	  
	  vector<const pat::PackedCandidate*> cands_Dstar{pionCand,&kaonCand,pion2Cand};
	  reco::VertexCompositeCandidate DstarCand = fit(cands_Dstar, pv, pdgId_Dstar,
							 /*pion2Cand->charge()*pdgId_Dstar*/
							 dca_Dstar, angleXY_Dstar, angleXYZ_Dstar);

	  // if ( applyCuts_ ) DstarCand.addDaughter( *softlepCands[0] );
	  
	  float diffMass_Dstar = DstarCand.mass() - D0Cand.mass();
	  if ( diffMass_Dstar < DstarDiffMin_ || diffMass_Dstar > DstarDiffMax_ ) continue;

	  int mc_Dstar = -5;
	  if (runOnMC){
	    // mc_Dstar = findMCmatch(aPatJet, cands_Dstar, pdgId_Dstar);
	  }
	  fill(mc_Dstar, DstarCand, cands_Dstar, angleXY_Dstar, angleXYZ_Dstar, dca_Dstar);
	}
      }
    }
  //   ++njet;
  // }
  
  auto mesonTable = make_unique<nanoaod::FlatTable>(cmVxCand->size(),"meson",false);
  // For SV we fill from here only stuff that cannot be created with the SimpleFlatTableProducer 
  mesonTable->addColumn<float>("dca",dca,"distance of closest approach cm",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<float>("angleXY",angleXY,"2D angle between vertex and tracks",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<float>("angleXYZ",angleXYZ,"3D angle between vertex and tracks",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<int>("nJet",nJet,"nJet of vertex cand",nanoaod::FlatTable::IntColumn);
  mesonTable->addColumn<int>("mcMatch",mcMatch,"mc matching",nanoaod::FlatTable::IntColumn);

  mesonTable->addColumn<float>("trk_normalizedChi2",trkChi2,"trk chi2/ndof",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<int>("trk_nHits",trknHits,"trk nHits",nanoaod::FlatTable::IntColumn);
  mesonTable->addColumn<float>("trk_pt",trkPt,"trk Pt",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<float>("trk_ipsigXY",trkipsigXY,"trk ipsigXY",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<float>("trk_ipsigZ",trkipsigZ,"trk ipsigZ",nanoaod::FlatTable::FloatColumn);

  mesonTable->addColumn<float>("lxy",lxy,"2D decay length in cm",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<float>("lxySig",lxySig,"2D decay length sig in cm",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<float>("l3D",l3D,"3D decay length in cm",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<float>("l3DSig",l3DSig,"3D decay length sig in cm",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<float>("jetDR",jetDR,"DR between jet",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<float>("legDR",legDR,"DR between leg",nanoaod::FlatTable::FloatColumn);
  mesonTable->addColumn<float>("diffMass",diffMass,"diffMass",nanoaod::FlatTable::FloatColumn);
  
  iEvent.put(move(mesonTable),"meson");
  iEvent.put(move(cmVxCand));
  
}

reco::VertexCompositeCandidate MesonProducer::fit(vector<const pat::PackedCandidate*> cands,
						   reco::Vertex pv, int pdgId,
						   float &dca, float &angleXY, float &angleXYZ)
{
  int charge = 0;
  vector<reco::TransientTrack> transientTracks;
  for (auto trk : cands){
    if (trk->bestTrack() == nullptr) continue;
    const TransientTrack transientTrack = trackBuilder_->build(trk->bestTrack());
    transientTracks.emplace_back(transientTrack);
    charge += trk->charge();
  }

  if (transientTracks.size() < 2) return reco::VertexCompositeCandidate();

  // impactPointTSCP DCA
  // measure distance between tracks at their closest approach
  dca = -1;
  ClosestApproachInRPhi cApp;
  reco::TransientTrack tt1 = transientTracks[0], tt2 = transientTracks[1];
  cApp.calculate(tt1.impactPointTSCP().theState(), tt2.impactPointTSCP().theState());  
  if (!cApp.status()) return reco::VertexCompositeCandidate();
  
  dca = cApp.distance();

  if (dca > tkDCACut_) return reco::VertexCompositeCandidate();
  
  // the POCA should at least be in the sensitive volume
  GlobalPoint cxPt = cApp.crossingPoint();
  if (sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) > 120. || abs(cxPt.z()) > 300.)
    return reco::VertexCompositeCandidate();

  // the tracks should at least point in the same quadrant
  TrajectoryStateClosestToPoint const & posTSCP = tt1.trajectoryStateClosestToPoint(cxPt);
  TrajectoryStateClosestToPoint const & negTSCP = tt2.trajectoryStateClosestToPoint(cxPt);
  if (!posTSCP.isValid() || !negTSCP.isValid()) return reco::VertexCompositeCandidate();
  if (posTSCP.momentum().dot(negTSCP.momentum()) < 0) return reco::VertexCompositeCandidate();
  
  KalmanVertexFitter m_kvf(true);
  
  TransientVertex tv = m_kvf.vertex(transientTracks);      
  if (!tv.isValid()) return reco::VertexCompositeCandidate();
  
  reco::Vertex theVtx = tv;
  // loose cut on chi2
  if (theVtx.normalizedChi2() > vtxChi2Cut_) return reco::VertexCompositeCandidate();
  
  GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());

  math::XYZTLorentzVector tlv;
  GlobalVector totalP;
  int i = 0;
  for (auto trk : tv.refittedTracks()){
    TrajectoryStateClosestToPoint const & tscp = trk.trajectoryStateClosestToPoint(vtxPos);
    GlobalVector mom = tscp.momentum();
    double mass = cands[i]->mass();
    double energy = sqrt(mom.mag2() + mass*mass);    
    const math::XYZTLorentzVector lv(mom.x(), mom.y(), mom.z(), energy);
    totalP += mom;
    tlv += lv;
    ++i;
  }
  math::XYZPoint referencePos = pv.position();

  // 2D pointing angle
  double dx = theVtx.x()-referencePos.x();
  double dy = theVtx.y()-referencePos.y();
  double px = totalP.x();
  double py = totalP.y();
  angleXY = (dx*px+dy*py)/(sqrt(dx*dx+dy*dy)*sqrt(px*px+py*py));
  if (angleXY > cosThetaXYCut_) return reco::VertexCompositeCandidate();
  
  // 3D pointing angle
  double dz = theVtx.z()-referencePos.z();
  double pz = totalP.z();
  angleXYZ = (dx*px+dy*py+dz*pz)/(sqrt(dx*dx+dy*dy+dz*dz)*sqrt(px*px+py*py+pz*pz));
  if (angleXYZ > cosThetaXYZCut_) return reco::VertexCompositeCandidate();

  reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());
  const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());

  reco::VertexCompositeCandidate secVert(charge, tlv, vtx, vtxCov, theVtx.chi2(), theVtx.ndof(), pdgId);
  for (auto trk : cands){
    if (trk->bestTrack() == nullptr) continue;
    secVert.addDaughter( *trk );    
  }
  
  return secVert;
}

MesonProducer::trackVars MesonProducer::getTrackVars(vector<const pat::PackedCandidate*> cands,
						       reco::Vertex pv)
{
  float normalizedChi2 = 0;
  int nHits = 100;
  float trkPt = 100;
  float ipsigXY = 100, ipsigXYtemp = 100;
  float ipsigZ = 100, ipsigZtemp = 100;
  math::XYZPoint referencePos = pv.position();
  
  for (auto trk : cands){
    if (trk->bestTrack() == nullptr) continue;
    const reco::Track * rtrk = trk->bestTrack();
    normalizedChi2 = (normalizedChi2 > rtrk->normalizedChi2()) ? normalizedChi2 : rtrk->normalizedChi2();
    nHits = (nHits < rtrk->numberOfValidHits()) ? nHits : rtrk->numberOfValidHits();
    trkPt = (trkPt < rtrk->pt()) ? trkPt : rtrk->pt();
    ipsigXYtemp = std::abs(rtrk->dxy(referencePos)/trk->dxyError());
    ipsigXY = (ipsigXY < ipsigXYtemp) ? ipsigXY : ipsigXYtemp;
    ipsigZtemp = std::abs(rtrk->dz(referencePos)/rtrk->dzError());
    ipsigZ = (ipsigZ < ipsigZtemp) ? ipsigZ : ipsigZtemp;    
  }

  trackVars tt;
  tt.normalizedChi2 = normalizedChi2;
  tt.nHits = nHits;
  tt.trkPt = trkPt;
  tt.ipsigXY = ipsigXY;
  tt.ipsigZ = ipsigZ;
  
  return tt;
}

MesonProducer::SVector3
MesonProducer::getDistanceVector(int dim, reco::VertexCompositeCandidate vertex,reco::Vertex pv)
{
  float z = 0.;
  if (dim == 3) z = vertex.vz() - pv.position().z();
  SVector3 distanceVector(vertex.vx() - pv.position().x(),
			  vertex.vy() - pv.position().y(),
			  z);
  return distanceVector;
}

pair<float, float> MesonProducer::getDistance(int dim, reco::VertexCompositeCandidate vertex,reco::Vertex pv)
{
  SMatrixSym3D totalCov = vertex.vertexCovariance() + pv.covariance();
  SVector3 distVecXYZ = getDistanceVector(dim, vertex, pv);
  float distMagXYZ = ROOT::Math::Mag(distVecXYZ);
  float sigmaDistMagXYZ = sqrt(ROOT::Math::Similarity(totalCov, distVecXYZ)) / distMagXYZ;

  return make_pair(distMagXYZ,sigmaDistMagXYZ);
}

int MesonProducer::findMCmatch(const pat::Jet & aPatJet, vector<const pat::PackedCandidate*> cands, int pdgid)
{
  if (!doMatch_) return -1;
  const reco::JetFlavourInfo & jetInfo =  aPatJet.jetFlavourInfo();
  if (abs(jetInfo.getHadronFlavour()) != 5) return -1;
  if (abs(jetInfo.getPartonFlavour()) != 5) return -1;
  
  const GenParticleRefVector & cHads = jetInfo.getcHadrons();
  int matches = 0;
  for( reco::GenParticleRefVector::const_iterator im = cHads.begin(); im!=cHads.end(); ++im) {
    const reco::GenParticle& part = **im;
    //cout << " cHads " << part.pt() << " " <<part.eta() << " " <<part.pdgId() << " " <<part.status()<< endl;

    if (abs(part.pdgId()) != pdgid) continue;
    
    for( unsigned int idx = 0 ; idx < part.numberOfDaughters() ; ++idx) {
      const reco::Candidate *dauCand = part.daughter(idx);
      // skip non stable particles
      if (dauCand->status() != 1) continue;
      //cout << " mc " << dauCand->pt() << " " <<dauCand->eta() << " " <<dauCand->pdgId() << " " <<dauCand->status() << endl;
      bool match = false;
      for (auto pc : cands){
	//cout << " pc " << pc->pt() << " " <<pc->eta() << " " <<pc->pdgId() << " " <<pc->status()<< endl;
	if (pc->pdgId() != dauCand->pdgId()) continue;
	if (reco::deltaR( *dauCand, *pc) < 0.1){
	  match = true;
	  //cout << " match is " << reco::deltaR( *dauCand, *pc) << endl;
	}
      }
      if (match) matches++;
    }
  }
  return matches;
}

void
MesonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MesonProducer);
