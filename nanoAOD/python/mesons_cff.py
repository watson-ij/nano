import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################### Tables for final output and docs ##########################
mesonTable = cms.EDProducer("MesonProducer",
  jetLabel = cms.InputTag("slimmedJets"),
  pfLabel = cms.InputTag("packedPFCandidates"),
  vertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
  mcLabel  = cms.InputTag("prunedGenParticles"),
  applySoftLeptonCut = cms.bool(True),
  doMCMatching = cms.bool(False),
  # -- cuts on initial track collection --
  # Track normalized Chi2 <
  tkChi2Cut = cms.double(100),
  # Number of valid hits on track >=
  tkNHitsCut = cms.int32(0),
  # Pt of track >
  tkPtCut = cms.double(0.),
  # Track impact parameter significance >
  tkIPSigXYCut = cms.double(100000000000),
  tkIPSigZCut = cms.double(100000000000),
  
  # -- cuts on the vertex --
  # Vertex chi2 <
  vtxChi2Cut = cms.double(10),
  # XY decay distance significance >
  vtxDecaySigXYCut = cms.double(-1),
  # XYZ decay distance significance >
  vtxDecaySigXYZCut = cms.double(-1.),
  
  # -- miscellaneous cuts --
  # POCA distance between tracks <
  tkDCACut = cms.double(100),
  # cos(angleXY) between x and p of V0 candidate >
  cosThetaXYCut = cms.double(100),
  # cos(angleXYZ) between x and p of V0 candidate >
  cosThetaXYZCut = cms.double(100),
)

mesonCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("mesonTable"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("meson"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(True), 
    variables = cms.PSet(P4Vars,
        x   = Var("vx()", float, doc = "secondary vertex X position, in cm",precision=10),
        y   = Var("vy()", float, doc = "secondary vertex Y position, in cm",precision=10),
        z   = Var("vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
        chi2= Var("vertexChi2()", float, doc = "chi2",precision=14),
        ndof= Var("vertexNdof()", int, doc = "number of degrees of freedom"),
        pdgId=Var("pdgId()", int, doc = "pdgId"),
    ),
)
mesonCandidateTable.variables.pt.precision=14
mesonCandidateTable.variables.phi.precision=14
mesonCandidateTable.variables.eta.precision=14
mesonCandidateTable.variables.mass.precision=14

#before cross linking
mesonSequence = cms.Sequence()
#after cross linkining
mesonTables = cms.Sequence(mesonTable+mesonCandidateTable)

