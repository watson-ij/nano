import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

kshortCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("slimmedKshortVertices"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("Kshort"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), 
    variables = cms.PSet(P4Vars,
        x   = Var("vx()", float, doc = "secondary vertex X position, in cm",precision=10),
        y   = Var("vy()", float, doc = "secondary vertex Y position, in cm",precision=10),
        z   = Var("vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
        chi2= Var("vertexChi2()", float, doc = "chi2",precision=14),
        ndof= Var("vertexNdof()", int, doc = "number of degrees of freedom"),
        pdgId=Var("pdgId()", int, doc = "pdgId"),
    ),
)
kshortCandidateTable.variables.pt.precision=14
kshortCandidateTable.variables.phi.precision=14
kshortCandidateTable.variables.eta.precision=14
kshortCandidateTable.variables.mass.precision=14

lambdaCandidateTable =  cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("slimmedLambdaVertices"),
    cut = cms.string(""),  #DO NOT further cut here, use vertexTable.svCut
    name = cms.string("Lambda"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), 
    variables = cms.PSet(P4Vars,
        x   = Var("vx()", float, doc = "secondary vertex X position, in cm",precision=10),
        y   = Var("vy()", float, doc = "secondary vertex Y position, in cm",precision=10),
        z   = Var("vz()", float, doc = "secondary vertex Z position, in cm",precision=14),
        chi2= Var("vertexChi2()", float, doc = "chi2",precision=14),
        ndof= Var("vertexNdof()", int, doc = "number of degrees of freedom"),
        pdgId=Var("pdgId()", int, doc = "pdgId"),
    ),
)
lambdaCandidateTable.variables.pt.precision=14
lambdaCandidateTable.variables.phi.precision=14
lambdaCandidateTable.variables.eta.precision=14
lambdaCandidateTable.variables.mass.precision=14

#before cross linking
v0Sequence = cms.Sequence()
#after cross linkining
v0Tables = cms.Sequence(kshortCandidateTable+lambdaCandidateTable)
