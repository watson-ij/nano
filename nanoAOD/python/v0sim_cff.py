import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

##################### Tables for final output and docs ##########################
v0simTable = cms.EDProducer("V0SimProducer",
  moduleLabelTk = cms.InputTag('g4SimHits'),
  moduleLabelVtx = cms.InputTag('g4SimHits'),
)

#after cross linkining
v0simTables = cms.Sequence(v0simTable)
