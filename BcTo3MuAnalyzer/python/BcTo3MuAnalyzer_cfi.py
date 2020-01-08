import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer("BcTo3MuAnalyzer",
    dimuons = cms.InputTag('slimmedMuons'),
    Trak = cms.InputTag('packedPFCandidates'),
    prunedGenParticles = cms.InputTag('prunedGenParticles'),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    genPUProtons = cms.InputTag("genPUProtons"),
    primaryVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    bslabel = cms.InputTag('offlineBeamSpot'),
    TriggerResults = cms.InputTag('TriggerResults', '','HLT'),
    OnlyBest = cms.bool(False),
    isMC = cms.bool(True),
    OnlyGen = cms.bool(False),
    )
 
