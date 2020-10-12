import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer("BcToJpsiTrkAnalyzer",
    dimuons = cms.InputTag('slimmedMuons'),
    Trak = cms.InputTag('packedPFCandidates'),
    prunedGenParticles = cms.InputTag('prunedGenParticles'),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    genPUProtons = cms.InputTag("genPUProtons"),
    primaryVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    bslabel = cms.InputTag('offlineBeamSpot'),
    triggerResults = cms.InputTag('TriggerResults', '','HLT'),
    triggerPrescales = cms.InputTag('patTrigger'),
    OnlyBest = cms.bool(False),
    isMC = cms.bool(True),
    OnlyGen = cms.bool(False),
    )
 
