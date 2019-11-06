import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'file:miniAOD_99.root'
      )
    )

process.triggerSelection = cms.EDFilter('TriggerResultsFilter',
    triggerConditions = cms.vstring('HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*',
      'HLT_Dimuon25_Jpsi_v*',
      'HLT_Dimuon25_Jpsi_v*',
      'HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*',
      'HLT_Dimuon0_Jpsi3p5_Muon2_v*',
      ),
    hltResults = cms.InputTag('TriggerResults', '', 'HLT'),
    l1tResults = cms.InputTag(''),
    throw = cms.bool(False)
    )

#TODO: Add my own modueles. 
process.load("RJPsiAnalyzers.BcTo3MuAnalyzer.slimmedMuonsTriggerMatcher_cfi")
process.load("RJPsiAnalyzers.BcTo3MuAnalyzer.BcTo3MuAnalyzer_cfi")

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('RootupleBcTo3Mu_miniAOD.root')
    )

process.p = cms.Path(process.triggerSelection*process.slimmedMuonsWithTriggerSequence*process.rootuple)
#process.p = cms.Path(process.triggerSelection*process.slimmedMuonsWithTriggerSequence)
