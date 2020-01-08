import FWCore.ParameterSet.Config as cms
from inputFilesList import files_jpsi_munu, files_jpsi_taunu
# dacayChannel= 'muon'
decayChannel= 'tau'


process = cms.Process("Rootuple")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
#process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(wantSummary = (cms.untracked.bool(True))
    )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5000))


inputFilesList = []
outputRootFileName = 'RootupleBcTo3Mu.root'
if(decayChannel == 'muon'):
  inputFilesList = files_jpsi_munu
  outputRootFileName = 'RootupleBcTo3Mu_muonChannel.root'
elif(decayChannel == 'tau'):
  inputFilesList = files_jpsi_taunu
  outputRootFileName = 'RootupleBcTo3Mu_tauChanel.root'

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        #files_jpsi_munu
        #files_jpsi_taunu
        inputFilesList
      )
    )

process.triggerSelection = cms.EDFilter('TriggerResultsFilter',
    triggerConditions = cms.vstring('HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*',
      'HLT_Dimuon25_Jpsi_v*',
      'HLT_Dimuon0_Jpsi3p5_Muon2_v*',
      'HLT_DoubleMu4_3_Jpsi_v*',
      'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*',
      'HLT_DoubleMu4_JpsiTrk_Displaced_v*',
      'HLT_DoubleMu4_Jpsi_Displaced_v*'
      ),
    hltResults = cms.InputTag('TriggerResults', '', 'HLT'),
    l1tResults = cms.InputTag(''),
    throw = cms.bool(False)
    )

process.load("RJPsiAnalyzers.BcTo3MuAnalyzer.BcTo3MuAnalyzer_cfi")

process.TFileService = cms.Service('TFileService',
    fileName = cms.string(outputRootFileName)
    )

#process.p = cms.Path(process.triggerSelection*process.slimmedMuonsWithTriggerSequence*process.rootuple)
process.p = cms.Path(process.triggerSelection*process.rootuple)
