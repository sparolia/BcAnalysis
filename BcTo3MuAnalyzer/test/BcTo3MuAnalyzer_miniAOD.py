import FWCore.ParameterSet.Config as cms
from inputFilesList import files_jpsi_munu
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
process.options = cms.untracked.PSet(wantSummary = (cms.untracked.bool(True)),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5000)) 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        #files_jpsi_munu
        '/store/mc/RunIIFall17MiniAODv2/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/40000/DA10C116-4FA9-E811-9B53-0CC47A78A426.root'
        #'/store/mc/RunIIFall18GS/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/GEN-SIM/102X_upgrade2018_realistic_v11-v2/230000/AE31B611-FB1A-334F-A8B9-9FDDA8405DFC.root'
        #'/store/mc/RunIIFall18pLHE/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/LHE/102X_upgrade2018_realistic_v11-v2/230000/1C5385EB-E947-BB43-B2BF-D7C7C471988D.root'
        #'/store/data/Run2018D/Charmonium/MINIAOD/PromptReco-v2/000/325/022/00000/B4173414-9341-6442-8A8C-39C79E9FE09B.root'
        #'/store/data/Run2018D/Charmonium/MINIAOD/PromptReco-v2/000/324/835/00000/116C58C4-7F69-594D-A988-2EAF8AD47266.root'
        #'/store/data/Run2018C/Charmonium/MINIAOD/PromptReco-v2/000/319/756/00000/EEF6CEC1-698B-E811-8081-02163E00AF5F.root',
      )

    )

process.triggerSelection = cms.EDFilter('TriggerResultsFilter',
    triggerConditions = cms.vstring('HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*',
      'HLT_Dimuon25_Jpsi_v*',
      'HLT_Dimuon0_Jpsi3p5_Muon2_v*',
      #'HLT_DoubleMu4_3_Jpsi_Displaced_v*',
      'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*',
      'HLT_DoubleMu4_JpsiTrk_Displaced_v*',
      'HLT_DoubleMu4_Jpsi_Displaced_v*' 
      ),
    hltResults = cms.InputTag('TriggerResults', '', 'HLT'),
    l1tResults = cms.InputTag(''),
    throw = cms.bool(False)
    )

#TODO: Add my own modueles. 
process.load("RJPsiAnalyzers.BcTo3MuAnalyzer.slimmedMuonsTriggerMatcher_cfi")
process.load("RJPsiAnalyzers.BcTo3MuAnalyzer.BcTo3MuAnalyzer_cfi")
#process.rootuple.dimuons = cms.InputTag('slimmedMuonsWithTrigger') 

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('RootupleBcTo3Mu_miniAOD.root')
    )

process.p = cms.Path(process.triggerSelection*process.slimmedMuonsWithTriggerSequence*process.rootuple)
#process.p = cms.Path(process.triggerSelection*process.rootuple)
#process.p = cms.Path(process.triggerSelection*process.slimmedMuonsWithTriggerSequence)
