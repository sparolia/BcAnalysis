import FWCore.ParameterSet.Config as cms
from inputFilesList import files_jpsi_munu, files_jpsi_taunu, files_jpsi_plusX
isSigChannel = True
isBkg = False

process = cms.Process("Rootuple")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(wantSummary = (cms.untracked.bool(True))
    )
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000000))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


inputFilesList = []
outputRootFileName = 'RootupleBcToJpsiTrk.root'
decayChannel = 'tau'
if(isSigChannel):
  decayChannel == 'tau'
  inputFilesList = files_jpsi_taunu
  outputRootFileName = 'RootupleBcToJpsiTrk_tauChannel.root'
else:
  decayChannel == 'muon'
  inputFilesList = files_jpsi_munu
  outputRootFileName = 'RootupleBcToJpsiTrk_muonChannel.root'

if(isBkg):
  decayChannel == 'jpsiX'
  inputFilesList = files_jpsi_plusX
  outputRootFileName = 'RootupleBcToJpsiTrk_bkg.root'

#outputRootFileName = 'RootupleBcTo3Mu_pionChannel.root'

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        #files_jpsi_munu
        #files_jpsi_taunu
        #'/store/mc/RunIIAutumn18MiniAOD/BcToJpsiMuNu_JpsiToMuMu_13TeV-TuneCP5-evtgen-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/230000/A9DE66AD-7301-004C-8D87-3FBE21CA63BD.root'
        #'file:mini-aod-4031228.root'
        #'file:miniAOD_99.root' 
        #'file:mini-aod-1829052.root'
        #'/store/data/Run2018A/Charmonium/MINIAOD/17Sep2018-v1/90000/FBB6E58B-3F6C-004A-A1E3-21AB715F7D2B.root'
        #'/store/user/cgalloni/BJpsiX_MuMu_270819/Autumn18_10_2_9_miniAOD/190827_143312/0000/miniAOD_327.root'
        inputFilesList
        #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/230000/FEE73ECF-9377-EA11-97F0-0023AEEEB703.root',
      )
    )

process.triggerSelection = cms.EDFilter('TriggerResultsFilter',
    triggerConditions = cms.vstring('HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*',
      'HLT_Dimuon25_Jpsi_v*',
      'HLT_Dimuon0_Jpsi3p5_Muon2_v*',
      'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*',
      'HLT_DoubleMu4_JpsiTrk_Displaced_v*',
      'HLT_DoubleMu4_Jpsi_Displaced_v*',
      '*'
      ),
    hltResults = cms.InputTag('TriggerResults', '', 'HLT'),
    l1tResults = cms.InputTag(''),
    throw = cms.bool(False)
    )

process.load("BcAnalysis.BcToJpsiTrkAnalyzer.BcToJpsiTrkAnalyzer_cfi")

process.rootuple.isMC = True


process.TFileService = cms.Service('TFileService',
    fileName = cms.string(outputRootFileName)
    )

#process.p = cms.Path(process.triggerSelection*process.slimmedMuonsWithTriggerSequence*process.rootuple)
process.p = cms.Path(process.triggerSelection*process.rootuple)
