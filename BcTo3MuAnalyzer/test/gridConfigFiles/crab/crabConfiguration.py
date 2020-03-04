from WMCore.Configuration import Configuration
from inputFilesList import files_jpsi_munu, files_jpsi_taunu, files_jpsi_plusX
config = Configuration()

config.section_("General")
#config.General.requestName = 'RJPsiNormalizationAnalysis'
#config.General.workArea = 'RJPsiNormalizationAnalysis'
#config.General.requestName = 'RJPsiSignalAnalysisCentralProd'
#config.General.workArea = 'RJPsiSignalAnalysisCentralProd'
config.General.requestName = 'RJPsiBackground1AnalysisBigSample'
config.General.workArea = 'RJPsiBackground1AnalysisBigSample'
config.General.transferLogs    = True
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '/afs/cern.ch/user/g/garamire/work/private/CMSPisa/RJPsiAnalysis/BcTo3MuReconstruction/CMSSW_10_2_9/src/RJPsiAnalyzers/BcTo3MuAnalyzer/test/gridConfigFiles/crab/BcTo3MuAnalyzer_miniAOD.py'

config.section_("Data")
config.Data.publication  = False
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 200
NJOBS = 26
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
#config.Data.userInputFiles = files_jpsi_munu 
#config.Data.userInputFiles = files_jpsi_taunu 
#config.Data.userInputFiles = files_jpsi_plusX
#config.Data.inputDataset = '/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
config.Data.inputDataset = '/BJpsiX_MuMu_270819/cgalloni-Autumn18_10_2_9_miniAOD-39a089a8e7301f392b8b059e430f83ef/USER'
config.Data.inputDBS = 'phys03'
#config.Data.outputPrimaryDataset = 'outputPrimaryDataset'
config.Data.outLFNDirBase = '/store/user/garamire'
config.Data.outputDatasetTag = 'outputDatasetTag'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.whitelist = ['T2_BE_*','T2_IT_*','T2_DE_*','T2_FR_*','T2_ES_*','T2_UK_*']
config.Site.storageSite = 'T2_IT_Pisa'

from datetime import datetime as dt
submitdate = dt.now().strftime('%Y%m%d')

