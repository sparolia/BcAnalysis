from WMCore.Configuration import Configuration
#from inputFilesList import files_jpsi_munu, files_jpsi_taunu, files_jpsi_plusX, files_centralproduction2017_jpsi_taunu, files_centralproduction2018_jpsi_taunu
config = Configuration()

config.section_("General")
#config.General.instance = 'preprod'
config.General.requestName = 'XXX_NAME_XXX'
config.General.workArea = 'XXX_NAME_XXX' 
config.General.transferLogs    = True
config.General.transferOutputs = True


config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '/afs/cern.ch/user/g/garamire/work/private/CMSPisa/RJPsiAnalysis/BcTo3MuReconstruction/CMSSW_9_4_17/src/RJPsiAnalyzers/BcTo3MuAnalyzer/test/gridConfigFiles/crab/BcTo3MuAnalyzer_miniAOD.py'

config.section_("Data")
config.Data.publication  = False
config.Data.inputDBS = 'XXX_INPUTDBS_XXX'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = XXX_UNITSPERJOB_XXX 
NJOBS = XXX_NJOBS_XXX
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
#config.Data.userInputFiles = XXX_INPUTFILES_XXX
config.Data.inputDataset = XXX_INPUTDATASET_XXX
config.Data.outLFNDirBase = '/store/user/garamire'
config.Data.outputDatasetTag = 'outputDatasetTag'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.whitelist = ['T2_BE_*','T2_IT_*','T2_DE_*','T2_FR_*','T2_ES_*','T2_UK_*']
config.Site.storageSite = 'T2_IT_Pisa'

from datetime import datetime as dt
submitdate = dt.now().strftime('%Y%m%d')

