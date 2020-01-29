from WMCore.Configuration import Configuration
from inputFilesList import files_jpsi_munu, files_jpsi_taunu
config = Configuration()

config.section_("General")
config.General.requestName = 'RJPsiNormalizationAnalysis'
config.General.workArea = 'RJPsiNormalizationAnalysis'
config.General.transferLogs    = True
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
#config.JobType.psetName    = '/gpfs/ddn/users/sanchez/RPJpsi/CMSSW_10_2_9/src/RJPsiAnalyzers/BcTo3MuAnalyzer/test/gridConfigFiles/crab/BcTo3MuAnalyzer_miniAOD.py'
config.JobType.psetName    = '/afs/cern.ch/user/g/garamire/work/private/CMSPisa/RJPsiAnalysis/BcTo3MuReconstruction/CMSSW_10_2_9/src/RJPsiAnalyzers/BcTo3MuAnalyzer/test/gridConfigFiles/crab/BcTo3MuAnalyzer_miniAOD.py'

config.section_("Data")
config.Data.publication  = False
config.Data.inputDBS = 'global'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
NJOBS = 24 
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.userInputFiles = files_jpsi_munu 
#config.Data.outputPrimaryDataset = 'outputPrimaryDataset'
#config.Data.outLFNDirBase = '/store/group/dpg_rpc/comm_rpc/Sandbox/garamire'
#config.Data.outLFNDirBase = '/store/group/phys_bphys/garamire/'
#config.Data.outLFNDirBase = 'gsiftp://eosuserftp.cern.ch/eos/user/g/garamire'
config.Data.outLFNDirBase = '/store/user/garamire'
config.Data.outputDatasetTag = 'outputDatasetTag'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.whitelist = ['T2_BE_*','T2_IT_*','T2_DE_*','T2_FR_*','T2_ES_*','T2_UK_*']
config.Site.storageSite = 'T2_IT_Pisa'

from datetime import datetime as dt
submitdate = dt.now().strftime('%Y%m%d')

