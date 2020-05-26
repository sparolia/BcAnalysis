# BcAnalysis
Code within the CMSSW framework for the recontruction of Bc+ candidates going to J/Psi + lepton + nu.

## Instructions to run the code:
```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_2_9
cd CMSSW_10_2_9/src/
cmsenv
voms-proxy-init -voms cms -valid 192:00
git clone https://github.com/gabriel-rmrz/BcAnalysis.git RJPsiAnalyzers
scram b
cd RJPsiAnalyzers/BcTo3MuAnalyzer/test
cmsRun BcTo3MuAnalyzer_miniAOD.py

```
