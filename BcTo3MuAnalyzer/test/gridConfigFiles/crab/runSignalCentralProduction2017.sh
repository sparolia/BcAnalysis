#! /bin/bash

dirName="jobs/signalCentralProduction2017"
for ((i=1 ; ;i++))
do
  if [[ ! -d "${dirName}_${i}" ]]
  then
    cd /afs/cern.ch/user/g/garamire/work/private/CMSPisa/RJPsiAnalysis/BcTo3MuReconstruction/CMSSW_9_4_17/src
    cmsenv
    cd -
    voms-proxy-init -voms cms
    mkdir -p "${dirName}_${i}"
    echo "${dirName}_${i} has been created as the working directory"
    cp template_crabConfiguration.py "${dirName}_${i}/crabConfiguration.py"
    cp inputFilesList.py "${dirName}_${i}/inputFilesList.py"
    cd "${dirName}_${i}"
    sed -i 's/XXX_NAME_XXX/RJPsiSignalAnalysis2017Sample/g' crabConfiguration.py
    sed -i 's/XXX_INPUTFILES_XXX/files_centralproduction2017_jpsi_taunu/g' crabConfiguration.py
    sed -i 's/XXX_NJOBS_XXX/10/g' crabConfiguration.py
    sed -i 's/XXX_UNITSPERJOB_XXX/10/g' crabConfiguration.py
    sed -i "s/XXX_INPUTDATASET_XXX/'\/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen\/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2\/MINIAODSIM'/g" crabConfiguration.py
    sed -i 's/XXX_INPUTDBS_XXX/global/g' crabConfiguration.py
    crab submit crabConfiguration.py
    cd -
    break
  fi
done

