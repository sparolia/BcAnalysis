#! /bin/bash

dirName="jobs/signalCentralProduction2018"
for ((i=1 ; ;i++))
do
  if [[ ! -d "${dirName}_${i}" ]]
  then
    cd /afs/cern.ch/user/g/garamire/work/private/CMSPisa/RJPsiAnalysis/BcTo3MuReconstruction/CMSSW_10_2_9/src
    cmsenv
    cd -
    voms-proxy-init -voms cms
    mkdir -p "${dirName}_${i}"
    echo "${dirName}_${i} has been created as the working directory"
    cp template_crabConfiguration.py "${dirName}_${i}/crabConfiguration.py"
    #cp inputFilesList.py "${dirName}_${i}/inputFilesList.py"
    cd "${dirName}_${i}"
    sed -i 's/XXX_NAME_XXX/RJPsiSignalAnalysis2018Sample/g' crabConfiguration.py
    sed -i 's/XXX_INPUTFILES_XXX/files_centralproduction2018_jpsi_taunu/g' crabConfiguration.py
    sed -i 's/XXX_NJOBS_XXX/86/g' crabConfiguration.py
    sed -i 's/XXX_UNITSPERJOB_XXX/1/g' crabConfiguration.py
    sed -i "s/XXX_INPUTDATASET_XXX/'\/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen\/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v4\/MINIAODSIM'/g" crabConfiguration.py
    sed -i 's/XXX_INPUTDBS_XXX/global/g' crabConfiguration.py
    crab submit crabConfiguration.py
    cd -
    break
  fi
done

