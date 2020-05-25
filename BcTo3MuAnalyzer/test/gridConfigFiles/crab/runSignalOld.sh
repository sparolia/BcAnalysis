#! /bin/bash

dirName="jobs/signalOld"
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
    sed -i 's/XXX_NAME_XXX/RJPsiSignalAnalysisOldSample/g' crabConfiguration.py
    sed -i 's/XXX_INPUTFILES_XXX/files_jpsi_taunu/g' crabConfiguration.py
    sed -i 's/XXX_NJOBS_XXX/24/g' crabConfiguration.py
    sed -i 's/XXX_UNITSPERJOB_XXX/10/g' crabConfiguration.py
    sed -i "s/XXX_INPUTDATASET_XXX/'\/BcJpsiTauNu_020519\/cgalloni-Fall18_10_2_9-MINIAODSIM_noDuplCheck_020519-092bfc61e82f18935ea11e32077a486f\/USER'/g" crabConfiguration.py
    sed -i 's/XXX_INPUTDBS_XXX/phys03/g' crabConfiguration.py
    crab submit crabConfiguration.py
    cd -
    break
  fi
done

