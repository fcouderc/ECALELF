#! /bin/bash

#source ../initCmsEnvCRAB2.sh

##createOnly -> submitOnly -> check
#option=--createOnly
#option=--submitOnly
option=--check
#To force the hadd of ntuple create the file "finished" in the res directory of your jobs
#touch ...../res/finished and then --chek again
#option=--submit

data(){
##DATI
where=caf
tag=config/reRecoTags/80X_dataRun2_Prompt_v9.py
#where=remoteGlidein #if dataset is at CERN, use caf => much faster
#json=/afs/cern.ch/work/g/gfasanel/80_ECALELF_ntuples_new/CMSSW_8_0_8/src/Calibration/EcalAlCaRecoProducers/json_files_2016/official_jsons/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON.txt
json=/afs/cern.ch/work/g/gfasanel/80_ECALELF_ntuples_new/CMSSW_8_0_8/src/Calibration/EcalAlCaRecoProducers/json_files_2016/diff_10August_4August.txt 
jsonName=GoldenJson_10August
runMin=277166

#

parseDatasetFile.sh alcareco_datasets.dat | grep DoubleEG-Run2016E| grep ${runMin}
#parseDatasetFile.sh alcareco_datasets.dat | grep DoubleEG-Run2016D| grep ${runMin}
#while [ "1" == "1" ];do 
#    #v1 has no runs in the json
   ./scripts/prodNtuples.sh `parseDatasetFile.sh alcareco_datasets.dat | grep DoubleEG-Run2016E| grep ${runMin}` --type MINIAOD -t ${tag} --scheduler=${where} --json=${json} --json_name=${jsonName} $option -s noSkim
#    sleep 5m
#done

}

mc(){
while [ "1" == "1" ];do  
    #where=remoteGlidein
    where=caf
    json=No_json_for_MC
    tag_MC=config/reRecoTags/80X_mcRun2_asymptotic_2016_miniAODv2_v0.py
    ./scripts/prodNtuples.sh `parseDatasetFile.sh alcareco_datasets.dat | grep DYJetsToEE_M-50_LTbinned_5To75` --type MINIAOD -t ${tag_MC} --json=${json} --json_name="LT" --scheduler=${where} ${option} --isMC
    #./scripts/prodNtuples.sh `parseDatasetFile.sh alcareco_datasets.dat | grep DYToEE_NNPDF30_13TeV-powheg-pythia8` --type MINIAOD -t ${tag_MC} --json=${json} --json_name="RunII-2016" --scheduler=${where} ${option} --isMC
    sleep 5m
done    

}

if [[ $1 = "data" ]]; then
    data
fi

if [[ $1 = "mc" ]]; then
    mc
fi


# crab -c prod_ntuples/MINIAODNTUPLE/80X_dataRun2_Prompt_v8/DoubleEG-Run2016B-PromptReco-v2_test/273150-273730/271036-273450_golden/ -submit 1
# crab -c prod_ntuples/MINIAODNTUPLE/80X_dataRun2_Prompt_v8/DoubleEG-Run2016B-PromptReco-v2_test/273150-273730/271036-273450_golden/ -status

#    if [[ $option == "--createOnly" ]]; then
#	crab -c prod_ntuples/MINIAODNTUPLE/80X_mcRun2_asymptotic_2016_v3/DYJets_amcatnlo-RunIISpring16-miniAODv1/allRange/RunII-2016/ -submit 1
#    fi

#    if [[ $option == "--submit" ]]; then
#	crab -c prod_ntuples/MINIAODNTUPLE/80X_mcRun2_asymptotic_2016_v3/DYJets_amcatnlo-RunIISpring16-miniAODv1/allRange/RunII-2016/ -submit 502-688
#    fi


#
# ##Se un job fallisce o ha finito, prenditi l'output
# #getoutput
# crab -c prod_ntuples/MINIAODNTUPLE/74X_dataRun2_Prompt_v4/DoubleEG-emanuele-ZElectron_Run2015D_v3_74X/256584-258158/246908-260627-Prompt_25ns-v1-golden_silver/ -status -getoutput
# # E controlla il res
# /res/