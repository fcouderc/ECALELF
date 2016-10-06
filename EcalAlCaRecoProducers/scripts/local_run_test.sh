#
#tagMC=config/reRecoTags/80X_mcRun2_asymptotic_2016_v3.py 
#file="root://xrootd.unl.edu//store/mc/RunIISpring16MiniAODv1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/008099CA-5501-E611-9AAE-24BE05BDCEF1.root"
##file="root://xrootd.unl.edu//store/mc/RunIISpring16MiniAODv1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/20000/0017320C-7BFC-E511-9B2D-0CC47A4C8E34.root"
#cmsRun python/alcaSkimming.py type=MINIAODNTUPLE files=file:$file maxEvents=10000 isCrab=0 tagFile=${tagMC} MC=1

#data
#json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-273450_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt
#json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt
#json=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt
#json=/afs/cern.ch/work/g/gfasanel/80_ECALELF_ntuples_new/CMSSW_8_0_8/src/Calibration/EcalAlCaRecoProducers/json_test_original.txt
tag_data=config/reRecoTags/80X_dataRun2_Prompt_v8.py
tag_MC=config/reRecoTags/80X_mcRun2_asymptotic_2016_miniAODv2_v0.py
#file="root://xrootd.unl.edu//store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/0227DB1C-E719-E611-872C-02163E0141F9.root"
#file="root://cms-xrd-global.cern.ch//store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v1/000/272/762/00000/361EFD18-C015-E611-A822-02163E0142F0.root"
file="root://cms-xrd-global.cern.ch//store/mc/RunIISpring16MiniAODv2/DYJetsToEE_M-50_LTbinned_5To75_5f_LO_13TeV-madgraph_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/02138B87-154E-E611-89DD-0CC47A745282.root"
cmsRun python/alcaSkimming.py type=MINIAODNTUPLE files=file:$file maxEvents=1000 isCrab=0 tagFile=${tag_MC} jsonFile=${json} doTree=1 doTreeOnly=1 isPrivate=0 MC=0 skim=
#if doTree=0 => nojsonFile applied
#type=MINIAODNTUPLE doTree=1 doTreeOnly=1  jsonFile=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274240_13TeV_PromptReco_Coll\
#isions16_JSON.txt isCrab=1 skim= tagFile=config/reRecoTags/80X_dataRun2_Prompt_v8.py isPrivate=0 MC=0
