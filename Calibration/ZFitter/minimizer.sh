#! /bin/bash

./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC --autoBin --smearerFit  --targetVariable=ptSum --targetVariable_min=0 --targetVariable_max=200 --targetVariable_binWidth=2 --corrEleType=HggRunEtaR9Et --smearEleType=stochastic &> tmp/debug_ptSum.txt

rm -r tmp/tmpFile-*.root
cp test/dato/fitres/params-ptSum_leading_scaleStep0-Et_25-trigger-noPF.txt output_ptSum_leading.txt
rm test/dato/fitres/histos_ptSum_leading_scaleStep0_Et_25_trigger_noPF.root
./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC --autoBin --smearerFit --initFile=output_ptSum_leading.txt --targetVariable=ptSum --targetVariable_min=0 --targetVariable_max=200 --targetVariable_binWidth=2 --plotOnly --nSmearToy=1 --outDirImgData=test/dato/img --outDirFitResData=test/dato/fitres

root -l -b
gROOT->ProcessLine(".L macro/plot_data_mc.C+");
PlotMeanHist("test/dato/fitres/histos_ptSum_leading_scaleStep0_Et_25_trigger_noPF.root");
gROOT->ProcessLine(".q");

##############################PtRatio

./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC --autoBin --smearerFit --targetVariable=ptRatio --targetVariable_min=0.5 --targetVariable_max=2 --targetVariable_binWidth=0.05 --configuration=leading --corrEleType=HggRunEtaR9Et --smearEleType=stochastic &> tmp/debug_ptRatio_leading.txt

rm -r tmp/tmpFile-*.root
cp test/dato/fitres/params-ptRatio_leading_scaleStep0-Et_25-trigger-noPF.txt output_ptRatio_leading.txt
rm test/dato/fitres/histos_ptRatio_leading_scaleStep0_Et_25_trigger_noPF.root
./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC --autoBin --smearerFit --initFile=output_ptRatio_leading.txt --targetVariable=ptRatio --configuration=leading --targetVariable_min=0.5 --targetVariable_max=2 --targetVariable_binWidth=0.05 --plotOnly --nSmearToy=1 --outDirImgData=test/dato/img --outDirFitResData=test/dato/fitres

root -l -b
gROOT->ProcessLine(".L macro/plot_data_mc.C+");
PlotMeanHist("test/dato/fitres/histos_ptRatio_leading_scaleStep0_Et_25_trigger_noPF.root");
gROOT->ProcessLine(".q");

####
sshfs gfasanel@lxplus6.cern.ch:/afs/cern.ch/user/g/gfasanel/new_version_ECALELF/CMSSW_7_5_0_pre4/src/Calibration/ZFitter/test/dato/img/ Scrivania/remote_dir/ 


