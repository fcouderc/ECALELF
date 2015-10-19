#!/bin/bash
#commonCut=Et_25-trigger-noPF #used in the category names
regionsFile=scaleStep0
regionDir=test/dato/fitres/toys/${regionsFile}

for scale in  0.98  0.985 0.99 0.995 1.00 1.005 1.01 1.015 1.02
do
    for constAlpha in  0.00-0.00 0.005-0.00 0.01-0.00 0.0105-0.00 0.02-0.00
#for scale in 0.99
#do
#    for constAlpha in  0.00-0.00
    do
	for det_region in EB EE
	do
	    const=`echo $constAlpha | cut -d '-' -f 1`
	    alpha=`echo $constAlpha | cut -d '-' -f 2`
	
	    echo "*********Fitting const:alpha:scale ${const}:${alpha}:${scale}"
	    workDir=${regionDir}/${const}_${alpha}_${scale}/${det_region}/
	    
	    ./script/haddTGraph.sh -o ${workDir}outProfile_ptRatio_pt2Sum_random_${scale}_${const}_${regionsFile}_Et_25_trigger_noPF.root ${workDir}*/outProfile_ptRatio_pt2Sum_random_${regionsFile}_Et_25_trigger_noPF.root # with this you have a root file with ALL the likelihood scan

	    echo "[STATUS] Making Likelihood plots for closure tests"
	    echo "Region is ${det_region}"
	    echo "{" > ${workDir}fitProfiles.C
	    echo "gROOT->ProcessLine(\".include $ROOFITSYS/include\");" >> ${workDir}fitProfiles.C
	    echo "gROOT->ProcessLine(\".L macro/macro_fit.C+\");" >> ${workDir}fitProfiles.C
	    echo "gROOT->ProcessLine(\".L macro/plot_data_mc.C+\");" >> ${workDir}fitProfiles.C
	    echo "FitProfile2(\"${workDir}outProfile_ptRatio_pt2Sum_random_${scale}_${const}_${regionsFile}_Et_25_trigger_noPF.root\",\"\",\"\",true,true,true);" >> ${workDir}fitProfiles.C
	    echo "}" >> ${workDir}/fitProfiles.C
	    root -l -b -q ${workDir}/fitProfiles.C
	done #barrel,endcap for
    done #constAlpha for
done #scale for

cd ${regionDir}
#ZFitter/test/dato/fitres/toys/scaleStep0
source fit_copy.sh
cd -
#Per mettere i fit nella www
