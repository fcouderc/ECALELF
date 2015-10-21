#!/bin/bash
source script/functions.sh

if [[ $1 = "ptRatio*pt2Sum" ]]; then
   validation_file=22Jan2012-runDepMCAll_checkMee.dat
   region=scaleStep0
   commonCut=Et_25-trigger-noPF
   outDirData=test/dato
   extension=ptRatio_pt2Sum
   InitFile=""
   Target=ptRatio*pt2Sum
   Min=0.5*0
   Max=2*200
   BinWidth=0.05*2
   Conf=random
   echo ${validation_file}
   echo ${region}
   echo ${extension}
fi

if [[ $1 = "ptRatio*pt2Sum*noCorr" ]]; then
   validation_file=22Jan2012-runDepMCAll_checkMee.dat
   region=scaleStep0
   commonCut=Et_25-trigger-noPF
   outDirData=test/dato
   extension=ptRatio_pt2Sum_withoutCorr
   InitFile=""
   Target=ptRatio*pt2Sum
   Min=0.5*0
   Max=2*200
   BinWidth=0.05*2
   Conf=random
   echo ${validation_file}
   echo ${region}
   echo ${extension}
fi

if [[ $1 = "ptRatio" ]]; then
   validation_file=22Jan2012-runDepMCAll_checkMee.dat
   region=scaleStep0
   commonCut=Et_25-trigger-noPF
   outDirData=test/dato
   extension=ptRatio
   InitFile=""
   Target=ptRatio
   Min=0.5
   Max=2
   BinWidth=0.05
   Conf=leading
   echo ${validation_file}
   echo ${region}
   echo ${extension}
fi

if [[ $1 = "ptRatio_random" ]]; then
   validation_file=22Jan2012-runDepMCAll_checkMee.dat
   region=scaleStep0
   commonCut=Et_25-trigger-noPF
   outDirData=test/dato
   extension=ptRatio_random
   InitFile=""
   Target=ptRatio
   Min=0.5
   Max=2
   BinWidth=0.05
   Conf=random
   echo ${validation_file}
   echo ${region}
   echo ${extension}
fi

if [[ $1 = "ptSum" ]]; then
   validation_file=22Jan2012-runDepMCAll_checkMee.dat
   region=scaleStep0
   commonCut=Et_25-trigger-noPF
   outDirData=test/dato
   extension=ptSum
   InitFile=""
   Target=ptSum
   Min=0
   Max=200
   BinWidth=2
   Conf=leading
   echo ${validation_file}
   echo ${region}
   echo ${extension}
fi

if [[ $1 = "" ]]; then
    echo "[ERROR] You must specify the targetVariable"
    exit 0
fi

#Be careful with onlyScale: you must use a NON ZERO SIGMA!

# you should run this script using screen 
# it creates a single job of 50 subjobs (array of jobs)
# Directory organization: 

 if [ ! -e "${outDirData}/${extension}/fitres" ];then mkdir ${outDirData}/${extension}/fitres -p; fi
 if [ ! -e "${outDirData}/${extension}/img" ];then mkdir ${outDirData}/${extension}/img -p; fi

############# SUBMITTING JOBS ###########################################

 for index in `seq 1 50` #Making directories
 do
     mkdir ${outDirData}/${extension}/${index}/fitres/ -p
     mkdir ${outDirData}/${extension}/${index}/img -p
 done
 
# echo ./bin/ZFitter.exe -f data/validation/${validation_file} --regionsFile=data/regions/${region}.dat --invMass_var=invMass_SC_regrCorrSemiParV5_ele --autoBin --smearerFit --targetVariable=${Target} --targetVariable_min=${Min} --targetVariable_max=${Max} --targetVariable_binWidth=${BinWidth} --configuration=${Conf} --corrEleType=HggRunEtaR9Et --smearEleType=stochastic --outDirImgData=${outDirData}/${extension}/\$LSB_JOBINDEX/img/ --outDirFitResData=${outDirData}/${extension}/\$LSB_JOBINDEX/fitres/ ${InitFile}

 if [[ $1 = "ptRatio*pt2Sum*noCorr" ]]; then 
     bsub -q 2nd\
            -oo ${outDirData}/${extension}/%I/fitres/${Target}-${region}-stdout.log\
            -eo ${outDirData}/${extension}/%I/fitres/${Target}-${region}-stderr.log\
            -J "${region} ${extension}[1-50]"\
            "cd $PWD; eval \`scramv1 runtime -sh\`; uname -a; echo \$CMSSW_VERSION;

./bin/ZFitter.exe -f data/validation/${validation_file} --regionsFile=data/regions/${region}.dat --invMass_var=invMass_SC_regrCorrSemiParV5_ele --autoBin --smearerFit --targetVariable=${Target} --targetVariable_min=${Min} --targetVariable_max=${Max} --targetVariable_binWidth=${BinWidth} --configuration=${Conf} --outDirImgData=${outDirData}/${extension}/\$LSB_JOBINDEX/img/ --outDirFitResData=${outDirData}/${extension}/\$LSB_JOBINDEX/fitres/ ${InitFile} || exit 1";
 else
     bsub -q 2nd\
            -oo ${outDirData}/${extension}/%I/fitres/${Target}-${region}-stdout.log\
            -eo ${outDirData}/${extension}/%I/fitres/${Target}-${region}-stderr.log\
            -J "${region} ${extension}[1-50]"\
            "cd $PWD; eval \`scramv1 runtime -sh\`; uname -a; echo \$CMSSW_VERSION;
     ./bin/ZFitter.exe -f data/validation/${validation_file} --regionsFile=data/regions/${region}.dat --invMass_var=invMass_SC_regrCorrSemiParV5_ele --autoBin --smearerFit --targetVariable=${Target} --targetVariable_min=${Min} --targetVariable_max=${Max} --targetVariable_binWidth=${BinWidth} --configuration=${Conf} --corrEleType=HggRunEtaR9Et --smearEleType=stochastic --outDirImgData=${outDirData}/${extension}/\$LSB_JOBINDEX/img/ --outDirFitResData=${outDirData}/${extension}/\$LSB_JOBINDEX/fitres/ ${InitFile} || exit 1"; 
 fi
 

exit 0
#Now, wait the end of your jobs
	
while [ "`bjobs -J \"${region} ${extension}\" | grep -v JOBID | grep -v found | wc -l`" != "0" ]; do /bin/sleep 2m; done

echo "job completed!"

##copiarci Likelihood fitter qui di seguito
./script/haddTGraph.sh -o ${outDirData}/${extension}/fitres/outProfile-eop-$region-${commonCut}.root ${outDirData}/${extension}/*/fitres/outProfile-eop-$region-${commonCut}.root

    