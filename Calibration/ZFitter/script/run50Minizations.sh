#!/bin/bash
source script/functions.sh

check_EoP=1
EoP_gain=2
EoP_single=3
EoP_gain_extra=4


if [ $1 -eq $check_EoP ]; then
   echo "check_EoP"
   validation_file=22Jan2012-runDepMCAll_v5.dat
   region=scaleStep4smearing_0a_EoP
   EoP_option=true
   commonCut=Et_25-trigger-noPF
   outDirData=test/dato
   extension=Double_Ele_scaleStep4_smearing_0a
   echo ${validation_file}
   echo ${region}
   echo ${EoP_option}
   echo ${extension}
fi


#: << 'MYCOMMENT'
if [ $1 -eq $EoP_gain ]; then
   echo "EoP_gain";
   validation_file=22Jan2012-runDepMCAll_v5.dat
   region=gainSwitch_6
   EoP_option=true
   commonCut=Et_25-trigger-noPF
   outDirData=test/dato
   extension=Double_Ele_gainSwitch_6
   initFile=--initFile=init_eop_gainSwitch_6.dat
   echo ${validation_file}
   echo ${region}
   echo ${EoP_option}
   echo ${extension}
fi
#MYCOMMENT


if [ $1 -eq $EoP_single ]; then
   echo "EoP_single";
   validation_file=22Jan2012-runDepMCAll_single.dat
   region=gainSwitch_7
   EoP_option=true
   commonCut=Et_25-trigger-noPF
   outDirData=test/dato
   extension=single_Ele_gainSwitch_7
   initFile=--initFile=init_eop_single.dat
   echo ${validation_file}
   echo ${region}
   echo ${EoP_option}
   echo ${extension}
fi


if [ $1 -eq $EoP_gain_extra ]; then
   echo "EoP_gain_extra";
   validation_file=22Jan2012-runDepMCAll_v5.dat
   region=gainSwitch_8
   EoP_option=true
   commonCut=Et_25-trigger-noPF
   outDirData=test/dato
   extension=Double_Ele_gainSwitch_8
   initFile=--initFile=init_eop_gainSwitch_8.dat
   echo ${validation_file}
   echo ${region}
   echo ${EoP_option}
   echo ${extension}
fi


#Be careful with onlyScale: you must use a NON ZERO SIGMA!

# you should run this script using screen 
# it creates a single job of 50 subjobs (array of jobs)
# Directory organization: 

############################My EoP case####################################

 if [ ! -e "${outDirData}/${extension}/fitres" ];then mkdir ${outDirData}/${extension}/fitres -p; fi
 if [ ! -e "${outDirData}/${extension}/img" ];then mkdir ${outDirData}/${extension}/img -p; fi

############# SUBMITTING JOBS ###########################################

        for index in `seq 1 50`
        do
            mkdir ${outDirData}/${extension}/${index}/fitres/ -p
            mkdir ${outDirData}/${extension}/${index}/img -p
        done

echo ./bin/ZFitter.exe -f data/validation/${validation_file} --corrEleType=HggRunEtaR9Et --smearEleType=stochastic --regionsFile=data/regions/${region}.dat --EoP=${EoP_option} --invMass_var=invMass_SC_regrCorrSemiParV5_ele --autoBin --outDirImgData=${outDirData}/${extension}/\$LSB_JOBINDEX/img/ --smearerFit ${initFile} --commonCut=Et_20
 
	bsub -R "rusage[mem=10000]" -q 2nd\
            -oo ${outDirData}/${extension}/%I/fitres/`basename ${outFile} .dat`-${region}-stdout.log\
            -eo ${outDirData}/${extension}/%I/fitres/`basename ${outFile} .dat`-${region}-stderr.log\
            -J "${region} ${extension}[1-50]"\
            "cd $PWD; eval \`scramv1 runtime -sh\`; uname -a; echo \$CMSSW_VERSION;

./bin/ZFitter.exe -f data/validation/${validation_file} --corrEleType=HggRunEtaR9Et --smearEleType=stochastic --regionsFile=data/regions/${region}.dat --EoP=${EoP_option} --invMass_var=invMass_SC_regrCorrSemiParV5_ele --autoBin --outDirImgData=${outDirData}/${extension}/\$LSB_JOBINDEX/img/ --outDirFitResData=${outDirData}/${extension}/\$LSB_JOBINDEX/fitres --smearerFit ${initFile} --commonCut=Et_20 || exit 1";

exit 0
#At this point use EoP_Fitter.sh
#Qui aspetta che i job siano completati
	while [ "`bjobs -J \"${region} ${extension}\" | grep -v JOBID | grep -v found | wc -l`" != "0" ]; do /bin/sleep 2m; done

    ./script/haddTGraph.sh -o ${outDirData}/${extension}/fitres/outProfile-eop-$region-${commonCut}.root ${outDirData}/${extension}/*/fitres/outProfile-eop-$region-${commonCut}.root

    