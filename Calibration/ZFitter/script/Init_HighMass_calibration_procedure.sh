#! /bin/bash
source script/functions.sh
source script/bash_functions_calibration_highMass.sh

##file=$(echo $1 |sed "s|data/validation/||"| sed "s|.dat||" ) #so you can tab the name of the file :-)
file=`basename $1 .dat`
if [[ $1 = "" ]]; then
    echo "you should specify the validation file"; exit;
fi

##add branch r9 -->Not yet there for highMass
./script/addBranch.sh data/validation/${file}.dat R9Eleprime

##make pileupHist: they are used to make pileupTrees, and also they are required for step1
pileupHist
##make pileupTree
pileupTrees

##Load Z correction
invMass_type=$2
scale_file=/afs/cern.ch/user/g/gfasanel/scratch1/CMSSW_7_4_15/src/Calibration/ZFitter/test/dato/December2015_Rereco_C_D_withPho/loose/invMass_SC_pho_regrCorr/table/step2-invMass_SC_pho_regrCorr-loose-Et_20-noPF-HggRunEtaR9.dat  
smear_file=/afs/cern.ch/user/g/gfasanel/scratch1/CMSSW_7_4_15/src/Calibration/ZFitter/test/dato/December2015_Rereco_C_D_withPho/loose/invMass_SC_pho_regrCorr/table/outFile-step4-invMass_SC_pho_regrCorr-loose-Et_20-noPF-HggRunEtaR9-smearEle.dat

echo -e "\e[0;31m scale correction file is \e[0m ${scale_file}"
echo -e "\e[0;31m smear correction file is \e[0m ${smear_file}"
echo -e "\e[0;31m Change those file paths in script/Init_HighMass_calibration_procedure if you do not agree \e[0m "


./bin/ZFitter.exe -f data/validation/${file}.dat --regionsFile=data/regions/scaleStep0.dat invMass_var=${invMass_type} --saveRootMacro --corrEleType=Zcorr --corrEleFile=${scale_file} 

./bin/ZFitter.exe -f data/validation/${file}.dat --regionsFile=data/regions/scaleStep0.dat invMass_var=${invMass_type} --saveRootMacro --smearEleType=Zcorr --smearEleFile=${smear_file}

mv tmp/scaleEle_Zcorr_d*-${file}.root friends/others/
mv tmp/smearEle_Zcorr_s*-${file}.root friends/others/
#
for tag in `grep "^d" data/validation/${file}.dat | grep selected | awk -F" " ' { print $1 } '`
do  
    echo "${tag} scaleEle_Zcorr friends/others/scaleEle_Zcorr_${tag}-${file}.root" >> data/validation/${file}.dat 
done
for tag in `grep "^s" data/validation/${file}.dat | grep selected | awk -F" " ' { print $1 } '`
do
    echo "${tag} smearEle_Zcorr friends/others/smearEle_Zcorr_${tag}-${file}.root" >> data/validation/${file}.dat 
done