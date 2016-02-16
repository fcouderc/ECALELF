#! /bin/bash
#1)Fai pileupHist
commonCut=Et_35-noPF

Categorize_region(){
./bin/ZFitter.exe -f data/validation/${file}.dat --regionsFile=data/regions/${region1}.dat --addBranch=smearerCat invMass_var=${invMass_type} --saveRootMacro

#categorie
mv tmp/smearerCat_${region1}_*${file}* friends/smearerCat/ 
for tag in `grep "^s" data/validation/${file}.dat | grep selected | awk -F" " ' { print $1 } '`
do
    echo "${tag} smearerCat_${region} friends/smearerCat/smearerCat_${region}_${tag}-${file}.root" >> data/validation/${file}.dat 
done
for tag in `grep "^d" data/validation/${file}.dat | grep selected | awk -F" " ' { print $1 } '`
do  
    echo "${tag} smearerCat_${region} friends/smearerCat/smearerCat_${region}_${tag}-${file}.root" >> data/validation/${file}.dat 
done
}


Test_1_job_highMass(){
echo "[INFO] you are correcting scale with Zcorr"
echo "[INFO] you are correcting smearing with Zcorr"
echo "[INFO] Using initParameters: ${initParameters1}"
./bin/ZFitter.exe -f data/validation/${file}.dat --regionsFile=data/regions/${region1}.dat --invMass_var=${invMass_type} --commonCut=${commonCut} --autoBin --smearerFit --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*300 --targetVariable_binWidth=0.02*6 --configuration=random --corrEleType=Zcorr --smearEleType=Zcorr ${initParameters1} #--plotOnly --profileOnly
#./script/fit.sh test/dato/fitres/outProfile-${region1}-${commonCut}.root
#./script/fit.sh test/dato/fitres/outProfile-${region2}-${commonCut}.root
#Likelihood_plot_dir=~/scratch1/www/Validation_ntuple_Paolo
#mv test/dato/img/outProfile-${region1}-${commonCut}-*.png ${Likelihood_plot_dir}
#mv test/dato/img/outProfile-${region2}-${commonCut}-*.png ${Likelihood_plot_dir}
}


Submit_50_jobs_highMass(){
echo "[INFo} you are submitting jobs: the output will be in ${outDirData}/${extension}/fitres"
#####Una volta tunati i job di prova, sottometti 50 job
 if [ -e "${outDirData}/${extension}/fitres" ];then echo "${outDirData}/${extension}/fitres already exists: be sure not to mix different job tasks---> EXIT"; exit;fi
 if [ ! -e "${outDirData}/${extension}/fitres" ];then mkdir ${outDirData}/${extension}/fitres -p; fi
 if [ ! -e "${outDirData}/${extension}/img" ];then mkdir ${outDirData}/${extension}/img -p; fi

############# SUBMITTING JOBS ###########################################

 for index in `seq 1 50` #Making directories
 do
     mkdir ${outDirData}/${extension}/${index}/fitres/ -p
     mkdir ${outDirData}/${extension}/${index}/img -p
 done
 
 bsub -q 2nd\
            -oo ${outDirData}/${extension}/%I/fitres/${region1}-stdout.log\
            -eo ${outDirData}/${extension}/%I/fitres/${region1}-stderr.log\
            -J "${region1} ${extension}[1-50]"\
            "cd $PWD; eval \`scramv1 runtime -sh\`; uname -a; echo \$CMSSW_VERSION;
./bin/ZFitter.exe -f data/validation/${file}.dat --regionsFile=data/regions/${region1}.dat --invMass_var=${invMass_type}  --commonCut=${commonCut} --autoBin --smearerFit --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*300 --targetVariable_binWidth=0.02*6 --configuration=random --corrEleType=Zcorr --smearEleType=Zcorr ${initParameters1} --outDirFitResData=${outDirData}/${extension}/\$LSB_JOBINDEX/fitres/;"  
}


Likelihood_Fitter_highMass(){
#check jobs
while [ "`bjobs -J \"${region1} ${extension}[1-50]\" | grep -v JOBID | grep -v found | wc -l`" != "0" ]; do /bin/sleep 2m; done
###########################Make Likelihood plots###########################
./script/haddTGraph.sh -o ${outDirData}/${extension}/fitres/outProfile_ptRatio_pt2Sum_random_${region1}_Et_35_noPF.root ${outDirData}/${extension}/*/fitres/outProfile_ptRatio_pt2Sum_random_${region1}_Et_35_noPF.root
echo "{" > tmp/fitProfiles_${invMass_type}_${extension}.C
echo "gROOT->ProcessLine(\".include $ROOFITSYS/include\");" >> tmp/fitProfiles_${invMass_type}_${extension}.C
echo "gROOT->ProcessLine(\".L macro/macro_fit.C+\");" >> tmp/fitProfiles_${invMass_type}_${extension}.C
echo "gROOT->ProcessLine(\".L macro/plot_data_mc.C+\");" >> tmp/fitProfiles_${invMass_type}_${extension}.C
echo "FitProfile2(\"${outDirData}/${extension}/fitres/outProfile_ptRatio_pt2Sum_random_${region1}_Et_35_noPF.root\",\"\",\"\",true,true,true);" >> tmp/fitProfiles_${invMass_type}_${extension}.C
echo "}" >> tmp/fitProfiles_${invMass_type}_${extension}.C
if [ ! -e "${outDirData}/${extension}/img" ];then mkdir test/dato/img/${extension}/img -p; fi
root -l -b -q tmp/fitProfiles_${invMass_type}_${extension}.C
if [ ! -e "~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/" ];then mkdir ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/ -p; cp ~/scratch1/www/index.php ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/; fi
mv ${outDirData}/${extension}/img/*.png ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/
mv ${outDirData}/${extension}/img/*.C ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/
echo "************Find png plots in ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/"
}

Plot_data_MC_highMass(){
##In teoria dovresti farti --plotOnly e poi...
if [ ! -e "~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/data_MC/" ]; then mkdir ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/data_MC -p; cp ~/scratch1/www/index.php ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/data_MC/; fi
echo "{" > tmp/plotter_data_MC.C
echo "gROOT->ProcessLine(\".L macro/plot_data_mc.C+\");" >> tmp/plotter_data_MC.C
echo "PlotMeanHist(\"${outDirData}/${extension}/1/fitres/histos_ptRatio_pt2Sum_random_${region1}_Et_35_noPF.root\");" >> tmp/plotter_data_MC.C
echo "}" >> tmp/plotter_data_MC.C
root -l -b -q ~/rootlogon.C tmp/plotter_data_MC.C
mv ${outDirData}/${extension}/1/./img/histos_ptRatio_pt2Sum_random_${region1}* ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/data_MC/
}

Write_down_dat_corr_highMass(){
##Scale corrections: per i dati -> 1/scala_MC
grep scale ${outDirData}/${extension}/img/outProfile_ptRatio_pt2Sum_random_${region1}_Et_35_noPF-FitResult-.config |  sed -r 's|[ ]+|\t|g;' | cut -f 1,3,5 | sed "s|-${commonCut}||g"> tmp/corr_MC.dat

categories=`grep scale tmp/corr_MC.dat | cut -f 1`
categories=`echo $categories | sed "s/scale_//g"`
categories=(${categories// / }) # array
scales=`grep scale tmp/corr_MC.dat | cut -f 2`
scales=(${scales// / }) # array
scales_data=()
errors=`grep scale tmp/corr_MC.dat | cut -f 3`
errors=(${errors// / }) # array

for scale in "${scales[@]}"
do
scales_data+=(`echo "1./$scale" | bc -l`)
done

if [ -e ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/scale_corrections_${file}_${extension}.dat ]; then rm ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/scale_corrections_${file}_${extension}.dat; fi
for i in "${!scales_data[@]}"; do 
echo ${categories[$i]} ${scales_data[$i]} ${errors[$i]} >> ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/scale_corrections_${extension}.dat
done

grep constTerm ${outDirData}/${extension}/img/outProfile_ptRatio_pt2Sum_random_${region1}_Et_35_noPF-FitResult-.config |  sed -r 's|[ ]+|\t|g;' | cut -f 1,3,5 | sed "s|-${commonCut}||g" > ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/smearing_corrections_${extension}.dat
}

Write_down_root_corr(){
Write_down_dat_corr
#qui, il file di regione non serve a niente, solo a far funzionare lo script che senza non ti fa fare niente
./bin/ZFitter.exe -f data/validation/${file}.dat --regionsFile=data/regions/${region1}.dat invMass_var=${invMass_type} --saveRootMacro --corrEleType=EtaR9_${extension} --corrEleFile=data_scale/scale_corrections_${file}_${extension}.dat 
./bin/ZFitter.exe -f data/validation/${file}.dat --regionsFile=data/regions/${region1}.dat invMass_var=${invMass_type} --saveRootMacro --smearEleType=stochastic_${extension} --smearEleFile=mc_smear/smearing_corrections_${file}_${extension}.dat

mv tmp/scaleEle_EtaR9_${extension}_d*-${file}.root friends/others/
mv tmp/smearEle_stochastic_${extension}_s*-${file}.root friends/others/

for tag in `grep "^s" data/validation/${file}.dat | grep selected | awk -F" " ' { print $1 } '`
do
    echo "${tag} smearEle_stochastic_${extension} friends/others/smearEle_stochastic_${extension}_${tag}-${file}.root" >> data/validation/${file}.dat 
done
for tag in `grep "^d" data/validation/${file}.dat | grep selected | awk -F" " ' { print $1 } '`
do  
    echo "${tag} scaleEle_EtaR9_${extension} friends/others/scaleEle_EtaR9_${extension}_${tag}-${file}.root" >> data/validation/${file}.dat 
done
}

closure_test(){
#senza initParameters, ma applicando le correzioni trovate
#Likelihood con minimo a 1 per la scala e zero per lo smearing
#Data/MC plot perfetti
./bin/ZFitter.exe -f data/validation/${file}.dat --regionsFile=data/regions/${region1}.dat --invMass_var=${invMass_type}  --commonCut=${commonCut} --autoBin --smearerFit --corrEleType=EtaR9_${extension} --smearEleType=stochastic_${extension} --profileOnly --plotOnly > test/dato/fitres/closure_step2_1.dat
./bin/ZFitter.exe -f data/validation/${file}.dat --regionsFile=data/regions/${region2}.dat --invMass_var=${invMass_type}  --commonCut=${commonCut} --autoBin --smearerFit --corrEleType=EtaR9_${extension} --smearEleType=stochastic_${extension} --profileOnly --plotOnly > test/dato/fitres/closure_step2_2.dat
if [ ! -e "test/dato/${extension}/closure_step2_${file}_${extension}/" ]; then mkdir test/dato/${extension}/closure_step2_${file}_${extension}/; fi
mv test/dato/fitres/closure_step2*.dat test/dato/${extension}/closure_step2_${file}_${extension}/ #output of the closure tests "job"
mv test/dato/fitres/outProfile-${region1}-${commonCut}.root test/dato/${extension}/closure_step2_${file}_${extension}/
mv test/dato/fitres/outProfile-${region2}-${commonCut}.root test/dato/${extension}/closure_step2_${file}_${extension}/
mv test/dato/fitres/histos-${region1}-${commonCut}.root test/dato/${extension}/closure_step2_${file}_${extension}/
mv test/dato/fitres/histos-${region2}-${commonCut}.root test/dato/${extension}/closure_step2_${file}_${extension}/
}

Prepare_Table(){
Write_down_dat_corr
echo "\begin{table}[htb]"> tmp/table_${file}_${extension}.tex
echo " \caption{Results for scale and smearing corrections: scale has to be applied to data, smearing to MC.}">> tmp/table_${file}_${extension}.tex
echo " \begin{center}">> tmp/table_${file}_${extension}.tex
echo "   \begin{tabular}{ccc}">> tmp/table_${file}_${extension}.tex
echo "     \hline">> tmp/table_${file}_${extension}.tex
echo "     Category & Scale & Smearing\\\\ \hline">> tmp/table_${file}_${extension}.tex #I want to write \\ -> you have to escape \ 

for i in "${!scales_data[@]}"; do
categories_name[$i]=`echo ${categories[$i]} | sed "s/_/-/g"`
printf "%s %s %.5lf %s %.5lf %s %.4lf %s %.4lf %s\n" "${categories_name[$i]}" "&" "${scales_data[$i]}" "$\pm$" "${errors[$i]}" "&" "${smearings[$i]}" "$\pm$" "${errors_smear[$i]}" "\\\\">> tmp/table_${file}_${extension}.tex
done


echo "   \hline">> tmp/table_${file}_${extension}.tex
echo "    \end{tabular}">> tmp/table_${file}_${extension}.tex
echo "  \end{center}">> tmp/table_${file}_${extension}.tex
echo "\end{table}">> tmp/table_${file}_${extension}.tex

mv tmp/table_${file}_${extension}.tex ${outDirData}/${extension}/fitres/table_${file}_${extension}.tex
}

Prepare_Slide_8_cat_profileOnly(){
Prepare_Table
echo "\documentclass[9pt, xcolor=dvipsnames]{beamer}" >tmp/slide_${file}_${extension}_profileOnly.tex
echo "\usepackage{hyperref}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{document}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex
less ${outDirData}/${extension}/fitres/table_${file}_${extension}.tex >>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex

####SCALE step2_1
echo "\begin{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region1}-${commonCut}-scale${categories[0]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region1}-${commonCut}-scale${categories[1]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\caption{scale parameter for ${categories_name[0]} (left) and ${categories_name[1]} (right)}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{figure}">>tmp/slide_${file}_${extension}_profileOnly.tex

echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region1}-${commonCut}-scale${categories[2]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region1}-${commonCut}-scale${categories[3]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\caption{scale parameter for ${categories_name[2]} (left) and ${categories_name[3]} (right)}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{figure}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex

####SCALE step2_2
echo "\begin{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region2}-${commonCut}-scale${categories[4]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region2}-${commonCut}-scale${categories[5]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\caption{scale parameter for ${categories_name[4]} (left) and ${categories_name[5]} (right)}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{figure}">>tmp/slide_${file}_${extension}_profileOnly.tex

echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region2}-${commonCut}-scale${categories[6]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region2}-${commonCut}-scale${categories[7]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\caption{scale parameter for ${categories_name[6]} (left) and ${categories_name[7]} (right)}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{figure}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex

####constTerm step2_1
echo "\begin{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region1}-${commonCut}-constTerm${categories[0]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region1}-${commonCut}-constTerm${categories[1]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\caption{constTerm parameter for ${categories_name[0]} (left) and ${categories_name[1]} (right)}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{figure}">>tmp/slide_${file}_${extension}_profileOnly.tex

echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region1}-${commonCut}-constTerm${categories[2]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region1}-${commonCut}-constTerm${categories[3]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\caption{constTerm parameter for ${categories_name[2]} (left) and ${categories_name[3]} (right)}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{figure}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex

####constTerm step2_2
echo "\begin{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region2}-${commonCut}-constTerm${categories[4]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region2}-${commonCut}-constTerm${categories[5]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\caption{constTerm parameter for ${categories_name[4]} (left) and ${categories_name[5]} (right)}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{figure}">>tmp/slide_${file}_${extension}_profileOnly.tex

echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\begin{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region2}-${commonCut}-constTerm${categories[6]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region2}-${commonCut}-constTerm${categories[7]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{center}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\caption{constTerm parameter for ${categories_name[6]} (left) and ${categories_name[7]} (right)}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{figure}">>tmp/slide_${file}_${extension}_profileOnly.tex
echo "\end{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex

echo "\begin{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex
##Sostituire gli _ con \_ per latex
path=$(echo $extension | sed "s|_|\\\_|g")
echo "Find data/MC plots in \href{https://gfasanel.web.cern.ch/gfasanel/${path}/data\_MC/}{https://gfasanel.web.cern.ch/gfasanel/${path}/data\_MC/}">>tmp/slide_${file}_${extension}_profileOnly.tex
#
echo "\end{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex

echo "\end{document}">>tmp/slide_${file}_${extension}_profileOnly.tex

cd tmp/
pdflatex slide_${file}_${extension}_profileOnly.tex
cd ../

mv tmp/slide_${file}_${extension}_profileOnly.tex ${outDirData}/${extension}/
mv tmp/slide_${file}_${extension}_profileOnly.pdf ${outDirData}/${extension}/
cp ${outDirData}/${extension}/slide_${file}_${extension}_profileOnly.tex ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/
cp ${outDirData}/${extension}/slide_${file}_${extension}_profileOnly.pdf ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/

}

Prepare_Slide_8_cat(){
Prepare_Table
echo "\documentclass[9pt, xcolor=dvipsnames]{beamer}" >tmp/slide_${file}_${extension}.tex
echo "\usepackage{hyperref}">>tmp/slide_${file}_${extension}.tex
echo "\begin{document}">>tmp/slide_${file}_${extension}.tex
echo "\begin{frame}">>tmp/slide_${file}_${extension}.tex
less ${outDirData}/${extension}/fitres/table_${file}_${extension}.tex >>tmp/slide_${file}_${extension}.tex
echo "\end{frame}">>tmp/slide_${file}_${extension}.tex

for index in `seq 0 3`
do
    echo "\begin{frame}">>tmp/slide_${file}_${extension}.tex
    echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}.tex
    echo "\begin{center}">>tmp/slide_${file}_${extension}.tex
    echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region1}-${commonCut}-scale${categories[$index]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}.tex
    echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region1}-${commonCut}-constTerm${categories[$index]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}.tex
    echo "\end{center}">>tmp/slide_${file}_${extension}.tex
    echo "\caption{scale parameter (left) and smearing parameter (right) for ${categories_name[$index]}}">>tmp/slide_${file}_${extension}.tex
    echo "\end{figure}">>tmp/slide_${file}_${extension}.tex
    echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}.tex
    echo "\begin{center}">>tmp/slide_${file}_${extension}.tex
    echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/data_MC/histos-${region1}-${commonCut}_${categories[$index]}-${commonCut}${categories[$index]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}.tex
    echo "\end{center}">>tmp/slide_${file}_${extension}.tex
    echo "\caption{Data/MC comparison for ${categories_name[$index]}}">>tmp/slide_${file}_${extension}.tex
    echo "\end{figure}">>tmp/slide_${file}_${extension}.tex
    echo "\end{frame}">>tmp/slide_${file}_${extension}.tex
done

for index in `seq 4 7`
do
    echo "\begin{frame}">>tmp/slide_${file}_${extension}.tex
    echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}.tex
    echo "\begin{center}">>tmp/slide_${file}_${extension}.tex
    echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region2}-${commonCut}-scale${categories[$index]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}.tex
    echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/outProfile-${region2}-${commonCut}-constTerm${categories[$index]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}.tex
    echo "\end{center}">>tmp/slide_${file}_${extension}.tex
    echo "\caption{scale parameter (left) and smearing parameter (right) for ${categories_name[$index]}}">>tmp/slide_${file}_${extension}.tex
    echo "\end{figure}">>tmp/slide_${file}_${extension}.tex
    echo "\begin{figure}[!htbp]">>tmp/slide_${file}_${extension}.tex
    echo "\begin{center}">>tmp/slide_${file}_${extension}.tex
    echo "\includegraphics[width=0.43\textwidth]{{/afs/cern.ch/user/g/gfasanel/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/data_MC/histos-${region2}-${commonCut}_${categories[$index]}-${commonCut}${categories[$index]}-${commonCut}}.png}">>tmp/slide_${file}_${extension}.tex
    echo "\end{center}">>tmp/slide_${file}_${extension}.tex
    echo "\caption{Data/MC comparison for ${categories_name[$index]}}">>tmp/slide_${file}_${extension}.tex
    echo "\end{figure}">>tmp/slide_${file}_${extension}.tex
    echo "\end{frame}">>tmp/slide_${file}_${extension}.tex
done

#echo "\begin{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex
###Sostituire gli _ con \_ per latex
#path=$(echo $extension | sed "s|_|\\\_|g")
#echo "Find data/MC plots in \href{https://gfasanel.web.cern.ch/gfasanel/${path}/data\_MC/}{https://gfasanel.web.cern.ch/gfasanel/${path}/data\_MC/}">>tmp/slide_${file}_${extension}.tex
##
#echo "\end{frame}">>tmp/slide_${file}_${extension}_profileOnly.tex
#
echo "\end{document}">>tmp/slide_${file}_${extension}.tex

cd tmp/
pdflatex slide_${file}_${extension}.tex
cd ../

mv tmp/slide_${file}_${extension}.tex ${outDirData}/${extension}/
mv tmp/slide_${file}_${extension}.pdf ${outDirData}/${extension}/
cp ${outDirData}/${extension}/slide_${file}_${extension}.tex ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/
cp ${outDirData}/${extension}/slide_${file}_${extension}.pdf ~/scratch1/www/RUN2_ECAL_Calibration/HighMass/${extension}/

}