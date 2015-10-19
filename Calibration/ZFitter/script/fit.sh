##Qui non fai haddTGraph per avere un .root file con tutti i vari profili che hai fatto
##Vuoi fittare un preciso .root file

echo "[STATUS] Fitting NLL"
echo "{" > tmp/my_fitProfile.C
echo "gROOT->ProcessLine(\".include $ROOFITSYS/include\");" >> tmp/my_fitProfile.C
echo "gROOT->ProcessLine(\".L macro/macro_fit.C+\");" >> tmp/my_fitProfile.C
echo "gROOT->ProcessLine(\".L macro/plot_data_mc.C+\");" >> tmp/my_fitProfile.C
echo "FitProfile2(\"${1}.root\",\"\",\"\",true,true,true);" >> tmp/my_fitProfile.C
echo "}" >> tmp/my_fitProfile.C
root -l -b -q tmp/my_fitProfile.C
