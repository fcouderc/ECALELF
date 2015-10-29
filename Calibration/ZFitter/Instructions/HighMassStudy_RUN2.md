* RUN1 file: data/validation/22Jan2012-runDepMCAll_checkMee.dat

* per lanciare i boost ho inserito in .bashrc 
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/cms/slc5_amd64_gcc434/external/boost/1.47.0/lib
```

### Da ZFitter: 
```
cd new_version_ECALELF/CMSSW_7_5_0_pre4/src/Calibration/ZFitter/
```
* Prima di incominciare proprio, chiarisci a te stesso quale variabile di energia vuoi utilizzare 

```
invMass_var=invMass_SC_regrCorrSemiParV5_ele 
```

* L'energia giusta va specificata A MANO in src/ElectronCategory_class.cc 

(FALLO, altrimenti categorizzi con l'energia sbagliata) 

* A mano devi anche dire in src/ElectronCategory_class.cc se smearEle(true) o smearEle(false)

* Aggiungere il branch ZPt (mi serve perche' ci categorizzo sopra) 
```
./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --addBranch=ZPt_energySCEle_regrCorrSemiParV5_ele --corrEleType=HggRunEtaR9Et --smearEleType=stochastic --saveRootMacro &> debug.txt
./script/hadder.sh
mv tmp/ZPt_energySCEle_* friends/other/
``` 
* Aggiungi i branch di ZPt nel config file 

*Categorizza* 

* Ti ricordo che l'energia giusta va specificata A MANO in src/ElectronCategory_class.cc
* RICOMPILA e LANCIA

./script/GenRootChain.sh non accetta gli smearing e allora:

```
./bin/ZFitter.exe -f data/validation/reference_25nsReco.dat --regionsFile=data/regions/scaleStep0.dat --addBranch=smearerCat invMass_var=invMass_SC_corr --saveRootMacro
./script/hadder.sh
```

* Controllare che la categorizzazione sia fatta bene: 

Devi associare i friend trees con tmp/load.C. A questo punto puoi fare plot in questo modo, chiedendo branch da tree diversi 

```
root -l tmp/d_chain.root tmp/load_singleFile.C 
#data->Draw("ZPta_energySCEle_regrCorrSemiParV5_ele","smearerCat[0]>0") 
data->Draw("invMass_SC_corr","smearerCat[0]>0") 

root -l tmp/s_chain.root tmp/load_singleFile.C
data->Draw("invMass_SC_corr","smearerCat[0]>0")
```
* Se ti soddisfa la categorizzazione, aggiungi i branch di categoria nel config file 

```	
mv tmp/smearerCat_scaleStep0* friends/smearerCat/ 
echo "s1      smearerCat_scaleStep0     friends/smearerCat/smearerCat_scaleStep0_s1-reference_25nsReco.root" >> data/validation/reference_25nsReco.dat
echo "d1      smearerCat_scaleStep0     friends/smearerCat/smearerCat_scaleStep0_d1-reference_25nsReco.root" >> data/validation/reference_25nsReco.dat
echo "d2      smearerCat_scaleStep0     friends/smearerCat/smearerCat_scaleStep0_d2-reference_25nsReco.root" >> data/validation/reference_25nsReco.dat
echo "d3      smearerCat_scaleStep0     friends/smearerCat/smearerCat_scaleStep0_d3-reference_25nsReco.root" >> data/validation/reference_25nsReco.dat
```

In generale se vuoi fare un plot dei tree associati ad un certo .dat file, ad esempio reference_25nsReco.dat
```
#Devi sapere se vuoi le correzioni o no
#--corrEleType=HggRunEtaR9Et --smearEleType=stochastic
./bin/ZFitter.exe -f data/validation/reference_25nsReco.dat --regionsFile=data/regions/scaleStep0.dat  --saveRootMacro
./script/hadder.sh
#Fai ora il load dei friend trees e plotta
root -l tmp/d_chain.root tmp/load_singleFile.C 
TString commonCut="(((((eleID[0] & 2)==2)&&((eleID[1] & 2)==2))&&((energySCEle_corr[0]/cosh(etaSCEle[0]) >= 25)&&(energySCEle_corr[1]/cosh(etaSCEle[1]) >= 25)))&&(HLTfire==1))&&(recoFlagsEle[0] > 1 && recoFlagsEle[1] > 1)";
TString commonCut="((((eleID[0] & 2)==2)&&((eleID[1] & 2)==2))&&((energySCEle_corr[0]/cosh(etaSCEle[0]) >= 25)&&(energySCEle_corr[1]/cosh(etaSCEle[1]) >= 25)))&&(recoFlagsEle[0] > 1 && recoFlagsEle[1] > 1)";
data->Draw("invMass_SC_corr","(smearerCat[0]>0)&&"+commonCut);
data->Draw("HLTfire");
```

### Qui sei pronto a riempire gli istogrammi con la target variable e minimizzare la likelihod
* Ti ricordo che in RooSmearer::SetSmearedHisto ho eliminato lo smoothing artificiale smearedHisto->Smooth() 

###Fai delle prove rapide prima di mandare i job
* Riempire solo gli istogrammi
```
#--corrEleType=HggRunEtaR9Et --smearEleType=stochastic 
./bin/ZFitter.exe -f data/validation/reference_25nsReco.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC_corr --commonCut=Et_25-noPF --autoBin --smearerFit --plotOnly --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*200 --targetVariable_binWidth=0.05*2 --configuration=random &> tmp/debug.txt  
root -l test/dato/fitres/histos_ptRatio_pt2Sum_random_scaleStep0_Et_25_noPF.root 
```

* Prova a fare un *rapido scan della likelihood*, senza la vera minimizzazione
```
#--corrEleType=HggRunEtaR9Et --smearEleType=stochastic
#--commonCut=Et_25-noPF
./bin/ZFitter.exe -f data/validation/reference_25nsReco.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC_corr --commonCut=Et_25-noPF --autoBin --smearerFit --plotOnly --profileOnly --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*200 --targetVariable_binWidth=0.05*2 --configuration=random &> tmp/debug.txt  
root -l test/dato/fitres/outProfile_scaleStep0_Et_25_noPF.root 
```
* Prova a vedere se riesce a trovare il *minimo della likelihood* da solo, senza dargli nessun initFile
```
./bin/ZFitter.exe -f data/validation/reference_25nsReco.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC_corr --commonCut=Et_25-noPF --autoBin --smearerFit --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*200 --targetVariable_binWidth=0.05*2 --configuration=random &> tmp/debug.txt  
```

* Plottare 1 solo profile della likelihood
```
./script/fit.sh path_to_profile (senza .root)

//Questo e' ormai obsoleto
root -l -b
.L tmp/fitOneProfile.C
fitOneProfile("profile.root","outDir")
```

*Se non trova da solo il minimo, guarda il profilo e passagli un initFile con parametri vicini al minimo, cosi' da imbeccarlo
```
cp test/dato/fitres/params-ptRatio_pt2Sum_random_scaleStep0-Et_25-noPF.txt init_RUN2.txt 
e modifichi init_RUN2.txt per avvicinare i parametri iniziali a quelli veri (allarghi i range se necessario) 
A questo punto rigiri con quell'initFile e cerchi di imboccarlo 
./bin/ZFitter.exe -f data/validation/reference_25nsReco.dat --regionsFile=data/regions/scaleStep0.dat --initFile=init_RUN2.txt --invMass_var=invMass_SC_corr --commonCut=Et_25-noPF --autoBin --smearerFit --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*200 --targetVariable_binWidth=0.05*2 --configuration=random &> tmp/debug.txt  
```
* Plotta il profile della likelihood
```
./script/fit.sh test/dato/fitres/outProfile_ptRatio_pt2Sum_random_scaleStep0_Et_25_noPF
source ./script/organizer_plots.sh

```

* Fai il plot Data/MC smearato per vedere se la minimizzazione della likelihood ha senso
```
cp test/dato/fitres/params-ptRatio_pt2Sum_random_scaleStep0-Et_25-noPF.txt output.txt
rm test/dato/fitres/histos_ptRatio_pt2Sum_random_scaleStep0_Et_25_noPF.root 
# --profileOnly
./bin/ZFitter.exe -f data/validation/reference_25nsReco.dat --regionsFile=data/regions/scaleStep0.dat --initFile=output.txt --invMass_var=invMass_SC_corr --commonCut=Et_25-noPF --autoBin --smearerFit --plotOnly --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*200 --targetVariable_binWidth=0.05*2 --configuration=random &> tmp/debug.txt  
root -l -b 
.L macro/plot_data_mc.C+ 
PlotMeanHist("test/dato/fitres/histos_ptRatio_pt2Sum_random_scaleStep0_Et_25_noPF.root") 
.q 
source ./script/organizer_plots.sh
```
### Minimizzare la likelihood (manda 50 job e fitta la likelihood media)
```
./script/run50Minimization.sh ptRatio (opp ptRatio_ranom opp ptSum)
./script/Likelihodd_fitter.sh ptRatio
```

### Plot Data/MC finali per vedere se la minimizzazione ha senso
```
cp test/dato/fitres/params-scaleStep0-Et_25-trigger-noPF.txt output_pt12.txt 
rm test/dato/fitres/histos_ptRatio_leading_scaleStep0_Et_25_trigger_noPF.root 
--plotOnly --profileOnly
root -l -b 
.L macro/plot_data_mc.C+ 
PlotMeanHist("test/dato/fitres/histos_ptRatio_leading_scaleStep0_Et_25_trigger_noPF.root") 
.q 
```



