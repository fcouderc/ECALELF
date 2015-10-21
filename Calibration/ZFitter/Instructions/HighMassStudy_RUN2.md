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
./bin/ZFitter.exe -f data/validation/run2_first.dat --regionsFile=data/regions/scaleStep0.dat --addBranch=smearerCat invMass_var=invMass_SC_corr --saveRootMacro
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
echo "s1      smearerCat_scaleStep0     friends/smearerCat/smearerCat_scaleStep0_s1-run2_first.root" >> data/validation/run2_first.dat
echo "d1      smearerCat_scaleStep0     friends/smearerCat/smearerCat_scaleStep0_d1-run2_first.root" >> data/validation/run2_first.dat
```

In generale se vuoi fare un plot dei tree associati ad un certo .dat file
```
./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --corrEleType=HggRunEtaR9Et --smearEleType=stochastic --saveR
ootMacro
./script/hadder.sh
```

Fai il load.C e plotta

### Qui sei pronto a riempire gli istogrammi con la target variable e minimizzare la likelihod
* Ti ricordo che in RooSmearer::SetSmearedHisto ho eliminato lo smoothing artificiale smearedHisto->Smooth() 

###Fai delle prove rapide prima di mandare i job
* Riempire solo gli istogrammi
```
./bin/ZFitter.exe -f data/validation/run2_first.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC_corr --autoBin --smearerFit --plotOnly --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*200 --targetVariable_binWidth=0.05*2 --configuration=random&> tmp/debug.txt #--corrEleType=HggRunEtaR9Et --smearEleType=stochastic 
root -l test/dato/fitres/histos_ptRatio_pt2Sum_random_scaleStep0_Et_25_trigger_noPF.root 
```

* Prova a fare un rapido scan della likelihood, senza la vera minimizzazione
```
./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC_regrCorrSemiParV5_ele --autoBin --smearerFit --plotOnly --profileOnly --targetVariable=ptRatio --targetVariable_min=0.5 --targetVariable_max=2 --targetVariable_binWidth=0.05 --corrEleType=HggRunEtaR9Et --smearEleType=stochastic &> tmp/debug.txt 
root -l test/dato/fitres/outProfile_scaleStep0_Et_25_trigger_noPF.root 
```
* Prova a vedere se riesce a trovare il minimo della likelihood da solo
```
./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC_regrCorrSemiParV5_ele --autoBin --smearerFit --targetVariable=ptRatio --targetVariable_min=0.5 --targetVariable_max=2 --targetVariable_binWidth=0.05 --configuration=leading --corrEleType=HggRunEtaR9Et --smearEleType=stochastic &> tmp/debug_ptRatio_leading.txt 
```
*Se non trova da solo il minimo, guarda il profilo e passagli un initFile con parametri vicini al minimo, cosi' da imbeccarlo
```
cp test/dato/fitres/params-scaleStep0-Et_25-noPF.txt init_pt12.txt 
e modifichi init_pt12.dat per avvicinare i parametri iniziali a quelli veri (allarghi i range se necessario) 
A questo punto rigiri con quell'initFile e cerchi di imboccarlo 
--initFile=
```
* Plottare 1 solo profile della likelihood
```
./script/fit.sh path_to_profile (senza .root)

//Questo e' ormai obsoleto
root -l -b
.L tmp/fitOneProfile.C
fitOneProfile("profile.root","outDir")
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



