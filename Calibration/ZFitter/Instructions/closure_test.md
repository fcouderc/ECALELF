1) Avere un initFile con scala e sigma da iniettare ai Toy

2) Solo i file con label "s" vengono utilizzati nel closure test: quindi solo il MC

3) In SmearingImporter si decide la divisione in dati vs MC (in base a %5)

3) in src/RooSmearer.cc Init puoi decidere come inizializzare la scala/sigma di quello che hai chiamato MC

```
if(mcToy && false){

do something; puoi randomizzare ad esempio la scala e la sigma ecc...
}
```

da ZFitter in ECALELF_run2 (link simbolico)
quindi da qui:
new_version_ECALELF/CMSSW_7_5_0_pre4/src/Calibration/ZFitter/

./script/justOnce/closure_test.sh 2nd

Se non specifichi la coda (2nd) il tutto viene fatto in locale

Viene fatto questo:

./bin/ZFitter.exe --invMass_var=invMass_SC_regrCorrSemiParV5_ele -f data/validation/closure.dat --regionsFile=data/regions/scaleStep0.dat --smearerFit --outDirFitResData=test/dato/fitres/toys/scaleStep0/0.00_0.00_1.1/1/ --noPU  --plotOnly --profileOnly --runToy --eventsPerToy=0 --nSmearToy=1 --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*200 --targetVariable_binWidth=0.05*2 --configuration=random --initFile=test/dato/fitres/toys/scaleStep0/0.00_0.00_1.1/init_mcToy.txt --autoBin > debug.txt 

Per ogni punto di scala/smearing vengono sottomessi 10 job, organizzati in cartelle da 1 a 10.
test/dato/fitres/toys/scaleStep0/


./script/justOnce/closure_test_fitter.sh>debug_fit.txt

Nota che prima viene fatto il fit e poi il plot e che quindi il range del plot puo' NON coincidere con il range del fit


Come hai fatto a scoprire il doppio minimo


./bin/ZFitter.exe --invMass_var=invMass_SC_regrCorrSemiParV5_ele -f data/validation/closure.dat --regionsFile=data/regions/scaleStep0.dat --smearerFit --outDirFitResData=test/dato/fitres/toys/scaleStep0/0.02_0.00_0.97/9/ --noPU  --plotOnly --profileOnly --runToy --eventsPerToy=0 --nSmearToy=5 --targetVariable=ptRatio --targetVariable_min=0.5 --targetVariable_max=2 --targetVariable_binWidth=0.05 --configuration=random --initFile=test/dato/fitres/toys/scaleStep0/0.02_0.00_0.97/init_mcToy.txt --autoBin> debug.txt 

Per fare i plot cumulativi: ./script/justOnce/closure_test.py