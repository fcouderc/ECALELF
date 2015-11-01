```
bsub -q 1nh -eo err.dat -oo out.dat -J "noCorr job" "cd $PWD; eval \`scramv1 runtime -sh\`; uname -a; echo \$CMSSW_VERSION;
./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC_regrCorrSemiParV5_ele --autoBin --smearerFit --plotOnly --profileOnly --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*200 --targetVariable_binWidth=0.05*2 --configuration=random --outDirFitResData=tmp/ --initFile=init_RUN1_noCorr.txt"
```

```
./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC_regrCorrSemiParV5_ele --autoBin --smearerFit --plotOnly --profileOnly --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*200 --targetVariable_binWidth=0.05*2 --configuration=random --outDirFitResData=tmp/ --corrEleType=HggRunEtaR9Et --smearEleType=stochastic --initFile=init_RUN1_Corr.txt
```