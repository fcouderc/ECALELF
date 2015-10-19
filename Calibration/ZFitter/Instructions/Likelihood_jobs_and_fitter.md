./script/run50Minizations.sh ptRatio*pt2Sum

```
./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --invMass_var=invMass_SC_regrCorrSemiParV5_ele --autoBin --smearerFit --targetVariable=ptRatio*pt2Sum --targetVariable_min=0.5*0 --targetVariable_max=2*200 --targetVariable_binWidth=0.05*2 --configuration=random --corrEleType=HggRunEtaR9Et --smearEleType=stochastic --outDirImgData=test/dato/ptRatio_pt2Sum/$LSB_JOBINDEX/img/ --outDirFitResData=test/dato/ptRatio_pt2Sum/$LSB_JOBINDEX/fitres/
```

Quando i job sono finiti

./script/Likelihood_fitter.sh ptRatio*pt2Sum

Cosi' fitti la media, ma puoi anche fittarne uno solo