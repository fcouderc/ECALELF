# Per creare una chain in tmp/ con dei particolari friend trees attaccati

Prendi il config file e modificalo inserendo i tree che vuoi, gli shift che vuoi, gli smearings che vuoi

./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --corrEleType=HggRunEtaR9Et --smearEleType=stochastic --saveRootMacro

./script/hadder.sh

* Nota che esplicitamente devi dichiarare solo smear e shift, gli altri friend trees vengono attaccati in maniera automatica.