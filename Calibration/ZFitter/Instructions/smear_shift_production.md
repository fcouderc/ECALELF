# Start from data/validation/22Jan2012-runDepMCAll_checkMee.dat, with no scale or smearing root file specified

./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --smearEleType=stochastic --smearEleFile=mc_smear/smearing_sigma_and_errors_stocastic_rd_mc.dat --corrEleType=HggRunEtaR9Et --corrEleFile=data_scale/step8-invMass_SC_regrCorrSemiParV5_pho-loose-Et_20-trigger-noPF-HggRunEtaR9Et.dat --saveRootMacro &> debug.txt
./script/hadder.sh

#check shifts and smearings, plotting the corrected invariantMass
root -l -b -q tmp/d_chain.root tmp/s1_chain.root tmp/load_dataMC.C tmp/massPlotter.C

Plots will be in tmp/Plots

* Il parsing dei file di smearing viene fatto in src/EnergyScaleCorrection_class.cc, da ReadSmearingFile
* Esistono due formati di file (mi pare): formato 1 (globe) e formato 2 (ecalelf)
* C'era un errore nel parsing e il file veniva identificato come formato 2, mentre invece era formato 1
* Ho commentato //err_phi nel metodo ReadSmearingFile e ha funzionato

* Se sei soddisfatto e il plot di massa invariante ti torna
# sposta i root file in friends/other
mv tmp/scaleEle_HggRunEtaR9Et_d1-22Jan2012-runDepMCAll_checkMee.root friends/other
mv tmp/smearEle_stochastic_s1-22Jan2012-runDepMCAll_checkMee.root friends/other

# Inserisci i root file prodotti con config file dal quale sei partito