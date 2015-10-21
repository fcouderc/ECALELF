### Start from a config file and do your trees in tmp/

./bin/ZFitter.exe -f data/validation/22Jan2012-runDepMCAll_checkMee.dat --regionsFile=data/regions/scaleStep0.dat --smearEleType=stochastic --corrEleType=HggRunEtaR9Et --saveRootMacro &> debug_2.txt
./script/hadder.sh

### Prototipo per 1D variable
root -l -b -q tmp/d_chain.root tmp/s1_chain.root tmp/load_dataMC.C tmp/massPlotter.C

### Prototipo per 2D variable
root -l -b -q tmp/d_chain.root tmp/s1_chain.root tmp/load_dataMC.C tmp/2DPlotter.C

### in load_dataMC.C viene caricato il file di stile

Plots will be in tmp/Plots

root -l tmp/d_chain.root tmp/s1_chain.root tmp/load_dataMC.C

data->Draw("invMass_SC")

data->Draw("invMass_SC>>h","invMass_SC>200")