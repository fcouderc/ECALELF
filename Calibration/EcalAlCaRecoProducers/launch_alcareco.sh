source initCmsEnv.sh
./scripts/prodAlcareco.sh `parseDatasetFile.sh alcareco_datasets.dat | grep 26May_evening` --skim ZSkim --type ALCARECO --isMC --jobname DY200_400 --createOnly
