#! /bin/bash
Uncomment(){
./script/Init_HighMass_calibration_procedure.sh data/validation/February2016_76_Rereco_HighMass.dat invMass_SC_pho_regrCorr
./script/calibration_highMass.sh data/validation/February2016_76_Rereco_HighMass.dat invMass_SC_pho_regrCorr CatOnly
##Test_job
./script/calibration_highMass.sh data/validation/February2016_76_Rereco_HighMass.dat invMass_SC_pho_regrCorr Test_job | tee debug.txt
}
##If you need initParameters
#mv test/dato/img/outProfile_ptRatio_pt2Sum_random_scaleStep0_Et_35_noPF-FitResult-.config init_Parameters/highMass_RUN2.dat
#E modifichi il dat file con i valori di inizializzazione che piu' ti soddisfano

####Sottometti 50 job
./script/calibration_highMass.sh data/validation/February2016_76_Rereco_HighMass.dat invMass_SC_pho_regrCorr jobs --initFile=init_Parameters/highMass_RUN2.dat 

 ####Likelihood plots, data_MC e dat file
./script/calibration_highMass.sh data/validation/February2016_76_Rereco_HighMass.dat invMass_SC_pho_regrCorr finalize_plots