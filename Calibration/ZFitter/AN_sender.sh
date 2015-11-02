#! /bin/bash

source="/afs/cern.ch/user/g/gfasanel/scratch1/www/Pt1Pt2/"

dest="/afs/cern.ch/user/g/gfasanel/CERN_documents/svnrepo/notes/AN-15-241/trunk/figs/energy_scale_factors/"


cp test/dato/ptRatio_pt2Sum/1/img/histos*.png ${dest}
cp test/dato/ptRatio_pt2Sum/1/fitres/*.png ${dest}
cp ${source}closure_test/EB/check*.png ${dest}
cp ${source}closure_test/EB/fit*.png ${dest}
cp ${source}closure_test/EE/check*.png ${dest}
cp ${source}closure_test/EE/fit*.png ${dest}
cp ${source}Run2/*.png ${dest}RUN2

