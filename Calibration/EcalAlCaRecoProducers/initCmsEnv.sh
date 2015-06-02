#!/bin/bash
eval `scramv1 runtime -sh`
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms -hours 192 -valid 192:00
#voms-proxy-init -voms cms -out $HOME/gpi.out
PATH=$PATH:/afs/cern.ch/project/eos/installation/pro/bin/
PATH=$PATH:$CMSSW_BASE/src/Calibration/ALCARAW_RECO/bin
export ECALELFINIT=y