#!/bin/bash

#export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/NewSkimDir/CMSSW_8_0_30/src
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/NewSkimDir/CMSSW_8_0_30/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/SkimSequence/BSUB/DIRNAME/ANALYZER.py .
cmsRun ANALYZER.py
rm ANALYZER.py 
exit 0
