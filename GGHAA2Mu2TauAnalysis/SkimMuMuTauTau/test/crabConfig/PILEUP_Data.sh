#!/bin/bash
#parse arguments
DATE=`date +TIME_%H-%M-%S_DATE_%y-%m-%d`
#parse arguments
if [ $# -ne 1 ]; then
  echo "Going to CHECK all directories in crab_projects."
  DIR=""
  FILENAME="CRAB_CHECK_ALL__${DATE}.txt"
else
  DIR="${1}"
  echo "Going to check the directories with \"${DIR}\" in their name."
  FILENAME="CRAB_CHECK_ALL__${DIR}_${DATE}.txt"
fi


export SCRAM_ARCH=slc6_amd64_gcc481
cd /afs/cern.ch/user/k/ktos/SkimDir/CMSSW_8_0_17/src
eval `scramv1 runtime -sh`
cd -
source ~/SETENV.sh

DATE=`date +TIME_%H-%M-%S__DATE_%y-%m-%d`
FILENAME="CRAB_PILEUP_ALL_${new_tag}_${DATE}.txt"

for i in crab_projects/*${DIR}*; do 
  cd $i/results/; 
  pileupCalc.py -i processedLumis.json --inputLumiJSON  /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt  --calcMode true --minBiasXsec 62900 --maxPileupBin 80 --numPileupBins 80  pileupHistogram${i##*SEP4}.root >> /afs/cern.ch/user/k/ktos/SkimDir/CMSSW_8_0_17/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/crabConfig/$FILENAME
  cd ../../../; 
done

