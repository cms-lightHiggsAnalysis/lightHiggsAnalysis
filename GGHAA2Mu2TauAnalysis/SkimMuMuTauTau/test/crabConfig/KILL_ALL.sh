#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Going to CHECK all directories in crab_projects."
  DIR=""
  FILENAME="CRAB_CHECK_ALL__${DATE}.txt"
else
  DIR="${1}"
  echo "Going to check the directories with \"${DIR}\" in their name."
  FILENAME="CRAB_CHECK_ALL__${DIR}_${DATE}.txt"
fi

#export SCRAM_ARCH=slc6_amd64_gcc481
#cd /afs/cern.ch/user/k/ktos/SkimDir/CMSSW_8_0_17/src
eval `scramv1 runtime -sh`
#cd -

DATE=`date +TIME_%H-%M-%S_DATE_%y-%m-%d`
FILENAME="CRAB_KILL_ALL_${DATE}.txt"
echo $FILENAME
for dir in crab_projects/*${DIR}*; do
  echo "$dir"
  crab kill  $dir >> $FILENAME 
done
