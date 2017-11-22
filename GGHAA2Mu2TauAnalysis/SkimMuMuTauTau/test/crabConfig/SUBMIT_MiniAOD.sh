#!/bin/bash
#parse arguments
if [ $# -ne 2 ]
    then
    echo "Usage: GIVE ME A NEW TAG and then a Old Tag!"
    exit 0
fi

new_tag="${1}"
old_tag="${2}"

#export SCRAM_ARCH=slc6_amd64_gcc481
eval `scramv1 runtime -sh`

DATE=`date +TIME_%H-%M-%S__DATE_%y-%m-%d`
FILENAME="CRAB_SUBMIT_ALL_${new_tag}_${DATE}.txt"

sed -i "s|${old_tag}|${new_tag}|g" crabConf*

echo "Starting"
echo "" >> $FILENAME
crab submit crabConfig_SingleMu1.py >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_SingleMu2.py >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_DoubleMu1.py >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_DoubleMu2.py >> $FILENAME

echo "DATA DONE"

echo "" >> $FILENAME
crab submit crabConfig_DYHighMass.py >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_DYLowMass.py >> $FILENAME




