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

#sed -i "s|${old_tag}|${new_tag}|g" crabConf*

echo "" >> $FILENAME
crab submit crabConfig_QCD_1000toInf.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_120to170.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_15to20.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_170to300.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_20to30.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_300to470.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_30to50.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_470to600.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_50to80.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_600to800.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_800to1000.py  >> $FILENAME
echo "" >> $FILENAME
crab submit crabConfig_QCD_80to120.py  >> $FILENAME
echo "" >> $FILENAME
