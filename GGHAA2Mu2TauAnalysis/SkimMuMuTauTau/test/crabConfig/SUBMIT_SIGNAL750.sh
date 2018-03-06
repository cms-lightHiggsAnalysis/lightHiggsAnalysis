#!/bin/bash
#parse arguments
if [ $# -ne 2 ]
    then
    echo "Usage: GIVE ME A NEW TAG and then a Old Tag!"
    exit 0
fi

new_tag="${1}"
old_tag="${2}"

eval `scramv1 runtime -sh`

echo "Set Environment"
DATE=`date +TIME_%H-%M-%S__DATE_%y-%m-%d`
FILENAME="CRAB_SUBMIT_ALL_${new_tag}_${DATE}.txt"

#sed -i "s|${old_tag}|${new_tag}|g" crabConf*
echo "Sed-ed"

mass=5
oldMass=21
while [ $mass -lt 22 ]; do
  oldMassString1="h750a$oldMass"
  massString1="h750a$mass"
  oldMassString2="M-${oldMass}_"
  massString2="M-${mass}_"
  sed -i "s|${oldMassString1}|${massString1}|g" crabConfig_H750.py
  sed -i "s|${oldMassString2}|${massString2}|g" crabConfig_H750.py
  echo "Submitting $massString1 $massString2"
  crab submit crabConfig_H750.py >> $FILENAME
  let oldMass=mass
  let mass=mass+2
  echo "" >> $FILENAME
done
