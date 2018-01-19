#!/bin/bash

#parse arguments
if [ $# -ne 4 ]
    then
    echo "Usage: ./generate.sh cfg_name script_name queue name_addon"
    exit 0
fi

cfg_name=$1
script_name=$2
queue=$3
name_addon=$4
input="/afs/cern.ch/user/k/ktos/NewSkimDir/CMSSW_8_0_30/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/SkimSequence/SummedWeightsFiles/SummedWeights_InputFile.txt"
while IFS= read -r line
do
  echo ""
  echo ""
  echo "line=$line"
  noSuffix="${line%.*}"
  dir="${noSuffix#INPUT_FILE_}"
  dir="${dir}${name_addon}"
  echo "dir=$dir"
  mkdir -p BSUB/$dir
  cd BSUB/$dir
  echo "python file= ${cfg_name}_${dir}.py"
  echo "Script name= ${script_name}_${dir}.sh"
  sed -e "s%TXT_INPUT_FILE%${line}%g" ../../${cfg_name}.py > ${cfg_name}_${dir}.py
  sed -e "s%ANALYZER%${cfg_name}_${dir}%g" -e "s%DIRNAME%${dir}%g" ../../${script_name}.sh > ${script_name}_${dir}.sh
  chmod u+x ${script_name}_${dir}.sh
  bsub -q $queue -J ${cfg_name}_${dir} < ${script_name}_${dir}.sh
  cd ../../
done <"$input"
exit 0
