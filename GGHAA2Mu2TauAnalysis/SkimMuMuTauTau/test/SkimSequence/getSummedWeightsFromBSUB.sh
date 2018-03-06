 #!/bin/bash

fileName='/afs/cern.ch/user/k/ktos/NewSkimDir/CMSSW_8_0_30/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/SkimSequence/SummedWeightsFiles/SummedWeightsValues.out'
for i in /afs/cern.ch/user/k/ktos/NewSkimDir/CMSSW_8_0_30/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/SkimSequence/BSUB/*
do 
  temp=${i%%_SumW}
  dirName=${temp#*BSUB/}
  echo "$dirName"
  valueLine=$(grep "value=" ${i}/LSF*/STDOUT | tail -1)
  value=${valueLine##*=}
  echo "  $value"
  echo "${dirName} $value" >> SummedWeightsFiles/SummedWeightsValues.out
done

