 #!/bin/bash

inputFilePrefix="INPUT_FILE_"
rootFilePrefix="root://eoscms/"
for i in /eos/cms/store/group/phys_higgs/HiggsExo/ktos/*
do 
  dirLs="$(eos ls $i | grep "JAN1")"
  if [ -z "$dirLs" ]; then
    echo "Skipped $i"
    continue
  fi
  echo "$dirLs" | while IFS= read -r line
  do 
    echo ""
    echo "i=$i"
    echo "line=$line"
    filename="${inputFilePrefix}${line}.txt"
    echo "filename=$filename"
    num=$(eos ls $i/$line)
    rootFileTotalPrefix="$rootFilePrefix$i/$line/$num/0000/"
    echo "rootFileTotalPrefix=$rootFileTotalPrefix"
    eos ls "$i/$line/$num/0000/" | grep Region > SummedWeightsFiles/${filename}
    sed -i -e "s|^|${rootFileTotalPrefix}|"   SummedWeightsFiles/${filename}
    done  
done

