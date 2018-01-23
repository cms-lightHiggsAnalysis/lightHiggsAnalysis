mapfile -t myArray < TableMC.txt
afterfix=".py"
for i in "${myArray[@]}"
  do
    JobName=${i}
    ReplaceJobName=${JobName//\//\\/}
    Name1=$(cut -d/ -f2 <<<"${JobName}")
    Middle="_"
    Name2=$(cut -d/ -f3 <<<"${JobName}")
    NameToKeep=${Name2: (-3)}
    SubVersion="Version8"
    Name=$Name1$Middle$NameToKeep$SubVersion
    echo $Name
    file="$Name"$afterfix
    
    sed "s/NAME/$Name/g" "crabMCTemplate.py" > $file
    sed -i "s/INPUTDATASET/$ReplaceJobName/g" $file 
    crab submit -c $file
  done
#crabDataTemplate.py or crabMCTemplate.py
#TableData.txt or TableMC.txt
#Comment or uncomment crab submit -c $file

