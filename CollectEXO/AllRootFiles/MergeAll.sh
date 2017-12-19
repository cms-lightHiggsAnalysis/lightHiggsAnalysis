mapfile -t myArray < HaddList.txt
for i in "${myArray[@]}"
  do
    TXTName=${i}
    Name="${TXTName%%.*}"
    echo $Name
    prefix="hadd "$Name".root"
    echo $prefix
    SHName=$Name".sh"
    cp ./SUB/$i  ./SHFiles/$SHName
    sed -i ':a;N;$!ba;s/\n/ /g' ./SHFiles/$SHName
    sed "s/^/$prefix /" ./SHFiles/$SHName >tmp
    mv tmp ./SHFiles/$SHName
    chmod 777 ./SHFiles/$SHName
    source ./SHFiles/$SHName
  done
