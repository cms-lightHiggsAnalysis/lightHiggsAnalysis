mapfile -t myArray < HaddList.txt
for i in "${myArray[@]}"
  do
    TXTName=${i}
    Name="${TXTName%%.*}"
    echo $Name
    prefix="Hadd "$Name".root"
    echo $prefix
    sed -i ':a;N;$!ba;s/\n/\t/g' ${i}
    sed "s/^/$prefix /" $i >tmp
    mv tmp $i
    SHName=$Name".sh"
    cp $i  $SHName
    chmod 777 $SHName
  done
