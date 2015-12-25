#!/bin/bash


start=$1

j=$start

l=1 #line number


imin=10000 # minimum shootpoint value
imax=0 # maximum shootpoint value
#cp acc.txt new/acc.txt 
for v in *shootpt.pdb ; do

    n=${v%%.*}
    if (($n>$imax)) ; then
	imax=$n
    fi
 
    if (($n<$imin)) ; then
	imin=$n
    fi 
done

echo "imin= $imin"
echo "imax= $imax"

l=1

for i in `seq $imin $imax`; do

     n=$((i+j))
#    echo "${i%%.*}"
#    n=${i%%.*}
#    n=$((n+j))
    echo $n

    echo "cp $i.shootpt.pdb new/$n.shootpt.pdb"
#    cp $i new/$n.${i#*.}
    
    echo "sed "${l}s/$i/$n/1" acc.txt > new.acc.txt"
    sed "${l}s/$i/$n/1" acc.txt > new.acc.txt
    cp new.acc.txt acc.txt

    l=$((l+1))

done

cp acc.txt new/acc.txt
cp basin_evals.txt new/basin_evals.txt
cp h_basin.txt new/h_basin.txt
