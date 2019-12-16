#!/bin/bash

hn=$1
tn=$2

counter=$1

mkdir mtoz
mkdir charge

cp ../MOPAC2012.exe mtoz
while [ "$counter" -le "$tn" ]

do

cd traj$counter

tail -n104 output.dat > temp1
awk '{print $1}' temp1 > x.dat
awk '{print $2}' temp1 > y.dat
awk '{print $3}' temp1 > z.dat

paste x.dat ../../column1.dat y.dat  ../../column1.dat z.dat ../../column1.dat > temp2
paste -d' ' ../../atoms_3a.dat temp2 > temp3
cat ../../mopac_header.dat temp3 > last$counter.dat

cp last$counter.dat ../mtoz

rm temp1 temp2 temp3 x.dat y.dat z.dat

cd ../mtoz
{
yes "" | ./MOPAC2012.exe last$counter.dat
} &> /dev/null

sed -n '646,749p' last$counter.out | awk '{print $2 "\t" $3}' > ../charge/mtoz$counter.dat
cd ../

let counter+=1

done