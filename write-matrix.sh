#!/bin/bash

hn=$1
tn=$2

counter=$1

while [ "$counter" -le "$tn" ]

do

#echo -e "r105$counter(rows)," >> rows-3b

echo -e "	r105$counter(i)= sqrt ((x(i,105)-x(i,$counter))**2+(y(i,105)-y(i,$counter))**2&
	 +(z(i,105)-z(i,$counter))**2)" >> matrix-3b
echo -e "	if (r105$counter(i)  .le. 3.0) then
        	write (17,*) 'r105$counter, ', r105$counter(i),  i*0.130E-14
        	endif" >> matrix-3b

 

let counter+=1

done

