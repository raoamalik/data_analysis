 grep -A104 "Q			P" output.dat | grep -v -- "^--$" | grep -vE "Q		 P" > qp.dat
 cut -c3-35 qp.dat > qcoord.dat

# cut -d" " -f2- qcoord.dat > qcoord2.dat
# cut -d" " -f2- qcoord2.dat > qcoord3.dat

