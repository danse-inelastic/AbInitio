#!/bin/bash
gnuplot -persist <<EOF
fit a*x+b "<grep -h temp_ion $*" u 0:2 via a,b
fit c "<grep -h temp_ion $*" u 0:2 via c
p "<grep -h temp_ion $*" u 2 w l, a*x+b, c
print "Tavg=",c,"   dT/dt=",a
EOF
