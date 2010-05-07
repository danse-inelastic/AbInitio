#!/bin/bash
gnuplot -persist <<EOF
p "<grep -h econst $*" u 2 w l, "<grep -h '<etotal>' $*" u 2 w l
EOF
