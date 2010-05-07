#!/bin/bash
gnuplot -persist <<EOF
plot "<grep 'etotal' $1" u 2 w l
EOF
