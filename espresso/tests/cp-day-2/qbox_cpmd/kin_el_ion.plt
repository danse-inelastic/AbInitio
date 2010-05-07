#!/bin/bash
gnuplot -persist <<EOF
p "<grep -h '<ekin_e>' $*" u 2 w l, "<grep -h 'ekin_ion' $*" u 2 w l
EOF

