#!/bin/bash
#usage: cp_xyz [number of ions] [cp output file] > output.xyz
grep -A$1 'ATOMIC_POSITIONS' $2 | grep -v '\--' | sed "s/ATOMIC_POSITIONS/$1\ncomment/g"
