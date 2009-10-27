#!/bin/bash

# Returns simulation directory with input, output files only to the local directory.
# Temporary directory will NOT be delivered!

ssh dexity@foxtrot.danse.us 'cd /home/dexity/espresso; tar czf Ni.tar.gz Ni/ni.scf.in Ni/ni.scf.out;'
scp dexity@foxtrot.danse.us:/home/dexity/espresso/Ni.tar.gz /home/dexity/temp/
cd /home/dexity/temp/
tar xzf Ni.tar.gz
rm Ni.tar.gz
ssh dexity@foxtrot.danse.us 'rm /home/dexity/espresso/Ni.tar.gz'