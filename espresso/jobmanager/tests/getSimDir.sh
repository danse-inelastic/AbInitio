#!/bin/bash

# Returns whole simulation directory which includes input, output files and temp directory
# to the local directory

ssh dexity@foxtrot.danse.us 'cd /home/dexity/espresso; tar czf Ni.tar.gz Ni;'
scp dexity@foxtrot.danse.us:/home/dexity/espresso/Ni.tar.gz /home/dexity/temp/
cd /home/dexity/temp/
tar xzf Ni.tar.gz
rm Ni.tar.gz
ssh dexity@foxtrot.danse.us 'rm /home/dexity/espresso/Ni.tar.gz'