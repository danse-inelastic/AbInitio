# run MD steps for a 64-atom silicon sample

# load sample with previously computed ground state
load si64.xml 

set wf_dyn PSDA
set ecutprec 3

set atoms_dyn MD
set dt 40

set thermostat LOWE
set th_temp 500
set th_time 1500

run 750 5

# save sample
save si64md_500.xml

