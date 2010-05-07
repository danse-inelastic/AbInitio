ch4.sys
set ecut 35
set xc PBE
set wf_dyn PSDA
set ecutprec 3
randomize_wf
set xc PBE
run 0 50
save ch4.xml
