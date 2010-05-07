# MD with constraints
# start from CH4 sample (with atoms displaced!)
load ch4.xml

# define constraints for bond lengths
constraint define distance d1 C H1 2.08
constraint define distance d2 C H2 2.08
constraint define distance d3 C H3 2.08
constraint define distance d4 C H4 2.08
constraint enforce

angle H1 C H2
angle H1 C H3
angle H1 C H4
angle H2 C H3
angle H2 C H4
angle H3 C H4

set wf_dyn PSDA
set ecutprec 3
run 0 20

save ch4_disp.xml

set atoms_dyn MD
set dt 20
run 200 5

