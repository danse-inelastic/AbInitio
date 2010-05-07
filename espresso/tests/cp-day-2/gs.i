# compute the ground state of a 64-atom silicon sample

# define the unit cell 
set cell    20.52 0 0  0 20.52 0  0 0 20.52

species silicon silicon.xml 
atom Si01 silicon   -10.260000   -10.260000   -10.260000
atom Si02 silicon    -7.695000    -7.695000    -7.695000
atom Si03 silicon    -5.130000    -5.130000   -10.260000
atom Si04 silicon    -2.565000    -2.565000    -7.695000
atom Si05 silicon   -10.260000    -5.130000    -5.130000
atom Si06 silicon    -7.695000    -2.565000    -2.565000
atom Si07 silicon    -5.130000   -10.260000    -5.130000
atom Si08 silicon    -2.565000    -7.695000    -2.565000
atom Si09 silicon   -10.260000   -10.260000     0.000000
atom Si10 silicon    -7.695000    -7.695000     2.565000
atom Si11 silicon    -5.130000    -5.130000     0.000000
atom Si12 silicon    -2.565000    -2.565000     2.565000
atom Si13 silicon   -10.260000    -5.130000     5.130000
atom Si14 silicon    -7.695000    -2.565000     7.695000
atom Si15 silicon    -5.130000   -10.260000     5.130000
atom Si16 silicon    -2.565000    -7.695000     7.695000
atom Si17 silicon   -10.260000     0.000000   -10.260000
atom Si18 silicon    -7.695000     2.565000    -7.695000
atom Si19 silicon    -5.130000     5.130000   -10.260000
atom Si20 silicon    -2.565000     7.695000    -7.695000
atom Si21 silicon   -10.260000     5.130000    -5.130000
atom Si22 silicon    -7.695000     7.695000    -2.565000
atom Si23 silicon    -5.130000     0.000000    -5.130000
atom Si24 silicon    -2.565000     2.565000    -2.565000
atom Si25 silicon   -10.260000     0.000000     0.000000
atom Si26 silicon    -7.695000     2.565000     2.565000
atom Si27 silicon    -5.130000     5.130000     0.000000
atom Si28 silicon    -2.565000     7.695000     2.565000
atom Si29 silicon   -10.260000     5.130000     5.130000
atom Si30 silicon    -7.695000     7.695000     7.695000
atom Si31 silicon    -5.130000     0.000000     5.130000
atom Si32 silicon    -2.565000     2.565000     7.695000
atom Si33 silicon     0.000000   -10.260000   -10.260000
atom Si34 silicon     2.565000    -7.695000    -7.695000
atom Si35 silicon     5.130000    -5.130000   -10.260000
atom Si36 silicon     7.695000    -2.565000    -7.695000
atom Si37 silicon     0.000000    -5.130000    -5.130000
atom Si38 silicon     2.565000    -2.565000    -2.565000
atom Si39 silicon     5.130000   -10.260000    -5.130000
atom Si40 silicon     7.695000    -7.695000    -2.565000
atom Si41 silicon     0.000000   -10.260000     0.000000
atom Si42 silicon     2.565000    -7.695000     2.565000
atom Si43 silicon     5.130000    -5.130000     0.000000
atom Si44 silicon     7.695000    -2.565000     2.565000
atom Si45 silicon     0.000000    -5.130000     5.130000
atom Si46 silicon     2.565000    -2.565000     7.695000
atom Si47 silicon     5.130000   -10.260000     5.130000
atom Si48 silicon     7.695000    -7.695000     7.695000
atom Si49 silicon     0.000000     0.000000   -10.260000
atom Si50 silicon     2.565000     2.565000    -7.695000
atom Si51 silicon     5.130000     5.130000   -10.260000
atom Si52 silicon     7.695000     7.695000    -7.695000
atom Si53 silicon     0.000000     5.130000    -5.130000
atom Si54 silicon     2.565000     7.695000    -2.565000
atom Si55 silicon     5.130000     0.000000    -5.130000
atom Si56 silicon     7.695000     2.565000    -2.565000
atom Si57 silicon     0.000000     0.000000     0.000000
atom Si58 silicon     2.565000     2.565000     2.565000
atom Si59 silicon     5.130000     5.130000     0.000000
atom Si60 silicon     7.695000     7.695000     2.565000
atom Si61 silicon     0.000000     5.130000     5.130000
atom Si62 silicon     2.565000     7.695000     7.695000
atom Si63 silicon     5.130000     0.000000     5.130000
atom Si64 silicon     7.695000     2.565000     7.695000
set ecut 6

# electron dynamics: preconditioned steepest descent with Anderson acceleration
set wf_dyn PSDA

# set preconditioning cutoff
set ecutprec 3 

# add randomi amplitudes to wavefunctions to break symmetry
randomize_wf

# compute ground state
run 10 5 10

#save sample
save si64.xml

