&CONTROL
    title = Ni
    calculation = 'nscf'
    pseudo_dir = './'
    tprnfor = .true.
    prefix = 'ni'
    outdir = 'temp/'
/
&ELECTRONS
    mixing_beta = 0.7
    conv_thr = 1.0d-8
/
&SYSTEM
    nspin = 2
    ecutwfc = 27.0
    occupations = 'smearing'
    celldm(1) = 6.65
    ibrav = 2
    starting_magnetization(1) = 0.5
    degauss = 0.02
    smearing = 'gauss'
    nat = 1
    ntyp = 1
    ecutrho = 300.0
/
K_POINTS AUTOMATIC
 4 4 4 1 1 1
ATOMIC_POSITIONS
 Say Hi! :)
