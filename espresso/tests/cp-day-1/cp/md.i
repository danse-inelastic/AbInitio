&CONTROL
  calculation = 'cp',
  title = 'Methane'
  restart_mode = 'restart',
  ndr = 51,
  ndw = 52,
  dt = 5,
  nstep  = 100,
  iprint = 10,
  pseudo_dir = './',
  outdir = '$ESPRESSO_TMPDIR'
/

&SYSTEM
  ibrav = 1,
  celldm(1) = 16.0
  nat  = 5,
  ntyp = 2,
  ecutwfc = 22,
  ecutrho = 220
  nr1b = 30, nr2b = 30, nr3b = 30
/

&ELECTRONS
  electron_dynamics = 'verlet'
  emass = 300,
  emass_cutoff = 4,
/

&IONS
  ion_dynamics = 'verlet'
/

ATOMIC_SPECIES
 C 12.0 C.pbe-rrkjus.UPF
 H 1.0 H.pbe-rrkjus.UPF

ATOMIC_POSITIONS (bohr)
 C 0.0 0.0 0.0 0 0 0
 H 1.2 1.2 1.2 
 H 1.2 -1.2 -1.2 
 H -1.2 1.2 -1.2 
 H -1.2 -1.2 1.2 

