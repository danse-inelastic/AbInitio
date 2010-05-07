&CONTROL
  calculation = 'cp',
  title = 'Methane'
  restart_mode = 'from_scratch',
  ndr = 51
  ndw = 51
  dt = 5,
  nstep  = 1,
  iprint = 1,
  pseudo_dir = './',
  outdir = '$ESPRESSO_TMPDIR'
/

&SYSTEM
  ibrav = 1,
  celldm(1) = 16.0
  nat  = 5,
  ntyp = 2,
  ecutwfc = 25,
  ecutrho = 250
  nr1b = 30, nr2b = 30, nr3b = 30
/

&ELECTRONS
  emass = 300.d0,
  emass_cutoff = 3.d0,
  orthogonalization = 'Gram-Schmidt',
  tcg = .true.,
  startingwfc = 'random'
  ampre = 0.02
/

&IONS
  ion_dynamics = 'none'
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

