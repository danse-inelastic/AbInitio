# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2010  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

"""
Fixtures for QE parser unit tests
"""

textPW = """
 &control
    calculation='scf'
    restart_mode='from_scratch',
    tprnfor = .true.
    prefix='ni',
    pseudo_dir = '',
    outdir=''
 /
 &system
    ibrav=2,
    celldm(1) =6.65,
    nat=  1,
    ntyp= 1,
    nspin=2,
    starting_magnetization(1)=0.5,
    degauss=0.02,
    smearing='gauss',
    occupations='smearing',
    ecutwfc =27.0
    ecutrho =300.0
 /
 &electrons
    conv_thr =  1.0d-8
    mixing_beta = 0.7
 /


ATOMIC_SPECIES
 Ni  26.98  Ni.pbe-nd-rrkjus.UPF

ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
K_POINTS AUTOMATIC
4 4 4 1 1 1
blah

"""

textCards = """
ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
K_POINTS AUTOMATIC
4 4 4 1 1 1
blah

"""

assertCards = """
ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
K_POINTS (automatic)
 4 4 4 1 1 1
 blah
"""

textProblem = """
CELL_PARAMETERS
   0.993162743  -0.000000000   0.000000000
  -0.496581371   0.860104165  -0.000000000
  -0.000000000  -0.000000000   4.345938530
ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
K_POINTS AUTOMATIC
4 4 4 1 1 1
blah

"""

textMatdyn = """
 &input
    asr='crystal',
    amass(1)=24.305, amass(2)=11.000,
    flfrc='mgb2666.fc'
 /
176
0.000000    0.000000    0.456392    0.000000
0.000000    0.000000    0.447264    0.009128
0.000000    0.000000    0.438137    0.018256
0.000000    0.000000    0.429009    0.027384
0.000000    0.000000    0.419881    0.036511
"""


assertMatdyn = """&INPUT
    asr = 'crystal',
    amass(1) = 24.305,
    amass(2) = 11.000,
    flfrc = 'mgb2666.fc',
/
176
0.000000    0.000000    0.456392    0.000000
0.000000    0.000000    0.447264    0.009128
0.000000    0.000000    0.438137    0.018256
0.000000    0.000000    0.429009    0.027384
0.000000    0.000000    0.419881    0.036511
"""


textDynmat = """
&input  fildyn='mgb2.dynG', asr='simple',
        q(1)=0.0, q(2)=0.0, q(3)=0.0 /
"""


assertDynmat = """&INPUT
    fildyn = 'mgb2.dynG',
    asr = 'simple',
    q(1) = 0.0,
    q(2) = 0.0,
    q(3) = 0.0,
/
"""


textPh  = """
 &inputph
  tr2_ph=1.0d-10,
  amass(1)=24.305,
  amass(2)=11.000,
  prefix='mgb2',
  outdir='/scratch/markovsk/mgb2'
  fildyn='mgb2.dynG',
 /

"""


# File: ph.mgb2.in
assertMgB2  = """Phonons of MgB2
&INPUTPH
    tr2_ph = 1.0d-10,
    reduce_io = .true.,
    amass(1) = 24.305,
    amass(2) = 11.000,
    prefix = 'mgb2',
    fildyn = 'mgb2.dyn',
    recover = .true.,
    outdir = 'temp/',
    trans = .true.,
    ldisp = .true.,
    nq1 = 6,
    nq2 = 6,
    nq3 = 6,
/
"""


textHeader  = """
&INPUTPH
   tr2_ph = 1.0d-12,
   prefix = 'si',
   epsil = .false.,
   trans = .true.,
   zue = .false.,
   outdir = '/scratch/si',
   amass(1) = 28.0855,
   fildyn = 'si.dyn_G',
   fildrho = 'si.drho_G',
/
0.0 0.0 0.0
"""

assertHeader = """&INPUTPH
    tr2_ph = 1.0d-12,
    prefix = 'si',
    epsil = .false.,
    trans = .true.,
    zue = .false.,
    outdir = '/scratch/si',
    amass(1) = 28.0855,
    fildyn = 'si.dyn_G',
    fildrho = 'si.drho_G',
/
0.0 0.0 0.0
"""

# This is not a problem text (just add spaces between commas)
textComma   = """&input
   asr='crystal',  dos=.true.
   amass(1)=26.982538, amass(2)=11.000,
   flfrc='mgalb4666.fc', fldos='mgalb4.666.phdos', nk1=28,nk2=28,nk3=28
/
"""

# Does not correctly parses text with comma: ",nk2 = 28,"; should be: "nk2 = 28,"
assertComma = """&INPUT
    asr = 'crystal',
    dos = .true.,
    amass(1) = 26.982538,
    amass(2) = 11.000,
    flfrc = 'mgalb4666.fc',
    fldos = 'mgalb4.666.phdos',
    nk1 = 28,
    ,nk2 = 28,
    ,nk3 = 28,
/
"""

# File: ni.scf.in
assertFile = """&CONTROL
    calculation = 'scf',
    restart_mode = 'from_scratch',
    tprnfor = .true.,
    prefix = 'ni',
    pseudo_dir = '',
    outdir = '',
/
&SYSTEM
    ibrav = 2,
    celldm(1) = 6.65,
    nat = 1,
    ntyp = 1,
    nspin = 2,
    starting_magnetization(1) = 0.5,
    degauss = 0.02,
    smearing = 'gauss',
    occupations = 'smearing',
    ecutwfc = 27.0,
    ecutrho = 300.0,
/
&ELECTRONS
    conv_thr = 1.0d-8,
    mixing_beta = 0.7,
/
ATOMIC_SPECIES
 Ni  26.98  Ni.pbe-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
K_POINTS (automatic)
 4 4 4 1 1 1
"""

textMain = """&CONTROL
    calculation='scf', restart_mode = 'from_scratch', /
ATOMIC_SPECIES
 Ni  26.98  Ni.pbe-nd-rrkjus.UPF
ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00

K_POINTS (automatic)
 4 4 4 1 1 1
"""

assertNL    = """&CONTROL
    title = hello,
/
"""

assertNL_space_3    = """&CONTROL
   title = hello,
/
"""

assertC_no_arg      = """ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
"""

assertC_arg      = """ATOMIC_POSITIONS (alat)
   Ni 0.00 0.00 0.00
"""


__date__ = "$Jul 29, 2010 12:37:03 PM$"


