#!/usr/bin/env python
#--docstring------------------------------
"""
User inputs for pyphon program.

The entries are used to generate the INPHON file used by Phon.

"""

# Options for symmetrization of force-constant matrix
lsymm    =    False
nti      =    50             # number of iterations for symmetrization

# number of ions types and masses
ntypes   =    1              # number of atom types
mass     =    (50.9415,)      # list of atom masses

# generate superlattice
lsuper   =    False          # flag for generation of supercell
ndim     =    (4,4,4)        # multipliers for supercell
disp     =    100            # inverse of displacement amplitude

# Reciprocal q points section for dispersion calculation
lrecip   =    True           # flag for calculation of dispersions
nd       =    5              # number of segments
npoints  =    25             # number of points per segment

qi = ((0.5,  -0.5,   0.5),   # tuple of segment start q points
      (0.0,   0.0,   0.5), 
      (0.25,  0.25,  0.25), 
      (0.0,   0.0,   0.5))  

qf = ((0.5,  -0.5,   0.5),   # tuple of segment end q points
      (0.0,   0.0,   0.5), 
      (0.25,  0.25,  0.25),
      (0.0,   0.0,   0.5), 
      (0.0,   0.0,   0.0))
 
# free energy calculation
lfree        =    False      # set to true for free energy calculation
temperature  =    300

# density of states
lgamma = False               # set this to True to force the Monckhorst-Pack grid to go through Gamma 
# Divisors for the Monckhorst-Pack grid
# If these are defined, the M-P grid will be calculated (lengthy)
# If QA is set to a negative number, the grid will be read from the file QPOINT
#  QA = -1; QB = 50 ; QC = 50

dosin     =    0     # start bin energy for DOS in THz
dosend    =    10    # end bin energy for DOS in THz
dosstep   =    0.1   # bin step for DOS in THz
dossmear  =    0.1   # smearing of the DOS in THz

# write force constant matrix
lforceout =    True  # set this to True to write the force-cosntant matrix

# cutoff in real space
rmax      =    15    # 

# inverse power calculation
linverse  =    False
alpha     =    5.86  
a         =    1.77

# verbosity
iprint    =    2     












