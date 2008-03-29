import numpy as np
from AbInitio.AbiPhon.AbiPhonCalc import AbiPhonCalc
from AbInitio.AbiPhon.PhonCalc import PhonCalc
from crystal.UnitCell import *
from crystal.Atom import *
from crystal.crystalIO import converters
from AbInitio.vasp.vasp import VASP


# create a unit cell:

uc = UnitCell( )
vectors = np.array([[2.9, 0, 0],[0, 2.9, 0],[0, 0, 2.9]])
uc.setCellVectors(vectors)

at1 = Atom(symbol='Fe', mass=57) ; pos1 = (0.0,0.0,0.0)
site1 = Site(pos1, at1)
at2 = Atom(symbol='Al') ; pos2 = (0.5,0.5,0.5)
site2 = Site(pos2, at2)

uc.addSite(site1, "Fe1" )
uc.addSite(site2, "Al1" )

print uc

### plot the unit cell ###
from ExcitationSlicer.plot3D import plotUnitCell
plotUnitCell(uc)

### create a VASP calculator:

v = VASP(pw=268, kpts=(2,2,2), xc='pawpbe', name='FeAl', vaspcmd='vasp')

loa = converters.unitCell2ListOfAtom(uc)
v._SetListOfAtoms(loa)
v.atoms()
v.incar['EDIFF']=1.0e-5
v.makePoscarFile()
v.makePotcarFile()

# First, run VASP on the original unit cell
v.GetPotentialEnergy()

# Render the structure using VTK renderering utility
#from ASE.Visualization.VTK import VTKPlotAtoms
#atomplot = VTKPlotAtoms(loa, parent=None)

# Now, we want to perform a Phon calculation.
# Before we can call phon(), we need to set up the INPHON input file

pc = PhonCalc(uc, name='FeAl', supersize=[2,2,2],  qgridsize=[20,20,20], dosmin=0.0, dosmax=20.0, amplitude=0.01)

pc.genSupercell((2,2,2))

suc = pc._supercell
plotUnitCell(suc)

pc.genPhonSupercell()
pc._makePosFiles()

pc.calcForces()
pc.calcPhonons()


