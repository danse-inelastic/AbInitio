import numpy as np
from numpy import random

#from inelastic.crystal.UnitCell import *
#from inelastic.crystal.crystalIO.converters import *
from crystal.UnitCell import *
from crystal.crystalIO.converters import *
from AbInitio.vasp.parsing.Structure import Structure
from AbInitio.vasp.parsing.SystemPM import * 

run=XMLSystemPM('vasprun.xml_eq')
struct=run.FINAL_STRUCTURE

uc = p4vaspStruct2UnitCell(struct)   
#kptgrid = uc.getFracMonkhorstPackGrid((2,2,2))

from AbInitio.AbiCalc import VaspCalc as VaspCalc
vc = VaspCalc.VaspCalc(unitcell = uc, kpts = (2,2,2), ekincutoff=320, name='V3Si', vaspcmd='mr vasp')

# this should be taken care of in the VaspCalc
#from crystal.crystalIO.converters import unitCell2ListOfAtom
#loa = unitCell2ListOfAtom(uc)
#vc._vasp._SetListOfAtoms(loa)

print vc.getPotEnergy()
print vc.getForces()
print vc.getStress()
