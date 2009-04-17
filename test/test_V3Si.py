import numpy as np
from numpy import random

#from inelastic.crystal.UnitCell import *
#from inelastic.crystal.crystalIO.converters import *
from crystal.UnitCell import *
#from crystal.crystalIO.converters import *
#from AbInitio.vasp.parsing.Structure import Structure
#from AbInitio.vasp.parsing.SystemPM import * 

#run=XMLSystemPM('vasprun.xml_eq')
#struct=run.FINAL_STRUCTURE

#uc = p4vaspStruct2UnitCell(struct)   
#kptgrid = uc.getFracMonkhorstPackGrid((2,2,2))

uc = UnitCell()
vectors = [(4.70898,0,0), (0,4.70898,0), (0,0,4.70898)]
#vectors = [(1,0,0), (0,1,0), (0,0,1)]
uc.setCellVectors(vectors)
at1 = Atom(symbol='Si')
pos1 = (0,0,0)
at2 = Atom(symbol='Si')
pos2 = (0.5,0.5,0.5)
at3 = Atom(symbol='V')
pos3 = (0.25, 0.50, 0.00)
at4 = Atom(symbol='V')
pos4 = (0.75, 0.50, 0.00)
at5 = Atom(symbol='V')
pos5 = (0.00, 0.25, 0.50)
at6 = Atom(symbol='V')
pos6 = (0.00, 0.75, 0.00)
at7 = Atom(symbol='V')
pos7 = (0.50, 0.00, 0.25)
at8 = Atom(symbol='V')
pos8 = (0.50, 0.00, 0.75)

site1 = Site(pos1, at1)
site2 = Site(pos2, at2)
site3 = Site(pos3, at3)
site4 = Site(pos4, at4)
site5 = Site(pos5, at5)
site6 = Site(pos6, at6)
site7 = Site(pos7, at7)
site8 = Site(pos8, at8)

sites = [site1, site2, site3, site4, site5, site6, site7, site8]
for i in range(8):
    uc.addSite(sites[i], sites[i].getAtom().symbol+str(i))

#all this unit cell stuff should be replaced with diffraction's structure, since that is what olivier is 
#really talking about

from AbInitio.AbiCalc import VaspCalc as VaspCalc
vc = VaspCalc.VaspCalc(unitcell = uc, kpts = (2,2,2), ekincutoff=320, name='V3Si', vaspcmd='vasp')

# this should be taken care of in the VaspCalc
#from crystal.crystalIO.converters import unitCell2ListOfAtom
#loa = unitCell2ListOfAtom(uc)
#vc._vasp._SetListOfAtoms(loa)

print vc.getPotEnergy()
print vc.getForces()
print vc.getStress()
