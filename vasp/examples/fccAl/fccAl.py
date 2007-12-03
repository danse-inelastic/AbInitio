from ASE import Atom,ListOfAtoms

#import VASP as vasp
#from VASP.parser2 import *
#from VASP import potcar
from vasp.vasp import VASP
from vasp.parsing.SystemPM import *

from Atom import *
from UnitCell import *
from io import converters

uc = UnitCell()

vectors = [(4.05, 0.0, 0.0),
           (0.0, 4.05, 0.0),
           (0.0, 0.0, 4.05)]

uc.setCellVectors(vectors)

at1=Atom(symbol='Al') ; pos1=(0.0,0.0,0.0)
at2=Atom(symbol='Al') ; pos2=(0.0,0.5,0.5)
at3=Atom(symbol='Al') ; pos3=(0.5,0.0,0.5)
at4=Atom(symbol='Al') ; pos4=(0.5,0.5,0.0)

site1 = Site(pos1, at1)
site2 = Site(pos2, at2)
site3 = Site(pos3, at3)
site4 = Site(pos4, at4)

uc.addSite(site1, 'Al1')
uc.addSite(site2, 'Al2')
uc.addSite(site3, 'Al3')
uc.addSite(site4, 'Al4')

print uc

loa = converters.unitCell2ListOfAtom(uc)

#atoms = ListOfAtoms([Atom('Al',[0,0,0]),
#                     Atom('Al',[0.5,0.5,0]),
#                     Atom('Al',[0.5,0,0.5]),
#                     Atom('Al',[0,0.5,0.5])])

#cellvectors = [[4.05, 0, 0],
#               [0.0, 4.05, 0],
#               [0.0, 0, 4.05]]
#atoms.SetUnitCell(cellvectors)

v = VASP(pw=150, kpts=(2,2,2), name='fccAl', vaspcmd='vasp')
loa.SetCalculator(v)

import ASE.Utilities.GeometricTransforms as ASEgeotrans
    
initvol = loa.GetUnitCellVolume()

vols, energies = [],[]

for f in [0.85, 0.9, 0.95, 1.0, 1.05, 1.1]:
    ASEgeotrans.SetUnitCellVolume(loa,f*initvol)
    v.SetName('%1.2f_eos' % f)
    print v.GetStress()
    vols.append(loa.GetUnitCellVolume())
    energies.append(loa.GetPotentialEnergy())

import pylab as pl

pl.plot(vols,energies,'ko ')
pl.show()

import ASE.Utilities.EquationOfState as ASEeos

eos = ASEeos.EquationOfState('Murnaghan',vols,energies)
# print the minimum volume, energy, bulk modulus and pressure
print eos         
g = eos.GetPlot()
eos.SavePlot('murn.png') #save the figure as a png
