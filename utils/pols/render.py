############################################################
# experimental stuff...

import numpy as np
from crystal import UnitCell as UC
from crystal.UnitCell import *
from crystal.Atom import Atom



### define the unit cell:

uc = UnitCell()
a = 4.709
vectors = [(a,0,0), (0,a,0), (0,0,a)]
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
pos6 = (0.00, 0.75, 0.50)
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

from crystal.crystalUtils.plotUC import plotUnitCell as plotUC

plot = plotUC(uc)

reps = plot.GetDictOfSpecies()
repSi = reps['Si']
repV = reps['V']
repSi.SetRadius(0.8)
repV.SetRadius(0.6)
plot.Update()


############################################################
### animate a mode:

# First, we want to pick the mode:

disp = np.real(pols[2][23])
qvec = qpts[3]

# now we set the displacement on each atom to the polarization vector:

for n in range(len(uc)):
    uc[n].getAtom().displacement = disp[n]

# generate supercell (use PhonCalc):
from AbInitio.AbiPhon.PhonCalc import PhonCalc

pc = PhonCalc(uc, name='V3Si', supersize=[2, 2, 2], abiCalc=None, dosmin=0.0, dosmax=50.0, dosstep=0.2, dossmear=0.2, temperature=300)

pc.genSupercell()
suc = pc._supercell

plot = plotUC(suc)
reps = plot.GetDictOfSpecies()
repSi = reps['Si']
repV = reps['V']
repSi.SetRadius(0.8)
repV.SetRadius(0.6)
plot.Update()


## calculate the sin(q.r) phase factor

natoms = len(suc)
scale = np.array(pc._supersize)

loa = plot.GetListOfAtoms()
refpos = loa.GetCartesianPositions()

phase = np.zeros((natoms))

for n in range(natoms):
    pos = suc[n].getPosition()
    pos = pos * scale
    phase[n] = 2*pi*np.dot(qvec, pos) + randphi
cosqr = np.cos(phase)
sinqr = np.sin(phase)



############################################################
# Render the distorted supercell
##
loa = plot.GetListOfAtoms()
#loa.SetCartesianVelocities(su)
loa.SetCartesianPositions(refpos+su) # get current positions to add to for superposition
plot.AddVelocities()
plot.GetAvatarDict()
vel= plot.GetAvatarDict()['vtkVelocities.#1']
vel.GetScale()
vel.SetScale(10)
vel.SetRadius(0.05)
plot.Update()

## check the mean square displacement in this mode:

np.sqrt((su * su).sum(1)).sum()/64

###


vels = np.zeros((natoms,3))
vels[:,2] = sinqr

plot.velocities._vectors = vels

sloa.SetCartesianVelocities(vels)
