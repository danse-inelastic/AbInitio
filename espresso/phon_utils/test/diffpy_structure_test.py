#!/usr/bin/env python
import math
import numpy

# To change this template, choose Tools | Templates
# and open the template in the editor.
from diffpy.Structure.structure import Structure
from diffpy.Structure.lattice import Lattice

myAtoms = ['Mg', 'B', 'B']
myLattice = Lattice()
theta = -30.*math.pi/180.
rotMatrix = [[math.cos(theta),-math.sin(theta),0],[math.sin(theta),math.cos(theta),0],[0,0,1]]
#myLattice.setLatPar(1., 1, 1.,gamma=120.0)#,baserot = rotMatrix)
print myLattice.base
a = 5.2
c = 6.3
base = [[1,0,0],[-1./2.,math.sqrt(3.)/2.,0.],[0,0,c/a]]
myLattice.setLatBase(base)
print myLattice.base
print "baserot:"
print myLattice.baserot


myLattice.setLatPar(5.2, 5.2, 6.3,gamma=120.0,baserot = myLattice.baserot)
print myLattice.base

#myLattice.setLatBase(numpy.array(base)*a)
#print myLattice
#print myLattice.abcABG()

myStruct = Structure(myAtoms, myLattice)

# fractional coordinates:
myStruct[0].xyz[:] = [0, 0, 0,]
myStruct[1].xyz[:] = [1, 2, 3,]

print myStruct

