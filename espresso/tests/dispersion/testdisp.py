#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2010  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To change this template, choose Tools | Templates
# and open the template in the editor.

    def _getBravaisType(self):
        centeringAndSystem2bravais = {
        ('P', 'CUBIC'):'simple cubic',
        ('F', 'CUBIC'):'face centered cubic',
        ('I', 'CUBIC'):'body centered cubic',
        ('P', 'HEXAGONAL'):'simple hexagonal',
        ('P', 'TRIGONAL'):'simple trigonal',
        ('P','TETRAGONAL'):'simple tetragonal',
        ('R','TETRAGONAL'):'simple tetragonal',
        ('I','TETRAGONAL'):'body centered tetragonal',
        ('P','ORTHORHOMBIC'):'simple orthorhombic',
        ('C','ORTHORHOMBIC'):'base centered orthorhombic',
        ('A','ORTHORHOMBIC'):'base centered orthorhombic',
        ('F','ORTHORHOMBIC'):'face centered orthorhombic',
        ('I','ORTHORHOMBIC'):'body centered orthorhombic',
        ('P','ORTHORHOMBIC'):'simple monoclinic',
        ('P','MONOCLINIC'):'base centered monoclinic',
        ('P','TRICLINIC'):'triclinic'
        }

defaultPoints = {
                ('P', 'CUBIC'):  { #Simple cubic symmetry points:
                                    'Gamma' : [0.0, 0.0, 0.0],
                                    'Gamma2': [1.0, 1.0, 1.0],
                                    'X'     : [0.5, 0.0, 0.0],
                                    'M'     : [0.5, 0.5, 0.0],
                                    'R'     : [0.5, 0.5, 0.5]
                                  },
                ('F', 'CUBIC'):  {   #Face Centered cubic symmetry points:
                                    'Gamma' : [0.0, 0.0, 0.0],
                                    'Gamma2': [1.0, 1.0, 1.0],
                                    'X'     : [0.5, 0.5, 0.0],
                                    'L'     : [0.5, 0.5, 0.5],
                                    'W'     : [0.5, 0.75, 0.25]
                                 },
                ('I', 'CUBIC'):
}




def getDefaultDispersionPlot(self, branches=None, npointspersegment=20):
    if branches is None:
        branches = range(self.dimension*self.nAtoms)

    matter = self.matter
    sg = matter.sg
    from math import pi

    if centeringAndSystem2bravais[(matter.sg.short_name[0], matter.sg.crystal_system)]

    if sg.number == 225: # fcc
        a = matter.lattice.a
        b = 2*pi/a
        Qpoints = [
            (0,0,0),
            (0,0,b),
            (0,b,b),
            (0,0,0),
            (b/2,b/2,b/2),
            ]

    elif sg.number == 229: # bcc
        a = matter.lattice.a
        b = 2*pi/a
        Qpoints = [
            (b/2,b/2,0),
            (0,0,0),
            (b,0,0),
            (b/2,b/2,b/2),
            (0,0,0),
            ]

    else:
        raise NotImplementedError

    return self.getDispersionPlot(Qpoints, branches, npointspersegment)



if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Jan 22, 2010 3:46:08 PM$"
