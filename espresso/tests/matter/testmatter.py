#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import matter
#from diffpy.Structure import Structure
from qecalc.qetask.qeparser.pwinput import PWInput

from matter import Structure, Lattice, Atom

testSCFString = """
"""

def testMatter():

    pwInput = PWInput(config = testSCFString)


    struct = matter.Structure()
    #struct.read('data/graphite.cif', format='cif')
    #struct.read('data/PbTe.cif', format='cif')
    struct.read('data/Ni.stru', format='pdffit')
    #struct.read('data/CdSe-wurtzite.stru', format='pdffit')


    print struct
    print struct.lattice.base

    #struct = Structure(filename='data/Ni.stru')
    #struct = Structure(filename='data/graphite.cif')
    #struct = Structure(filename='data/PbTe.cif')
    #struct = Structure(filename='data/CdSe-wurtzite.stru')
    #struct = Structure(pbte)


    # does not work well in matter:
    #s = pbte.writeStr(format='cif')

    at1 = Atom('V', [0., 0., 0.])
    at2 = Atom('V', [0.5, 0., 0.])
    at3 = Atom('V', [0., 0.5, 0.])
    at4 = Atom('V', [0., 0., 0.5])
    at5 = Atom('V', [0.5, 0.5, 0.])
    at6 = Atom('V', [0., 0.5, 0.5])
    at7 = Atom('V', [0.5, 0., 0.5])
    at8 = Atom('V', [0.5, 0.5, 0.5])

    at9 = Atom('V', [0.25, 0.25, 0.25])
    at10 = Atom('Fe', [0.75, 0.25, 0.25])
    at11 = Atom('V', [0.75, 0.75, 0.25])
    at12 = Atom('Fe', [0.25, 0.75, 0.25])

    at13 = Atom('Fe', [0.25, 0.25, 0.75])
    at14 = Atom('V', [0.75, 0.25, 0.75])
    at15 = Atom('Fe', [0.75, 0.75, 0.75])
    at16 = Atom('V', [0.25, 0.75, 0.75])

    struct2 = Structure( [ at1, at2, at3, at4, at5, at6, at7, at8, at9, at10, at11, at12, at13, at14, at15, at16], lattice = Lattice(2, 2, 2, 90, 90, 90) )

    #print struct

    massList = [50., 55.]
    psList  = ['ps1', 'ps2']

    #massList = [1, 2, 3, 4, 5, 6,1, 2, 3, 4, 5, 6]
    #psList  = ['ps1', 'ps2', 'ps2', 'ps3', 'ps4','ps1', 'ps2', 'ps2', 'ps3', 'ps4']

#    pwInput.structure.load(source = 'diffpy', ibrav = 0, structure = struct, \
#                           massList = massList, psList = psList)

    #pwInput.structure.load(source = 'diffpy', structure = struct )

    #print pwInput.structure.atomLabels()

    pwInput.structure.load(source = 'diffpy', structure = struct, ibrav = 2, \
                           massList = massList, psList = psList)


#    pwInput.structure.setStructureFromDiffpyStructure(struct, \
#                                                            massList = massList,\
#                                                            psList = psList)

#    pwInput.structure.setReducedStructureFromDiffpyStructure(struct, ibrav = 2, \
#                                                            massList = massList,\
#                                                            psList = psList)                                                            
    # this will update string s:                                                   
    #s = ''
    #s = pwInput.structure.toString(string = s)
    
    #pwInput.removeNamelist('system')
    #pwInput.removeCard('atomic_species')    
    #pwInput.removeCard('atomic_positions')
    #pwInput.removeCard('cell_parameters') 
    #pwInput.structure.parseInput()
    
    s = pwInput.structure.toString()

    #pwInput.structure.atomicPositionsType = 'alat'

    #s = pwInput.structure.toString()
    
    #pwInput.structure.save('scf.in')

    print s


if __name__ == "__main__":
    
    testMatter()

__author__="Nikolay Markovskiy"
__date__ ="$Jan 13, 2010 11:46:07 AM$"
