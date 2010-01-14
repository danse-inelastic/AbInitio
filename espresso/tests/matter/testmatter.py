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
from diffpy.Structure import Structure
from qecalc.qetask.qeparser.pwinput import PWInput

testSCFString = """
"""

def testMatter():

    pwInput = PWInput(config = testSCFString)


    struct = matter.Structure()
    #struct.read('data/graphite.cif', format='cif')
    #struct.read('data/PbTe.cif', format='cif')
    struct.read('data/Ni.stru', format='pdffit')
    struct.read('data/CdSe-wurtzite.stru', format='pdffit')


    #print struct
    #print struct.lattice.base

    #struct = Structure(filename='data/Ni.stru')
    #struct = Structure(filename='data/graphite.cif')
    #struct = Structure(filename='data/PbTe.cif')
    #struct = Structure(filename='data/CdSe-wurtzite.stru')
    #struct = Structure(pbte)


    #s = pbte.writeStr(format='cif')


    massList = [1, 2, 3, 4, 5, 6,1, 2, 3, 4, 5, 6]
    psList  = ['ps1', 'ps2', 'ps2', 'ps3', 'ps4','ps1', 'ps2', 'ps2', 'ps3', 'ps4']

    pwInput.structure.setStructureFromDiffpyStructure(struct, \
                                                            massList = massList,\
                                                            psList = psList)

#    pwInput.structure.setReducedStructureFromDiffpyStructure(struct, ibrav = 2, \
#                                                            massList = massList,\
#                                                            psList = psList)                                                            
                                                            
    s = ''
    s = pwInput.structure.save(string = s)

    print s


if __name__ == "__main__":
    
    testMatter()

__author__="Nikolay Markovskiy"
__date__ ="$Jan 13, 2010 11:46:07 AM$"
