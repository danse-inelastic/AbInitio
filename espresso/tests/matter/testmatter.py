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
from qecalc.pwcalc import PWCalc
from diffpy.Structure import Structure

configString = """
[pw.x]
pwInput: scf.in
"""

scfString = """
&CONTROL
    calculation = 'scf',
    wf_collect = .true.,
    prefix = 'si',
    pseudo_dir = './',
    verbosity = 'high',
    outdir = '/scratch/markovsk/d3paratest',
/

&SYSTEM
    ibrav = 2,
    celldm(1) = 10.20,
    nat = 2,
    ntyp = 1,
    ecutwfc = 24.0,
/

&ELECTRONS
    mixing_beta = 0.7,
    conv_thr = 1.0d-12,
/

ATOMIC_SPECIES
 Si  28.086  Si.pz-vbc.UPF

ATOMIC_POSITIONS (crystal)
 Si 0.00 0.00 0.00
 Si 0.25 0.25 0.25

K_POINTS (automatic)
 4 4 4  1 1 1
"""


testSCFString = """
"""

def testMatter():

    pwcalc = PWCalc( configString = configString )
    pwcalc.pw.syncSetting()




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
    #struct = matter.Structure(filename='data/graphite.cif')
    #struct = Structure(pbte)


    #s = pbte.writeStr(format='cif')

    #print s

#    print struct

#    print struct.lattice.a
#    print struct.lattice.b
#    print struct.lattice.c

#    print struct.lattice.alpha
#    print struct.lattice.beta
#    print struct.lattice.gamma

    massList = [1, 2, 3, 4, 5, 6,1, 2, 3, 4, 5, 6]
    psList  = ['ps1', 'ps2', 'ps2', 'ps3', 'ps4','ps1', 'ps2', 'ps2', 'ps3', 'ps4']

    pwcalc.pw.input.structure.setStructureFromDiffpyStructure(struct, ibrav = 4, \
                                                            massList = massList,\
                                                            psList = psList)

    pwcalc.pw.input.structure.save('qqq.in')
    #print pwcalc.pw.input.structure.diffpy()


if __name__ == "__main__":
    
    testMatter()

__author__="Nikolay Markovskiy"
__date__ ="$Jan 13, 2010 11:46:07 AM$"
