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

# This test provides an example of how to perform DOS task (for DOS plotting),
# Phonons on a Grid task (for VNF experiment) and export data obtained into
# IDF format

import Polarizations, Omega2, DOS

from qecalc.qetask.qeparser.pwinput import PWInput
from qecalc.qetask.matdyntask import MatdynTask
from qeutils import kmesh


def genQgridinfo(filename, nqGrid, recipLattice):

    s = ''
    for i in range(3):
        s = s + 'b%d = %# .8f , %# .8f , %# .8f\nn%d = %d\n'%(i+1, \
        recipLattice[i,0], recipLattice[i,1], recipLattice[i,2], i+1, nqGrid[i])

    open(filename, 'w').write(s)

def testIDF():


# matdyn setting:
    settingString = """
[matdyn.x]
matdynInput = matdyn.in
flfrc = mgb2666.fc
"""

#*********************Initialize matdyn*****************************************
    matdyn = MatdynTask( configString = settingString)
    matdyn.syncSetting()

#*******************************************************************************
    print 'generating "Phonons on the Grid"'

# scf.in with fields relevant to structure information
    pwString = """
&SYSTEM
    ibrav = 4
    celldm(1) = 5.78552800736
    celldm(3) = 1.13577682118
    nat = 3
    ntyp = 2
/

ATOMIC_SPECIES
 Mg 24.305 mg_6.ncpp
 B  11.0   B.pbe-n-van_ak.UPF

ATOMIC_POSITIONS ALAT
 Mg     0.  0.          0.
 B      0.5 0.28867513  0.56788841
 B      0.  0.57735027  0.56788841
"""

    nqGrid = [4,4,4]

    #initialize structure:
    pwInput = PWInput(config = pwString)

    pwInput.parse()

    qpoints = kmesh.kMeshCart(nqGrid,pwInput.structure.lattice.reciprocalBase())

    #update qpoints and launch matdyn
    matdyn.input.qpoints.set(qpoints)
    matdyn.input.save()
    matdyn.launch()
    
    matdyn.output.parse()
    Pols, Freqs, qPoints = matdyn.output.property('multi phonon')

    # convert to Hz**2 and save
    Omega2.write( (Freqs*29979245800.0*2.*3.14159265)**2,'Omega2.idf','')
    Polarizations.write(Pols, 'Polarizations.idf','Polarizations')    

    idfPolData = Polarizations.read('Polarizations.idf')
    idfOmega2Data = Omega2.read('Omega2.idf')


    #save lattice/grid information and make it angstrem compatible, multiply by 2pi:
    genQgridinfo('Qgridinfo', nqGrid, \
                  pwInput.structure.lattice.diffpy().\
                                  reciprocal().base*2.0*3.14159265*1.889725989)
    print idfPolData
    print idfOmega2Data


#*******************************************************************************
    print 'generating "Phonon DOS"'
    nqGrid = [10,10,10]

    matdyn.input.qpoints.setAutomatic(nqGrid)
    matdyn.input.save()
    matdyn.launch()
    matdyn.output.parse()

    axis, dos = matdyn.output.property('phonon dos')

    # save DOS in THz
    DOS.write(axis*0.0299792458, dos, filename = 'DOS', comment = '')

    idfDOSData = DOS.read('DOS')
    print idfDOSData
    #import numpy
    #import numpy.linalg as numalg
    #print 2 * 3.14159265 * numalg.inv(numpy.transpose(pwInput.structure.lattice.diffpy().base/1.889725989))
#*********************Cleaning**************************************************
    import os
    os.system('cat ./Qgridinfo')
    os.system('rm Omega2.idf Polarizations.idf matdyn.modes  matdyn.out \
                                          matdyn.dos matdyn.freq DOS Qgridinfo')

if __name__ == "__main__":

    testIDF()
    

__author__="Nikolay Markovskiy"
__date__ ="$Jan 15, 2010 3:17:39 PM$"
