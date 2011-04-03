#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

import idf.Polarizations
import idf.Omega2

from qecalc.qetask.qeparser.pwinput import PWInput
from qecalc.qetask.matdyntask import MatdynTask
#from qeutils import kmesh

C       = 29979245800.0
TWO_PI  = 2.*3.14159265
TO_THZ  = 0.0299792458
A2B     = 1.889725989   # Angstroms to bohrs
PI      = 3.14159265

def genQgridinfo(filename, nqGrid, recipLattice):

    s = ''
    for i in range(3):
        s = s + 'b%d = %# .8f , %# .8f , %# .8f\nn%d = %d\n'%(i+1, \
        recipLattice[i,0], recipLattice[i,1], recipLattice[i,2], i+1, nqGrid[i])

    open(filename, 'w').write(s)

def testIDF():
    settingString = """
[matdyn.x]
flvec = matdyn.modes
flfrq = matdyn.freq
fldos = matdyn.dos
"""

    # Initialize matdyn
    matdyn = MatdynTask( configString = settingString)
    matdyn.syncSetting()
    matdyn.output.parse()
    Pols, Freqs, qPoints = matdyn.output.property('multi phonon')

    # convert to Hz**2 and save
    idf.Omega2.write( (Freqs*C*TWO_PI)**2,'Omega2','')
    idf.Polarizations.write(Pols, 'Polarizations','Polarizations')

    idfPolData = idf.Polarizations.read('Polarizations')
    idfOmega2Data = idf.Omega2.read('Omega2')

    pwInput = PWInput(filename = "pw.in")
    pwInput.parse()

    # XXX
    #   - Why nqGrid?
    #   - Is there a way to get reciprocal lattice from matter?
    nqGrid = [4,4,4]

    #save lattice/grid information and make it angstrem compatible, multiply by 2pi:
    genQgridinfo('Qgridinfo', nqGrid, \
                  pwInput.structure.lattice.matter().\
                                  reciprocal().base*2.0*PI*A2B)
    print idfPolData
    print idfOmega2Data
    import os
    os.system('cat ./Qgridinfo')


if __name__ == "__main__":
    testIDF()


__date__ = "$Jan 24, 2010 4:34:49 AM$"


#    pwString = """
#&SYSTEM
#    ibrav = 4
#    celldm(1) = 5.78552800736
#    celldm(3) = 1.13577682118
#    nat = 3
#    ntyp = 2
#/
#
#ATOMIC_SPECIES
# Mg 24.305 mg_6.ncpp
# B  11.0   B.pbe-n-van_ak.UPF
#
#ATOMIC_POSITIONS ALAT
# Mg     0.  0.          0.
# B      0.5 0.28867513  0.56788841
# B      0.  0.57735027  0.56788841
#"""

#    nqGrid = [4,4,4]

#    #Initialize structure:
#    pwInput = PWInput(config = pwString)
#
#    pwInput.parse()
#
#    qpoints = kmesh.kMeshCart(nqGrid,pwInput.structure.lattice.reciprocalBase())
#
#    #update qpoints and launch matdyn
#    matdyn.input.qpoints.set(qpoints)
#    matdyn.input.save()
#    matdyn.launch()

    #matdyn.output.toString()
