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

from qecalc.qetask.qeparser.pwinput import PWInput
from qecalc.qetask.qeparser.matdyninput import MatdynInput
from qecalc.qetask.matdyntask import MatdynTask
from qeutils import kmesh

def testMatdyn():
    matdyn = MatdynTask( configString = "")

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

    matdynString = """
&INPUT
    asr = 'crystal',
    flfrc = '/home/dexity/espresso/qejobs/643E2QQI/default.fc',
    dos = .true.,
    nk1 = 16,
    nk2 = 16,
    nk3 = 16,
    amass(1) = 26.98,
/
"""
    # Initialize structure (PW input):
    pwInput = PWInput(config = pwString)
    pwInput.parse()

    # Initialize MATDYN input
    matdynInput = MatdynInput(config = matdynString)
    matdynInput.parse()
    print matdynInput.toString()

    # Pass matdyn input to matdyn task
    matdyn.input    = matdynInput
    nl      = matdynInput.namelist("input")

    # Populate grid
    nqGrid  = [int(nl.param("nk1")), int(nl.param("nk2")), int(nl.param("nk3"))]
    qpoints = kmesh.kMeshCart(nqGrid, pwInput.structure.lattice.reciprocalBase())

    # Generate grid point
    matdyn.input.qpoints.set(qpoints)   # Generate Q-mesh
    open("matdyn.2.in", "w").write(matdyn.input.toString())
    #print matdyn.input.toString()


if __name__ == "__main__":
    testMatdyn()

__date__ = "$Jan 24, 2010 10:55:56 PM$"


