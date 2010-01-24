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

import idf.Polarizations, idf.Omega2, idf.DOS   # Change imports as you need

from qecalc.qetask.qeparser.pwinput import PWInput
from qecalc.qetask.matdyntask import MatdynTask
#from qeutils import kmesh

C       = 29979245800.0
TWO_PI  = 2.*3.14159265

def testIDF():
    settingString = """
[matdyn.x]
flvec = matdyn.modes
flfrq = matdyn.freq
fldos = matdyn.dos
"""

# No flfrc or matdynInput are needed for conversion to DOS
#flfrc = mgb2666.fc
#matdynInput = matdyn.in

#*********************Initialize matdyn*****************************************
    matdyn = MatdynTask( configString = settingString)
    matdyn.syncSetting()
    matdyn.output.parse()
    axis, dos = matdyn.output.property('phonon dos')

    # save DOS in THz
    idf.DOS.write(axis*0.0299792458, dos, filename = 'DOS.3', comment = '')

    idfDOSData = idf.DOS.read('DOS.3')
    print idfDOSData

if __name__ == "__main__":
    testIDF()


__date__ = "$Jan 24, 2010 4:12:16 AM$"


