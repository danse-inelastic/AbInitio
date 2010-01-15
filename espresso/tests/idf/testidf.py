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


import Polarizations, Omega2

from qecalc.qetask.matdyntask import MatdynTask

if __name__ == "__main__":

    settingString = """
[matdyn.x]
flvec = matdyn.modes
"""
    matdyn = MatdynTask( configString = settingString)
    matdyn.output.parse()

    Pols, Omegas, qPoints = matdyn.output.property('multi phonon')

    #print Omegas

    Omega2.write( ((Omegas/4.1357)**2)*1.e+24,'Omega2.idf','')

    Polarizations.write(Pols, 'Polarizations.idf','Polarizations')
    idfdata = Polarizations.read('Polarizations.idf')
    idfdataOmega = Omega2.read('Omega2.idf')
    print idfdata
    print idfdataOmega
    import os
    os.system('rm Omega2.idf Polarizations.idf')
  #  print qPoints
  #  print len(qPoints[:,1])
    

__author__="Nikolay Markovskiy"
__date__ ="$Jan 15, 2010 3:17:39 PM$"
