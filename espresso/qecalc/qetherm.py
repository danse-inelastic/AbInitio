#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Brent Fultz
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To change this template, choose Tools | Templates
# and open the template in the editor.

from qecalc import QECalc
from thermo.thermodyn import PhononThermodynamics
import numpy
import os

#MgB2
e_total = [-15.379301119999999, -15.3792974, -15.379287189999999, -15.37927107, -15.37924881, -15.379220439999999, -15.379185870000001, -15.379145449999999, -15.37909943, -15.379047269999999, -15.378989519999999, -15.378925730000001, -15.37885619, -15.37878098]

fc_name = 'mgb2666.fc'
matdynfldos = 'mgb2666.phdos'

indexRange = [0,2,4,6,8,10,12]
prcntVol = 2.0*numpy.array(indexRange)/1000.0

temperatures = [1, 10, 50,100,300,450,720]

if __name__ == "__main__":
    qecalc = QECalc("config.ini")
    #Generate phonon doses from different fc files
    for i in indexRange:
        os.system('cp ' + str(i) + '_' + fc_name + ' ' + fc_name)
        qecalc.matdynLauncher()
        os.system('cp ' + matdynfldos + ' ' + str(i) + '_' + matdynfldos )

    #plot free energies:
    for temperature in temperatures:
        lines = ''
        for i,v in zip(indexRange, prcntVol):
            axis, dos = qecalc.getPhononDOS(str(i) + '_' + matdynfldos)
            phonTherm = PhononThermodynamics(axis, dos)
            Cv = phonTherm.Cv(temperature)
            phonon = phonTherm.freeEnergy(temperature)
            totalFreeEnergy = e_total[i] + phonon
            lines = lines + '%f    %f    %f    %f    %f\n'%(v, \
            totalFreeEnergy, e_total[i], phonon, Cv)
        file = open(str(temperature) + '_free_energy.out', 'w')
        file.write(lines)
        file.close()

__author__="Nikolay Markovskiy"
__date__ ="$Oct 12, 2009 3:01:51 PM$"

