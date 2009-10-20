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
# This example will find optimized lattice parameters at given temperatures
# grom QHA. only .fc and properly configured matdyn.in are required. Phonon
# free energy is fitted by a strait line. Total electron energy is fitted by
# forth order polynomial. Calculated values are saved to files as is at original
# volumes only. Eyvaz's QHA program is not required
#

from qecalc import QECalc
from thermo.thermodyn import PhononThermodynamics
import numpy
import os
from volufit import ValueFit
from voluphon import VoluPhon


#MgB2

a_range = numpy.array([5.7878685053999996, 5.7908517081057544, 5.7936410403163601, 5.7963685196991213, 5.7990925545725851, 5.8018743864089419, 5.8044582699567719, 5.8073022753593788, 5.8100511286079417, 5.812768996674941, 5.8155050491948899, 5.8182582903036675, 5.8210314800616452, 5.8236598575075487])

c_a_range = numpy.array([1.13684331855, 1.1373574367579258, 1.1377934839829309, 1.1382631316319463, 1.1387319742719577, 1.1391639469724812, 1.1397097036224211, 1.1400994348409152, 1.1405424190711311, 1.1410008641822587, 1.1414458227913356, 1.1418778930838909, 1.1422954623649801, 1.1427955331216353])

e_total = numpy.array([-15.379301119999999, -15.3792974, -15.379287189999999, -15.37927107, -15.37924881, -15.379220439999999, -15.379185870000001, -15.379145449999999, -15.37909943, -15.379047269999999, -15.378989519999999, -15.378925730000001, -15.37885619, -15.37878098])

fc_name = 'mgb2666.fc'
matdynfldos = 'mgb2666.phdos'

indexRange = [0,2,4,8,10,12]
prcntVol = 2.0*numpy.array(indexRange)/1000.0

# Volume range after fitting of energies
nVolumePoints = 1000.0
finePrcntVolMax = 42.0/1000.0
finePrcntVol = numpy.linspace(prcntVol[0], finePrcntVolMax, nVolumePoints)

temperatures = [1, 10, 50,100,150,200,250,300,450,720]

class FreeEnergy():
# total free energy class at some temperature
    def __init__(self, prcntVol, e_total, phonon_energy):
        self.prcntVol = prcntVol
        self.fitEnergy = ValueFit(prcntVol,e_total[indexRange], 'polynom', 4, True)
        self.fitPhononEnergy = ValueFit(prcntVol,phonon, 'polynom', 1, True)
    def totalFreeEnergy(self,v, *params):
        return self.fitEnergy.fittedValue(v) + self.fitPhononEnergy.fittedValue(v)

    def minimizeTotalFreeEnergy(self,minVol, maxVol):
        import scipy.optimize
        brentOut=scipy.optimize.brent(self.totalFreeEnergy,(),
           (minVol, maxVol),\
                          tol = 1.e-5, full_output = 1)
#        print '\n\nBrentOut:'
#        print brentOut
        optPrcntVolume = brentOut[0]
        optFreeEnergy = brentOut[1]
        return optPrcntVolume




if __name__ == "__main__":
    qecalc = QECalc("config.ini")
    voluPhon = VoluPhon("config.ini", prcntVol)
    voluPhon.setA(a_range[indexRange], 'polynom', 1, True)
    c_range = a_range*c_a_range
    voluPhon.setC(c_range[indexRange], 'polynom', 1, True)
    # Generate phonon doses from different fc files
    for i in indexRange:
        os.system('cp ' + str(i) + '_' + fc_name + ' ' + fc_name)
        qecalc.matdynLauncher()
        os.system('cp ' + matdynfldos + ' ' + str(i) + '_' + matdynfldos )

    #plot free energies:
    for temperature in temperatures:
        lines = ''
        phonon = []
        for i,v in zip(indexRange, prcntVol):
            axis, dos = qecalc.getPhononDOS(str(i) + '_' + matdynfldos)
            phonTherm = PhononThermodynamics(axis, dos)
            Cv = phonTherm.Cv(temperature)
            phonon.append(phonTherm.freeEnergy(temperature))
            totalFreeEnergy = e_total[i] + phonon[-1]
            lines = lines + '%f    %f    %f    %f    %f\n'%(v, \
            totalFreeEnergy, e_total[i], phonon[-1], Cv)
        file = open(str(temperature) + '_free_energy.out', 'w')
        file.write(lines)
        file.close()

        freeEnergy = FreeEnergy(prcntVol, e_total, phonon)
        optPrcntVolume = freeEnergy.minimizeTotalFreeEnergy(finePrcntVol[0], finePrcntVol[-1])
        print "\nOptimized lattice parameters:"
        a = voluPhon.a.fittedValue(optPrcntVolume)
        c = voluPhon.c.fittedValue(optPrcntVolume)
        print "Temperature = %fK"%temperature
        print 'Optimized %Volume = ', optPrcntVolume*100.0
        print 'Optimized a = ', a
        print 'Optimized c = ', c
        print '\n'
        lines = ''
        for v in finePrcntVol:
            e = freeEnergy.fitEnergy.fittedValue(v)
            p = freeEnergy.fitPhononEnergy.fittedValue(v)
            totalFreeEnergy = e + p
            lines = lines + '%f    %f    %f    %f\n'%(v, \
            totalFreeEnergy, e, p)
        file = open(str(temperature) + '_fitted_free_energy.out', 'w')
        file.write(lines)
        file.close()




__author__="Nikolay Markovskiy"
__date__ ="$Oct 12, 2009 3:01:51 PM$"

