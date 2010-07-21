#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2010  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from qecalc.qetask.matdyntask import MatdynTask
from qecalc.qetask.pwtask import PWTask
from qeutils.phdispersion import PHDispersion

import numpy

configString = """
[matdyn.x]
# input/output files relevant to multiple phonon calculation
matdynInput:   matdyn.in
matdynOutput:  matdyn.out
flfrc:         mgb2666.fc
"""

if __name__ == "__main__":
    matdyn = MatdynTask(configString = configString)
    pw = PWTask(configString = configString)
    pw.input.parse()

    matdyn.input.parse()
    matdyn.syncSetting()


    phDisp = PHDispersion(pw.input.structure.lattice, matdyn)

    #phDisp.launch('M', 'Gamma', 'A', 'L', 200, 200, 200)
    nP = 100
    phDisp.launch('Gamma', 'K', 'M', 'Gamma', 'A','H', 'L', 'A', nP, nP, nP, nP , nP , nP, nP)
    #phDisp.launch('Gamma', 'H', nP)
    #phDisp.launch('Gamma', 'K', 'M', 'Gamma', nP , nP, nP)
    #phDisp.launch( 'Gamma', 'K',  nP )

   # phDisp.solveIndex((0,0,0), [-0.1,0.2,0.4])

#    matdyn.output.parse()
#
#    Pol, Omegas, qPoints = matdyn.output.property('multi phonon')
#    Pol2 = Pol.reshape(Pol.shape[0], Pol.shape[1], Pol.shape[2]*Pol.shape[3])
#    print Pol2
#    print Omegas
#    print qPoints
#    print Pol2.shape
#    #print numpy.angle( Pol2[5,3,:], Pol2[6,3,:] )
#    a = Pol2[40,3,:]
#    b = Pol2[41,3,:]
#    print a
#    print b
#    print (abs(a)-abs(b)).sum()
#    print abs((a.conj()*b).sum())
#
#    thrsh_scl = 0.8
#    #dispersion = [list(Omegas[0])]
#    idxs = [range(Omegas.shape[1])]
#    for k in range(1, Omegas.shape[0]):
#        omg = []
#        idx = []
#        for i, o1 in zip( idxs[k-1], Omegas[k-1]):
#            for j, o2 in enumerate(Omegas[k]):
#                if ( abs((Pol2[k,j,:].conj()*Pol2[k-1,i,:]).sum()) >  thrsh_scl):
#                    #omg.append(o2)
#                    idx.append(j)
#                    break
#        if len(omg) != 9:
#            print "OMG!!! ", k, Omegas[k-1], Omegas[k],  omg, '\n'
#        #dispersion.append(omg)
#        idxs.append(idx)
#    dispersion = []
#    for k, idx in enumerate(idxs):
#        dispersion.append( Omegas[k, idx])

    #print dispersion
    #print dispersion[0]
    #phDisp.dispersion =  numpy.array(dispersion)

    phDisp.dispersion = phDisp.dispersion*0.1239
    phDisp.solveAllCrossings(thrsh_scl = 0.75)
    phDisp.plot()
    #print phDisp.dispersion
