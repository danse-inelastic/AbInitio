#! /usr/bin/env python

__doc__ = """Test of S(Qvec, E) calculator to plot the scattering intensity
along the [100] phonon dispersions in B2 FeAl."""

import sys
import numpy
import pylab
from crystal.UnitCell import *
from crystal.Atom import *
import AbInitio.kernelGenerator.SqeCalculator

def test_disp100(nq, ne):
    """Calculates the scattering intensity for phonon dispersions
    for B2 FeAl along q // [100].
    nq = number of q points to calculate.
    ne = number of e points to calculate."""

    uc = UnitCell( )
    at1=Atom(symbol='Fe', mass=57) ; pos1=(0.0,0.0,0.0)
    at2=Atom(symbol='Al') ; pos2=(0.5,0.5,0.5)
    site1 = Site(pos1, at1)
    site2 = Site(pos2, at2)
    uc.addAtom( at1, pos1, "Fe1" )
    uc.addAtom( at2, pos2, "Al1" )
    print uc

    kptlist = uc.getMonkhorstPackGrid((20,20,20)).reshape(8000,3)
    sqecalc = AbInitio.kernelGenerator.SqeCalculator.SqeCalculator(uc, kpoints=kptlist)

    sqecalc.readIDFeigenvectors(filename='polarizations.idf')
    sqecalc.readEigenvaluesFromIDFomega2(filename='omega2s.idf')

    sqecalc._DebyeWallerCalculator._energies = sqecalc._energies
    sqecalc._DebyeWallerCalculator._polvecs = sqecalc._polvecs

    estart = 0.0
    deltae = 50.0 / ne
    sqecalc._etransferTol = deltae

    deltaqx = 3.0 / nq
    sqecalc._qtransferTolRadius = deltaqx
    qstart = numpy.array([0.0, 0.0, 0.0])
    deltaq = numpy.array([deltaqx, 0.0, 0.0])

    sqe = numpy.zeros((nq,ne), dtype='float')

    for iq in range(nq):
        for ie in range(ne):
            qtransfer = qstart + iq * deltaq
            etransfer = estart + ie * deltae
            print "Q= %s , E= %s" % (qtransfer, etransfer)
            sqe[iq,ie] = sqecalc.calcSqeCohCreateAllmodes(qtransfer, etransfer)
            print "S(Q,E) = ", sqe[iq,ie]

    pylab.imshow(sqe)
    pylab.show()
    end = raw_input()
    return

if __name__ == '__main__':
       #try:
       nq = int(sys.argv[1])
       ne = int(sys.argv[2])
       test_disp100(nq, ne)
#except:
#    print "Usage:", sys.argv[0], "num_Q_points, num_E_points"

# End of file

