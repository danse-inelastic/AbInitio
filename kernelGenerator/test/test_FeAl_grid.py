from UnitCell import *
from Atom import *
import pylab

uc = UnitCell( )
at1=Atom(symbol='Fe', mass=57) ; pos1=(0.0,0.0,0.0)
at2=Atom(symbol='Al') ; pos2=(0.5,0.5,0.5)

site1 = Site(pos1, at1)
site2 = Site(pos2, at2)

uc.addAtom( at1, pos1, "Fe1" )
uc.addAtom( at2, pos2, "Al1" )

print uc


###

import AbInitio.kernelGenerator.SqeCalculator
import numpy

kptlist = uc.getMonkhorstPackGrid((40,40,40)).reshape(64000,3)

sqecalc = AbInitio.kernelGenerator.SqeCalculator.SqeCalculator(uc, kpoints=kptlist)

#sqecalc.readIDFeigenvectors(filename='Polarizations_FeAl_8000k')
sqecalc.readIDFeigenvectors(filename='pols_FeAl222_mp40.idf')


# create some arbitrary energies for now:
#sqecalc._energies = numpy.zeros((8000,6), dtype='float')
#for ik in range(len(sqecalc._kpts)):
#    for imode in range(6):
#        sqecalc._energies[ik][imode] = numpy.dot(sqecalc._kpts[ik], sqecalc._kpts[ik])*(imode/3 + 1)

sqecalc.readEigenvaluesFromIDFomega2(filename='omega2_FeAl222_mp40.idf')

sqecalc._DebyeWallerCalculator._energies = sqecalc._energies
sqecalc._DebyeWallerCalculator._polvecs = sqecalc._polvecs

sqecalc._etransferTol = 0.5
sqecalc._qtransferTolRadius = 0.25

#loop
qstart = numpy.array([0.0, 0.0, 0.0])
deltaq = numpy.array([0.15, 0.0, 0.0])

estart = 0.0
deltae = 1.0

sqe = numpy.zeros((20,50), dtype='float')

for iq in range(20):
    for ie in range(50):
        qtransfer = qstart + iq * deltaq
        etransfer = estart + ie * deltae
        sqe[iq,ie] = sqecalc.calcSqeCohCreateAllmodes(qtransfer, etransfer)
        print iq, ie, sqe[iq,ie]

pylab.imshow(sqe)
pylab.show()

end = raw_input()

#########

#mp40 grid

mp40=uc.getMonkhorstPackGrid((40,40,40))
sqecalc._kpts= mp40
sqecalc._numkpts = 64000

sqecalc.readIDFeigenvectors(filename='pols_FeAl222_mp40.idf')
sqecalc.readEigenvaluesFromIDFomega2(filename='omega2_FeAl222_mp40.idf')

estart = 0.0
deltae = 0.5

sqegrid = numpy.zeros((40,40,40,40), dtype='float')

for kx in range(40):
    for ky in range(40):
        for kz in range(40):
            qtransfer = mp40[kx,ky,kz]
            

