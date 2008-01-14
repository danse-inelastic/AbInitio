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

kptlist = uc.getMonkhorstPackGrid((20,20,20)).reshape(8000,3)

sqecalc = AbInitio.kernelGenerator.SqeCalculator.SqeCalculator(uc, kpoints=kptlist)

sqecalc.readIDFeigenvectors(filename='Polarizations_FeAl_8000k')

# create some arbitrary energies for now:
sqecalc._energies = numpy.zeros((8000,6), dtype='float')
for ik in range(len(sqecalc._kpts)):
    for imode in range(6):
        sqecalc._energies[ik][imode] = numpy.dot(sqecalc._kpts[ik], sqecalc._kpts[ik])*(imode/3 + 1)

sqecalc._DebyeWallerCalculator._energies = sqecalc._energies
sqecalc._DebyeWallerCalculator._polvecs = sqecalc._polvecs

sqecalc._etransferTol = 3.0
sqecalc._qtransferTolRadius = 0.5

#loop
qstart = numpy.array([-3.0, 0.0, 0.0])
deltaq = numpy.array([0.6, 0.0, 0.0])

estart = 0.0
deltae = 5.0

sqe = numpy.zeros((10,10), dtype='float')

for iq in range(10):
    for ie in range(10):
        qtransfer = qstart + iq * deltaq
        etransfer = estart + ie * deltae
        sqe[iq,ie] = sqecalc.calcSqeCohCreateAllmodes(qtransfer, etransfer)
        print sqe[iq,ie]

pylab.imshow(sqe)
pylab.show()

end = raw_input()



