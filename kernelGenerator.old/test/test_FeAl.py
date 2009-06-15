from UnitCell import *
from Atom import *

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

qtransfer = sqecalc._kpts[4248]
sqecalc._qtransfer = qtransfer

etransfer = 2.0
sqecalc._etransfer = etransfer

sqecalc._DebyeWallerCalculator._energies = sqecalc._energies
sqecalc._DebyeWallerCalculator._polvecs = sqecalc._polvecs


sqecalc._etransferTol = 0.5
sqecalc._qtransferTolRadius = 0.5

print sqecalc.calcSqeCohCreateAllmodes(qtransfer, etransfer)


