from crystal import UnitCell as UC
from crystal.UnitCell import *
from crystal.Atom import Atom

uc = UnitCell()
#vectors = [(4.70898,0,0), (0,4.70898,0), (0,0,4.70898)]
# better to use fractional coords for now (for unitCell <> loa translations...)
vectors = [(1,0,0), (0,1,0), (0,0,1)]
uc.setCellVectors(vectors)

at1 = Atom(symbol='Si')
pos1 = (0,0,0)
at2 = Atom(symbol='Si')
pos2 = (0.5,0.5,0.5)
at3 = Atom(symbol='V')
pos3 = (0.25, 0.50, 0.00)
at4 = Atom(symbol='V')
pos4 = (0.75, 0.50, 0.00)
at5 = Atom(symbol='V')
pos5 = (0.00, 0.25, 0.50)
at6 = Atom(symbol='V')
pos6 = (0.00, 0.75, 0.00)
at7 = Atom(symbol='V')
pos7 = (0.50, 0.00, 0.25)
at8 = Atom(symbol='V')
pos8 = (0.50, 0.00, 0.75)

site1 = Site(pos1, at1)
site2 = Site(pos2, at2)
site3 = Site(pos3, at3)
site4 = Site(pos4, at4)
site5 = Site(pos5, at5)
site6 = Site(pos6, at6)
site7 = Site(pos7, at7)
site8 = Site(pos8, at8)

sites = [site1, site2, site3, site4, site5, site6, site7, site8]

for i in range(8):
    uc.addSite(sites[i], sites[i].getAtom().symbol+str(i))


from AbInitio.AbiPhon.parsing.phonParsers import parseQpoints,parseFastPhon2IDF,plotQpoints

qpts, wts = parseQpoints()

#parseFastPhon2IDF(inputfilename='phon.out_mp20',
#                  polarizationsfile='polarizations.idf',
#                  omega2sfile='omega2s.idf',
#                  D=3)

#from inelastic.idf.Polarizations import read as readPols
#polsdata = readPols(filename="polarizations.idf")
#pols = polsdata[1]

#from inelastic.idf.Omega2 import read as readOm2s
#om2sdata = readOm2s(filename='omega2s.idf')
#om2s = om2sdata[1]

# S(Q,E) calculator:

import AbInitio.kernelGenerator.SqeCalculator

sqecalc = AbInitio.kernelGenerator.SqeCalculator.SqeCalculator(uc,kpoints=qpts)

sqecalc._qtransferTolRadius = 0.25
sqecalc._etransferTol = 1.0

sqecalc.readIDFeigenvectors(filename='polarizations.idf')
sqecalc.readEigenvaluesFromIDFomega2(filename='omega2s.idf')

sqecalc._DebyeWallerCalculator.getDWFactorForAtom(0, [0.5,0,0], 300.0)
