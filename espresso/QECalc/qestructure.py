from diffpy.Structure.structure import Structure
from qelattice import QELattice
import numpy
from configParser import *


class QEStructure():
    def __init__(self, fname):
        """the structure is initialized from PWSCF config file"""
        self.filename = fname
        self.qeConf = QEConfig(fname)
        self.qeConf.parse()
        self.setStructureFromPWSCF()

    
    def setStructureFromPWSCF(self):
        """ Loads structure from PWSCF config file"""
        self.qeLattice = QELattice(fname = self.filename)
        self.structure = Structure(lattice = self.qeLattice.primitiveLattice)
        nat  = int(self.qeConf.namelist('system').param('nat'))
        ntyp  = int(self.qeConf.namelist('system').param('ntyp'))
        atomicLines = self.qeConf.card('atomic_positions').getLines()
        self.atomicPositionsType = self.qeConf.card('atomic_positions').argument()
        if self.atomicPositionsType == 'bohr' or self.atomicPositionsType == 'angstrom':
            raise NotImplementedError
        if self.atomicPositionsType == None:
            self.atomicPositionsType = 'alat'
        print self.atomicPositionsType
        for line in atomicLines:
            if '!' not in line:
                words = line.split()
                coordsConstraints = [float(w) for w in words[1:]]
                atomSymbol = words[0]
                if self.atomicPositionsType == 'alat':
                    coords = self.qeLattice.primitiveLattice.fractional(numpy.array(coordsConstraints[0:3])*self.qeLattice.a)
                if self.atomicPositionsType == 'crystal':
                    coords = numpy.array(coordsConstraints[0:3])
                if len(words) > 4:
                    constraint = [coordsConstraints[3], coordsConstraints[4],
                                  coordsConstraints[5]]
                self.structure.addNewAtom(atomSymbol, xyz = numpy.array(coords[0:3]))
#    def writeStructureToPWSCF(self):
#        """Writes structure into PWSCF config file"""
#        if self.qeLattice.latticeType == 'celldm'

        

        

if __name__ == '__main__':
    myStruct = QEStructure('zro2.scf.in')
    print myStruct.structure