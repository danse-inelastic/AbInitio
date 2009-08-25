from diffpy.Structure.lattice import Lattice, cosd
from math import sqrt, degrees
import numpy
from configParser import *


class QELattice():
    def __init__(self, ibrav = 1,a = 1 ,b = 1,c = 1,
                 cBC = 0.,cAC = 0. ,cAB = 0., fname = None ):
#        Lattice.__init__(self)
        self.filename = fname
        self.latticeType = 'celldm'
        self.qeConfig = None   # should be none if nothing to parse
        self.primitiveLattice = Lattice()
        self.standardLattice = Lattice()
        if self.filename != None:
            self.setLatticeFromPWSCF(self.filename)
        else:
            self.setLattice(ibrav ,a ,b , c, cBC ,cAC ,cAB )
#        print self.getQEBaseFromParCos(ibrav ,a ,b ,c ,cBC ,cAC ,cAB )[1]

    def setLatticeFromPrimitiveVectors(self, ibrav, vectors):
        """ Will extract conventional lattice parameters from primitive vectors.
            'vectors' - is a list with primitive vectors (in QE notation),
            including lattice parameters. For example from PWSCF output"""
        from numpy import dot
        # default values:
        a = b = c = 1.0
        cBC = cAC = cAB = 0.0
        v = numpy.array(vectors, dtype = float)
        if ibrav == 0:
            raise NotImplementedError
        # sc simple cubic:
        if ibrav == 1:
            a = v[0,0]

        if ibrav == 2:
            a = b = c = sqrt( 2.0*dot(v[0,:],v[0,:]) )

        if ibrav == 3:
            a = b = c = 2.0*sqrt( dot(v[0,:],v[0,:])/3.0)

        if ibrav == 4:
            a = b = sqrt( dot(v[0,:],v[0,:]))
            c = sqrt( dot(v[2,:],v[2,:]))
            cAB = cosd(120.0)

        if ibrav == 5:
            a = b = c = sqrt( dot(v[0,:],v[0,:]))
            cBC = cAC = cAB = dot(v[0,:],v[2,:])/a**2

        if ibrav == 6:
            a = b = sqrt( dot(v[0,:],v[0,:]))
            c = sqrt( dot(v[2,:],v[2,:]))

        if ibrav == 7:
            a = b = v[1,0] - v[2,0]
            c = v[1,2] + v[2,2]

        if ibrav == 8:
            a = v[0,0]
            b = v[1,1]
            c = v[2,2]

        if ibrav == 9:
            a = v[0,0] - v[1,0]
            b = v[0,1] + v[1,1]
            c = v[2,2]

        if ibrav == 10:
            a = v[2,0] - v[0,0] - v[1,0]
            b = v[2,1] - v[0,1] + v[1,1]
            c = v[0,2] - v[1,2] + v[2,2]

        if ibrav == 11:
            a = v[0,0] - v[1,0]
            b = v[1,1] - v[2,1]
            c = v[0,2] - v[2,2]

        if ibrav == 12:
            a = v[0,0]
            b = sqrt( dot(v[1,:],v[1,:]))
            cAB = v[1,0]/b
            c = v[2,2]

        if ibrav == 13:
            a = v[0,0] + v[2,0]
            b = sqrt( dot(v[1,:],v[1,:]))
            c = v[2,2] - v[0,2]
            cAB = v[1,0]/b

        if ibrav == 14:
            a = v[0,0]
            b = sqrt( dot(v[1,:],v[1,:]))
            cAB = v[1,0]/b
            c = sqrt( dot(v[2,:],v[2,:]))
            cAC = v[2,0]/c
            cBC = v[2,1]*sqrt(1.0 - cAB**2)/c + cAC*cAB

        self.setLattice(ibrav, a, b, c, cBC, cAC, cAB)

    def setLatticeFromPWSCF(self, fname = None):
        if fname != None:
            self.filename = fname
            self.qeConf = QEConfig(fname)
            self.qeConf.parse()
        if 'ibrav' in self.qeConf.namelists['system'].params:
            ibrav  = int(self.qeConf.namelist('system').param('ibrav'))
            if ibrav >= 0:
                a, b, c, cBC, cAC, cAB = self.getLatticeParamsFromPWSCF(
                                                              ibrav, fname)
            else:
                raise NotImplementedError("ibrav should be integer >= 0")
        else:
            raise NotImplementedError("config file should have ibrav defined")
        self.setLattice(ibrav, a, b, c, cBC, cAC, cAB)

    def setLattice(self, ibrav = None, a = None, b = None, c = None,
                   cBC = None, cAC = None, cAB = None):
        from math import acos
        if [ibrav, a,b,c,cBC,cAC,cAB] == 7*[None]:
            return None
        if ibrav is not None: self.ibrav = ibrav
        if a is not None: self.a = a
        if b is not None: self.b = b
        if c is not None: self.c = c
        if cBC is not None: self.cBC = cBC
        if cAC is not None: self.cAC = cAC
        if cAB is not None: self.cAB = cAB
        qeBaseTuple = self.getQEBaseFromParCos(self.ibrav, self.a, self.b,
                                           self.c, self.cBC, self.cAC, self.cAB)
        qeBase = numpy.array(qeBaseTuple[1], dtype = float)*qeBaseTuple[0]
        print 'Found "' + qeBaseTuple[2] + '" cell'
        print 'Setting the base vectors according to QE conventions:'
        print qeBase
        self.primitiveLattice.setLatBase(qeBase)
        alpha = degrees(acos(self.cBC))
        beta = degrees(acos(self.cAC))
        gamma = degrees(acos(self.cAB))
        self.standardLattice.setLatPar(self.a,self.b,self.c,alpha,beta,gamma)
        print "Standard Lattice:"
        print self.standardLattice.base

    def getLatticeParams(self):
        return [self.a, self.b,self.c, self.cBC, self.cAC, self.cAB]
        
    def getLatticeParamsFromPWSCF(self, ibrav, fname):
        qeConf = QEConfig(fname)
        qeConf.parse()
        cBC = 0.0
        cAC = 0.0
        cAB = 0.0
        if 'celldm(1)' in qeConf.namelists['system'].params:
            self.latticeType = 'celldm' # celldm(i), i=1,6
            a = float(qeConf.namelist('system').param('celldm(1)'))
            
            if ibrav == 0:
                # lattice is set in the units of celldm(1)
                # need to parse CELL_PARAMETERS
                cellParLines = qeConf.card('cell_parameters').getLines()
                cellParType = qeConf.card('cell_parameters').argument()
                if cellParType == 'cubic' or cellParType == None:
                    self.latticeType = 'generic cubic'
                else:
                    if cellParType == 'hexagonal':
                        self.latticeType = 'generic hexagonal'
            if ibrav > 0 and ibrav < 4:
                return a, a, a, cBC, cAC, cAB
            if ibrav == 4:
                cAB = cosd(120.0)
            if ibrav == 4 or ibrav == 6 or ibrav == 7:
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                return a, a, c_a*a, cBC, cAC, cAB
            if ibrav == 5:
                cAB = float(qeConf.namelist('system').param('celldm(4)'))
                return a, a, a, cAB, cAB, cAB
            if ibrav > 7 and ibrav < 12:
                b_a = float(qeConf.namelist('system').param('celldm(2)'))
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                return a, b_a*a, c_a*a, cBC, cAC, cAB
            if ibrav == 12 or ibrav == 13:
                b_a = float(qeConf.namelist('system').param('celldm(2)'))
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                cAB = float(qeConf.namelist('system').param('celldm(4)'))
                return a, b_a*a, c_a*a, cBC, cAC, cAB
            if ibrav == 14:
                b_a = float(qeConf.namelist('system').param('celldm(2)'))
                c_a = float(qeConf.namelist('system').param('celldm(3)'))
                cBC = float(qeConf.namelist('system').param('celldm(4)'))
                cAC = float(qeConf.namelist('system').param('celldm(5)'))
                cAB = float(qeConf.namelist('system').param('celldm(6)'))
                return a, b_a*a, c_a*a, cBC, cAC, cAB
        else:
            if ibrav == 0:
                print "Should specify celldm(1) if use 'generic' lattice"
                raise NotImplementedError
            a = float(qeConf.namelist('system').param('A'))
            self.latticeType = 'traditional'   # A, B, C, cosAB, cosAC, cosBC
            if ibrav > 0 and ibrav < 4:
                return a, a, a, cBC, cAC, cAB
            if ibrav == 4:
                cAB = cosd(120.0)
            if ibrav == 4 or ibrav == 6 or ibrav == 7:
                c = float(qeConf.namelist('system').param('C'))
                return a, a, c, cBC, cAC, cAB
            if ibrav == 5:
                cAB = float(qeConf.namelist('system').param('cosAB'))
                return a, a, a, cAB, cAB, cAB
            if ibrav > 7 and ibrav < 12:
                b = float(qeConf.namelist('system').param('B'))
                c = float(qeConf.namelist('system').param('C'))
                return a, b, c, cBC, cAC, cAB
            if ibrav == 12 or ibrav == 13:
                b = float(qeConf.namelist('system').param('B'))
                c = float(qeConf.namelist('system').param('C'))
                cAB = float(qeConf.namelist('system').param('cosAB'))
                return a, b, c, cBC, cAC, cAB
            if ibrav == 14:
                b = float(qeConf.namelist('system').param('B'))
                c = float(qeConf.namelist('system').param('C'))
                cBC = float(qeConf.namelist('system').param('cosBC'))
                cAC = float(qeConf.namelist('system').param('cosAC'))
                cAB = float(qeConf.namelist('system').param('cosAB'))
                return a, b, c, cBC, cAC, cAB

    def saveLatticeToPWSCF(self):
        if self.qeConf == None:
            raise NotImplementedError("writeLatticeToPWSCF: qeConf was not properly initialized")
        ibrav = self.ibrav
        if self.latticeType == 'celldm':
            self.qeConf.namelist('system').addParam('celldm(1)', self.a)
            self.qeConf.namelist('system').addParam('celldm(2)', self.b/self.a)
            self.qeConf.namelist('system').addParam('celldm(3)', self.c/self.a)
            if self.ibrav < 14:
                self.qeConf.namelist('system').addParam('celldm(4)', self.cAB)
            else:
                self.qeConf.namelist('system').addParam('celldm(4)', self.cBC)
                self.qeConf.namelist('system').addParam('celldm(5)', self.cAC)
                self.qeConf.namelist('system').addParam('celldm(6)', self.cAB)
        else:
            if self.latticeType == 'traditional':
                self.qeConf.namelist('system').addParam('A', self.a)
                self.qeConf.namelist('system').addParam('B', self.a)
                self.qeConf.namelist('system').addParam('C', self.a)
                self.qeConf.namelist('system').addParam('cosAB', self.cAB)
                self.qeConf.namelist('system').addParam('cosAC', self.cAC)
                self.qeConf.namelist('system').addParam('cosBC', self.cBC)

        self.qeConf.save(self.filename)

#        a = float(qeConf.namelistParameter('system', 'celldm(1)'))


    def getQEBaseFromParAngles(self, ibrav = 1, a = 1, b = 1, c = 1,
                               alpha = 90.,beta = 90 ,gamma = 90):
        cBC = cosd(alpha)
        cAC = cosd(beta)
        cAB = cosd(gamma)
        return self.getQEBaseFromParCos(ibrav, a, b, c, cBC,cAC,cAB)



    def getQEBaseFromParCos( self, ibrav = 1, a = 1, b = 1, c = 1,
                                    cBC = 0.,cAC = 0. ,cAB = 0.):
        c_a = float(c)/a
        # description dictionary of QE base vectors:
        QEBase = {
        # sc simple cubic:
        1 : (a, [[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]],
                 'Simple Cubic'),
        # fcc face centered cubic:
        2 : (a/2., [[-1, 0, 1],
                    [0, 1, 1],
                    [-1, 1, 0]],
                    'Face Centered Cubic'),
        # bcc body entered cubic:
        3 : (a/2., [[1, 1, 1],
                    [-1, 1, 1],
                    [-1, -1, 1]],
                    'Body Centered Cubic'), 
        # simple hexagonal and trigonal(p):
        4 : (a, [[1, 0, 0],
                 [-0.5, sqrt(3.0)/2.0, 0.],
                 [0,    0,          c_a]],
                 'Simple Hexagonal or Trigonal(P)'),
        # trigonal(r):
        5 : (a, [[sqrt((1.-cAB)/2.),-sqrt((1.-cAB)/6.), sqrt((1.+2.*cAB)/3.)],
                 [0, 2.*sqrt((1.-cAB)/6.),  sqrt((1.+2.*cAB)/3.)],
                 [-sqrt((1.-cAB)/2.), -sqrt((1.-cAB)/6.), sqrt((1.+2.*cAB)/3.)]],
                 'Trigonal(R)'),
        # simple tetragonal (p):
        6 : (a, [[1, 0, 0],
                 [0, 1, 0.],
                 [0, 0, c_a]],
                 'Simple Tetragonal(P)'),
        # body centered tetragonal (i):
        7 : (a/2., [[1, -1, c_a],
                    [1,  1, c_a],
                    [-1, -1, c_a]],
                    'Body Centered Tetragonal (I)'),
        # simple orthorhombic (p):
        8 : (1.0, [[a, 0., 0.],
                    [0., b, 0.],
                    [0., 0., c]],
                    'Simple Orthorhombic (P)'),  
        # bco base centered orthorhombic:
        9:  (1.0,  [[a/2., b/2., 0.],
                    [-a/2., b/2., 0.],
                    [0., 0., c]],
                    'Base Centered Orthorhombic'),
        # face centered orthorhombic:
        10: (1.0,  [[a/2., 0., c/2.],
                    [a/2., b/2., 0.],
                    [0., b/2., c/2.]],
                    'Face Centered Orthorhombic' ),
        # body centered orthorhombic:
        11: (1.0,  [[a/2., b/2., c/2.],
                    [-a/2., b/2., c/2.],
                    [-a/2., -b/2., c/2.]],
                    'Body Centered Orthorhombic'),
        # monoclinic (p):
        12: (1.0,  [[a, 0, 0],
                    [b*cAB, b*sqrt(1.0 - cAB**2), 0],
                    [0, 0, c]],
                    'Monoclinic (P)'),
        # base centered monoclinic:
        13: (1.0,  [[a/2., 0, -c/2.],
                    [b*cAB, b*sqrt(1.0 - cAB**2), 0],
                    [a/2., 0, c/2.]],
                    'Base Centered Monoclinic'),
        # triclinic:
        14: (1.0,  [[a, 0, 0],
                    [b*cAB, b*sqrt(1.0 - cAB**2), 0],
                    [c*cAC, c*( cBC-cAC*cAB )/sqrt(1.-cAB**2), c*sqrt( 1. + 
                    2.*cBC*cAC*cAB - cBC**2 - cAC**2 - cAB**2)/sqrt(1.-cAB**2)]],
                    'Triclinic')
                    
        }
        return QEBase[ibrav]


if __name__ == '__main__':

    qeLattice = QELattice(fname = 'zro2_2.scf.in')
    qeLattice.writeLatticeToPWSCF()
#    qeLattice2 = QELattice()
#    qeLattice2.setLatticeFromPrimitiveVectors(qeLattice.ibrav, qeLattice.primitiveLattice.base )
    #print qeLattice.primitiveLattice.base
    #testLattice = Lattice(5,5,5,90,90,90)
    #print testLattice.base
