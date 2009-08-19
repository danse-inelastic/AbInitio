from diffpy.Structure.lattice import Lattice
from math import sqrt, cos, radians
import numpy
from configParser import *


class QELattice(Lattice):
    def __init__(self, ibrav = 1,a = 1 ,b = 1,c = 1,cBC = 0.,cAC = 0. ,cAB = 0., fname = None ):
        Lattice.__init__(self)
        self.filename = fname
        self.ibrav = ibrav
        if self.filename != None:
            qeBaseTuple = self.getLatticeFromFile(fname)
        else:
            qeBaseTuple = self.getQEBaseFromParCos(ibrav ,a ,b ,c ,cBC ,cAC ,cAB )
        qeBase = numpy.array(qeBaseTuple[1], dtype = float)*qeBaseTuple[0]

        print 'Found "' + qeBaseTuple[2] + '" cell from PWSCF input'
        print 'Setting the base vectors according to QE conventions:'
        print qeBase
        self.setLatBase(qeBase)

        print self.abcABG()
#        print self.getQEBaseFromParCos(ibrav ,a ,b ,c ,cBC ,cAC ,cAB )[1]

    def getLatticeFromFile(self, fname):
        qeConf = QEConfig(fname)
        qeConf.parse()
        if 'ibrav' in qeConf.namelists['system'].params:
            self.ibrav  = int(qeConf.namelistParameter('system', 'ibrav'))
            if self.ibrav > 0:                
                a, b, c, cBC, cAC, cAB = self.getLatticeParamsFromFile(self.ibrav, fname)
                return self.getQEBaseFromParCos(self.ibrav, a, b, c, cBC, cAC, cAB)
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError

    def getLatticeParamsFromFile(self, ibrav, fname):
        qeConf = QEConfig(fname)
        qeConf.parse()
        cBC = 0.0
        cAC = 0.0
        cAB = 0.0
        a = float(qeConf.namelistParameter('system', 'celldm(1)'))
        if ibrav > 0 and ibrav < 4:
            return a, a, a, cBC, cAC, cAB
        if ibrav == 4 or ibrav == 6 or ibrav == 7:
            print "qweqweqwe"
            c_a = float(qeConf.namelistParameter('system', 'celldm(3)'))
            print c_a
            print a, a, c_a*a, cBC, cAC, cAB
            return a, a, c_a*a, cBC, cAC, cAB
        if ibrav == 5:
            cAB = float(qeConf.namelistParameter('system', 'celldm(4)'))
            return a, a, a, cBC, cAC, cAB
        if ibrav > 7 and ibrav < 12:
            b_a = float(qeConf.namelistParameter('system', 'celldm(2)'))
            c_a = float(qeConf.namelistParameter('system', 'celldm(3)'))
            return a, b_a*a, c_a*a, cBC, cAC, cAB
        if ibrav == 12 or ibrav == 13:
            b_a = float(qeConf.namelistParameter('system', 'celldm(2)'))
            c_a = float(qeConf.namelistParameter('system', 'celldm(3)'))
            cAB = float(qeConf.namelistParameter('system', 'celldm(4)'))
            return a, b_a*a, c_a*a, cBC, cAC, cAB
        if ibrav == 14:
            b_a = float(qeConf.namelistParameter('system', 'celldm(2)'))
            c_a = float(qeConf.namelistParameter('system', 'celldm(3)'))
            cBC = float(qeConf.namelistParameter('system', 'celldm(4)'))
            cAC = float(qeConf.namelistParameter('system', 'celldm(5)'))
            cAB = float(qeConf.namelistParameter('system', 'celldm(6)'))
            return a, b_a*a, c_a*a, cBC, cAC, cAB


    def getQEBaseFromParAngles(self, ibrav = 1, a = 1, b = 1, c = 1, alpha = 90.,beta = 90 ,gamma = 90):
        cBC = cos(radians(alpha))
        cAC = cos(radians(beta))
        cAB = cos(radians(gamma))
        return self.getQEBaseFromParCos(ibrav, a, b, c, cBC,cAC,cAB)



    def getQEBaseFromParCos( self, ibrav = 1, a = 1, b = 1, c = 1, cBC = 0.,cAC = 0. ,cAB = 0.):
        c_a = float(c)/a
        # description dictionary of QE base vectors:
        QEBase = {
        
        1 : (a, [[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]],
                 'Simple Cubic'),  # sc simple cubic
                 
        2 : (a/2., [[-1, 0, 1],
                    [0, 1, 1],
                    [-1, 1, 0]],
                    'Face Centered Cubic'),  # fcc face centered cubic

        3 : (a/2., [[1, 1, 1],
                    [-1, 1, 1],
                    [-1, -1, 1]],
                    'Body Centered Cubic'),  # bcc body entered cubic

        4 : (a, [[1, 0, 0],
                 [-0.5, sqrt(3.0)/2.0, 0.],
                 [0,    0,          c_a]],
                 'Simple Hexagonal or Trigonal(P)'),  # simple hexagonal and trigonal(p)

        5 : (a, [[sqrt((1.-cAB)/2.),-sqrt((1.-cAB)/6.), sqrt((1.+2.*cAB)/3.)],
                 [0, 2.*sqrt((1.-cAB)/6.),  sqrt((1.+2.*cAB)/3.)],
                 [-sqrt((1.-cAB)/2.), -sqrt((1.-cAB)/6.), sqrt((1.+2.*cAB)/3.)]],
                 'Trigonal(R)'),  # trigonal(r)

        6 : (a, [[1, 0, 0],
                 [0, 1, 0.],
                 [0, 0, c_a]],
                 'Simple Tetragonal(P)'),  # simple tetragonal (p)

        7 : (a/2., [[1, -1, c_a],
                    [1,  1, c_a],
                    [-1, -1, c_a]],
                    'Body Centered Tetragonal (I)'),  # body centered tetragonal (i)

        8 : (1.0, [[a, 0., 0.],
                    [0., b, 0.],
                    [0., 0., c]],
                    'Simple Orthorhombic (P)'),  # simple orthorhombic (p)

        9:  (1.0,  [[a/2., b/2., 0.],
                    [-a/2., b/2., 0.],
                    [0., 0., c]],
                    'Base Centered Orthorhombic'),  # bco base centered orthorhombic

        10: (1.0,  [[a/2., 0., c/2.],
                    [a/2., b/2., 0.],
                    [0., b/2., c/2.]],
                    'Face Centered Orthorhombic' ), # face centered orthorhombic

        11: (1.0,  [[a/2., b/2., c/2.],
                    [-a/2., b/2., c/2.],
                    [-a/2., -b/2., c/2.]],
                    'Body Centered Orthorhombic'),  # body centered orthorhombic

        12: (1.0,  [[a, 0, 0],
                    [b*cAB, b*sqrt(1.0 - cAB**2), 0],
                    [0, 0, c]],
                    'Monoclinic (P)'),  # monoclinic (p)

        13: (1.0,  [[a/2., 0, -c/2.],
                    [b*cAB, b*sqrt(1.0 - cAB**2), 0],
                    [a/2., 0, c/2.]],
                    'Base Centered Monoclinic'),  # base centered monoclinic

        14: (1.0,  [[a, 0, 0],
                    [b*cAB, b*sqrt(1.0 - cAB**2), 0],
                    [c*cAC, c*( cBC-cAC*cAB )/sqrt(1.-cAB**2), c*sqrt( 1. + 
                    2.*cBC*cAC*cAB - cBC**2 - cAC**2 - cAB**2)/sqrt(1.-cAB**2)]],
                    'Triclinic')  # triclinic
                    
        }
        
        return QEBase[ibrav]




qeLattice = QELattice(fname = 'ni.scf.in')

