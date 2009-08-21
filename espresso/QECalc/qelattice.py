from diffpy.Structure.lattice import Lattice, cosd
from math import sqrt, degrees
import numpy
from configParser import *


class QELattice():
    def __init__(self, ibrav = 1,a = 1 ,b = 1,c = 1,
                 cBC = 0.,cAC = 0. ,cAB = 0., fname = None ):
#        Lattice.__init__(self)
        self.filename = fname
        self.primitiveLattice = Lattice()
        self.standardLattice = Lattice()
        if self.filename != None:
            self.setLatticeFromFile(fname)
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
            a = b = c = 2.0*sqrt( dot(v[0,:],v[0,:])/3)

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

    def setLatticeFromFile(self, fname):
        qeConf = QEConfig(fname)
        qeConf.parse()
        if 'ibrav' in qeConf.namelists['system'].params:
            ibrav  = int(qeConf.namelistParameter('system', 'ibrav'))
            if ibrav > 0:                
                a, b, c, cBC, cAC, cAB = self.getLatticeParamsFromFile(
                                                              ibrav, fname)
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError
        self.setLattice(ibrav, a, b, c, cBC, cAC, cAB)

    def setLattice(self, ibrav, a, b, c, cBC, cAC, cAB):
        from math import acos
        self.ibrav = ibrav
        self.a = a
        self.b = b
        self.c = c
        self.cBC = cBC
        self.cAC = cAC
        self.cAB = cAB
        qeBaseTuple = self.getQEBaseFromParCos(ibrav, a, b, c, cBC, cAC, cAB)
        qeBase = numpy.array(qeBaseTuple[1], dtype = float)*qeBaseTuple[0]
        print 'Found "' + qeBaseTuple[2] + '" cell from PWSCF input'
        print 'Setting the base vectors according to QE conventions:'
        print qeBase
        self.primitiveLattice.setLatBase(qeBase)
        alpha = degrees(acos(cBC))
        beta = degrees(acos(cAC))
        gamma = degrees(acos(cAB))
        self.standardLattice.setLatPar(a,b,c,alpha,beta,gamma)
        print "Standard Lattice:"
        print self.standardLattice.base


    def getLatticeParamsFromFile(self, ibrav, fname):
        qeConf = QEConfig(fname)
        qeConf.parse()
        cBC = 0.0
        cAC = 0.0
        cAB = 0.0
        if 'celldm' in qeConf.namelists['system'].params:
            a = float(qeConf.namelistParameter('system', 'celldm(1)'))
            if ibrav > 0 and ibrav < 4:
                return a, a, a, cBC, cAC, cAB
            if ibrav == 4:
                cAB = cosd(120.0)
            if ibrav == 4 or ibrav == 6 or ibrav == 7:
                c_a = float(qeConf.namelistParameter('system', 'celldm(3)'))
                return a, a, c_a*a, cBC, cAC, cAB
            if ibrav == 5:
                cAB = float(qeConf.namelistParameter('system', 'celldm(4)'))
                return a, a, a, cAB, cAB, cAB
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
        else:
            a = float(qeConf.namelistParameter('system', 'A'))
            if ibrav > 0 and ibrav < 4:
                return a, a, a, cBC, cAC, cAB
            if ibrav == 4:
                cAB = cosd(120.0)
            if ibrav == 4 or ibrav == 6 or ibrav == 7:
                c = float(qeConf.namelistParameter('system', 'C'))
                return a, a, c, cBC, cAC, cAB
            if ibrav == 5:
                cAB = float(qeConf.namelistParameter('system', 'cosAB'))
                return a, a, a, cAB, cAB, cAB
            if ibrav > 7 and ibrav < 12:
                b = float(qeConf.namelistParameter('system', 'B'))
                c = float(qeConf.namelistParameter('system', 'C'))
                return a, b, c, cBC, cAC, cAB
            if ibrav == 12 or ibrav == 13:
                b = float(qeConf.namelistParameter('system', 'B'))
                c = float(qeConf.namelistParameter('system', 'C'))
                cAB = float(qeConf.namelistParameter('system', 'cosAB'))
                return a, b, c, cBC, cAC, cAB
            if ibrav == 14:
                b = float(qeConf.namelistParameter('system', 'B'))
                c = float(qeConf.namelistParameter('system', 'C'))
                cBC = float(qeConf.namelistParameter('system', 'cosBC'))
                cAC = float(qeConf.namelistParameter('system', 'cosAC'))
                cAB = float(qeConf.namelistParameter('system', 'cosAB'))
                return a, b, c, cBC, cAC, cAB


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




qeLattice = QELattice(fname = 'zro2.scf.in')
qeLattice2 = QELattice()
qeLattice2.setLatticeFromPrimitiveVectors(qeLattice.ibrav, qeLattice.primitiveLattice.base )
#print qeLattice.primitiveLattice.base
#testLattice = Lattice(5,5,5,90,90,90)
#print testLattice.base
