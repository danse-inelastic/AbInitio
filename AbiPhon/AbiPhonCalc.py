# Olivier Delaire


__doc__ = """A module that defines an interface for a first-principles phonon calculator."""

import numpy as np
from crystal.UnitCell import Site, UnitCell
from crystal.Atom import Atom


class AbiPhonCalc:
    """A phonon calculator based on a first-principles calculation."""

    def __init__(self, unitcell=None, supersize=[1,1,1], abiCalc=None, qpts=[8,8,8]):
        self._unitcell = unitcell
        self._abicalc = abicalc
        self._supersize = supersize
        self._qpts = qpts
        self._weights = None
        self._supercell = None
        self._supercellReady = False
        self._amplitude = 0.01  #amplitude of displacement, fractional

        self._energies = None # phonon energies
        self._polVecs = None  # phonon polarization vectors
        
        pass # end of __init__

    def setUnitCell(self, unitcell):
        """sets the unit cell."""
        self._unitcell = unitcell

    def getUnitCell(self):
        """returns the unit cell."""
        return self._unitcell

    def setAbiCalc(self, abicalc):
        """sets the ab-initio calculation engine."""
        self._abicalc = abicalc

    def setSuperSize(self, supersize=[1,1,1]):
        """sets the supercell size in each dimension, eg:
        phoncalc.setSupercell (supercell=[2,2,2]),
        to calculate a 2x2x2 supercell."""
        self._supersize = supersize
        
    def getSuperSize(self):
        """returns the size of the supercell."""
        return self._supersize

    def setAmplitude(self, amplitude):
        """Sets the atomic displacement amplitude."""
        self._amplitude = amplitude

    def getAmplitude(self, amplitude):
        """Returns the displacement amplitude."""
        return self._amplitude

    def setQptsWts(self, qpts, weights):
        """Set the Q-points at which the phonons are calculated,
        and their weights, in case of a symmetry-reduced list of Q-points."""
        self._qpts = qpts
        self._weights = weights

    def getQptsWts(self):
        """Returns the Q-points and theit weights."""
        return (self._qpts, self._weights)

    def genSupercell(self, superdims=None):
        """Generates a supercell of dimensions given by superdims."""
        if superdims is None:
            superdims = self._supersize
        if superdims is None:
            raise ValueError, 'supercell should be integer multiple of unit cell.'
        # generate a supercell with multiplied lattice vectors:
        supercell = UnitCell()
        cellvectors = self._unitcell.getCellVectors()
        supercellvectors = cellvectors * np.array(superdims)
        print supercellvectors
        # sa1 = a1 * dim1; sa2 = a2 * dim2; sa3 = a3 * dim3
        supercell.setCellVectors(supercellvectors)
        # Add the images of all the atoms:
        for i0 in range(superdims[0]):
            for i1 in range(superdims[1]):
                for i2 in range(superdims[2]):
                    for site in self._unitcell:
                        pos = site.getPosition()
                        cart = self._unitcell.fractionalToCartesian(pos)
                        newcart = (cart
                                   + i0 * cellvectors[0]
                                   + i1 * cellvectors[1]
                                   + i2 * cellvectors[2])
                        newpos = supercell.cartesianToFractional(newcart)
                        newsite = Site(newpos, Atom(Z=site.getAtom().Z))
                        supercell.addSite(newsite, '')
        self._supercell = supercell
        self._supercellReady = True
        return

    def calcDisplacements(self):
        """Calculates atomic displacements,
        from which the full interatomic force-constants tensor
        can be obtained."""
        raise NotImplementedError

    def setDisplacements(self, displacements):
        """Sets the displacements for which to calculate the forces.
        It expects a list of displacements [ disp0, disp1, ..., dispN],
        where dispN = (atnumN, (dxN, dyN, dzN)),
        with atnumN the number of the atom to be displaced,
        and (dxN, dyN, dzN) the displacement vector in fractional coords.        
        """
        self._disp = displacements
        return

    def getDisplacements(self):
        """Returns the displacements for which are calculated."""
        return self._disp

    def calcForces(self):
        """Calculates the forces on the ions for all configurations,
        using the first-principles calculator."""
        raise NotImplementedError


    def _calcForcesForDisp(self, dispNum):
        """helper function to calculate the forces on a particular
        atomic displacement number, dispN.
        It expects dispN between 0 and len(self._disp)."""
        if self._disp is None:
            raise ValueError, 'No displacements found in self._disp.'
        if not self._supercellReady:
            raise ValueError, 'Supercell is not generated.'
        try:
            disp = self._disp[dispNum]
            print disp
        except:
            raise ValueError, 'Invalid displacement.'
        distorted = self._supercell.__copy__()
        atnum = disp[0]
        dispvec = np.array(disp[1])
        oldvec = distorted._sites[atnum].getPosition()
        newvec = np.array(oldvec) + np.array(dispvec)
        distorted._sites[atnum].setPosition(newvec)
        print "Distorted supercell: ", distorted
        
        self._abicalc.setUnitCell(distorted)
        forces = self._abicalc.getForces()
        return forces

    def calcPhonons(self):
        """Calculate the phonons from the interatomic force-constants
        obtained from forces on ions in displaced configurations.
        This calculates the phonon energies and polarization vectors,
        at the q-points specified in the calculator."""
        raise NotImplementedError

    def getPolVecs(self):
        """Gets the phonon polarization vectors."""
        raise NotImplementedError

    def getEnergies(self):
        """Gets the phonon energies."""
        return self._energies

    def readFracQptsWtsIDF(self,filename='FractionalQs.idf'):
        """read fractional q-points coordinates from IDF format file"""
        from idf.FractionalQs import read
        data = read(filename=filename)
        self._qpts = data[1]
        self._weights = data[2]
        return
        
 
    def writePolVecsIDF(self, filename='Polarizations.idf', comment=''):
        """Writes the phonon polarization vectors to file in 'idf' format."""
        from idf.Polarizations import write
        write(self._polvecs, filename=filename, comment=comment)
        return

    def writeOmega2IDF(self, filename='Omega2.idf', comment=''):
        """Writes the phonon frequencies squared to file in 'idf' format."""
        from idf.Omega2 import write
        write(self._energies, filename=filename, comment=comment)
        return

    def writeFracQptsWtsIDF(self, filename='FractionalQs.py', comment=''):
        """Writes the phonon Q-points to file."""
        from idf.FractionalQs import write
        write(self._qpts, self._weights, filename=filename, comment=comment)
        return

    pass # End of class AbiPhononCalc
