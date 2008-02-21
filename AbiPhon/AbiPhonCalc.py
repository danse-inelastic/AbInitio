# Olivier Delaire


__doc__ = """A module that defines an interface for a first-principles phonon calculator."""


class AbiPhonCalc:
    """A phonon calculator based on a first-principles calculation."""

    def __init__(self, unitcell, supersize=[1,1,1], qpts=None, abiCalc=None):
        self._unitcell = unitcell
        self._supersize = supersize
        self._qpts = qpts
        self._weights = None
        self._abiCalc = abiCalc
        self._superCell = None
        self._superCellReady = False
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

    def setQpts(self, qpts, weights):
        """Set the Q-points at which the phonons are calculated,
        and their weights, in case of a symmetry-reduced list of Q-points."""
        self._qpts = qpts
        self._weights = weights

    def getQpts(self):
        """Returns the Q-points and theit weights."""
        return (self._qpts, self._weights)

    def generateSupercell(self, supersize=None):
        """Generates a supercell from the crystal unit cell and the supercell size."""
        # calculate the supercell based on the unit cell and the supercell size
        raise NotImplementedError

    def calcDisplacements(self):
        """Calculates atomic displacements,
        from which the full interatomic force-constants tensor
        can be obtained."""
        raise NotImplementedError

    def setDisplacements(self, displacements):
        """Sets the displacements for which to calculate the forces."""
        raise NotImplementedError

    def getDisplacements(self):
        """Returns the displacements for which are calculated."""
        raise NotImplementedError

    def calcForces(self):
        """Calculates the forces on the ions for all configurations,
        using the first-principles calculator."""
        raise NotImplementedError


    def _calcForcesForConfig(self, config):
        """helper function to calculate the forces on a particular
        atomic configuration."""
        raise NotImplementedError

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
        raise NotImplementedError

    def writePolVecs(self, filename):
        """Writes the phonon polarization vectors to file in 'idf' format."""
        raise NotImplementedError

    def writeEnergies(self, filename):
        """Writes the phonon energies to file."""
        raise NotImplementedError

    def writeQpts(self, filename):
        """Writes the phonon Q-points to file."""
        raise NotImplementedError

    pass # End of class AbiPhononCalc
