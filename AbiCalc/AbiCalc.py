__doc__ = """A module to implement a common interface to first-principles calculation engines."""

class AbiCalc:
    """Abstract class for interface to first principle calculator."""

    def __init__(self, atomicConfig=None):
        self._atomicConfig = atomicConfig


    def setAtomicConfig(self, atomicConfig):
        self.atomicConfig = atomicConfig

    def getAtomicConfig(self):
        return self.atomicConfig

    def calculate(self):
        """Launches the calculation engine."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    def getPotEnergy(self):
        """Returns the potential energy for the current atomic configuration."""        # this should only be implemented in the derived classes
        raise NotImplementedError

    def getForces(self):
        """Returns the forces on the nuclei for the current atomic configuration."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    def optimizeGeometry(self):
        """optimizes the geometry of the atomic configuraiton to minimize the forces on the nuclei."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

        
    pass # end of class AbiCalc



### We need a class to handle the case of first-principles engines working with periodic systems:

class PeriodicAbiCalc(AbiCalc):
    """A first principles calculation interface for periodic systems and calculation engines."""

    def __init__(self, unitCell=None):
        AbiCalc.__init__(self, unitCell)

    def getStress(self):
        """Returns the stresses on the unit cell."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    pass # end of class PeriodicAbiCalc


### We need separate classes to handle plane-wave type calculations vs. localized absis-set calculations:


class PlaneWaveAbiCalc(PeriodicAbiCalc):
    """A class defining an interface to plane-waves, periodic first-principles calcultion engines."""

    def __init__(self, unitCell=None, kpts=None, ekincutoff=None):
        self._kpts = kpts
        self._ekincutoff = ekincutoff
        PeriodicAbiCalc.__init__(self, unitCell)

    def getBandEnergies(self):
        """Returns all the band energies, for all the k-points."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    def getBandEnergiesAtPoint(self, kpt):
        """Returns the band energies at a specific k-point."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    pass # end of class PlaneWaveAbiCalc


class BasisSetAbiCalc(PeriodicAbiCalc):
    """A class defining an interface to local basis-set, periodic first-principles calculation engines."""

    def __init__(self, unitCell=None, basisSet=None):
        self._basisSet = basisSet
        PeriodicAbiCalc.__init__(self, unitCell)

    def getBandEnergies(self):
        """Returns all the band energies, for all the k-points."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    def getBandEnergiesAtPoint(self, kpt):
        """Returns the band energies at a specific k-point."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    pass # end of class BasisSetAbiCalc
    
