__doc__ = """A module to implement a common interface to first-principles calculation engines."""

class FirstPrincipCalc:
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

        
    pass # end of class FirstPrincipCalc



### We need a class to handle the case of first-principles engines working with periodic systems:

class PeriodicFirstPrincipCalc(FirstPrincipCalc):
    """A first principles calculation interface for periodic systems and calculation engines."""

    def __init__(self, unitCell=None):
        FirstPrincipCalc.__init__(unitCell)

    def getStress(self):
        """Returns the stresses on the unit cell."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    pass # end of class PeriodicFirstPrincipCalc


### We need separate classes to handle plane-wave type calculations vs. localized absis-set calculations:


class PlaneWaveFirstPrincipCalc(PeriodicFirstPrincipCalc):
    """A class defining an interface to plane-waves, periodic first-principles calcultion engines."""

    def __init__(self, unitCell=None, kpts=None, ekincutoff=None):
        self._kpts = kpts
        self._ekincutoff = ekincutoff
        PeriodicFirstPrincipCalc.__init__(unitCell)

    def getBandEnergies(self):
        """Returns all the band energies, for all the k-points."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    def getBandEnergiesAtPoint(self, kpt):
        """Returns the band energies at a specific k-point."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    pass # end of class PlaneWaveFirstPrincipCalc


class BasisSetFirstPrincipCalc(PeriodicFirstPrincipCalc):
    """A class defining an interface to local basis-set, periodic first-principles calcultion engines."""

    def __init__(self, unitCell=None, basisSet=None):
        self._basisSet = basisSet
        PeriodicFirstPrincipCalc.__init__(unitCell)

    def getBandEnergies(self):
        """Returns all the band energies, for all the k-points."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    def getBandEnergiesAtPont(self, kpt):
        """Returns the band energies at a specific k-point."""
        # this should only be implemented in the derived classes
        raise NotImplementedError

    pass # end of class BasisSetFirstPrincipCalc
    
