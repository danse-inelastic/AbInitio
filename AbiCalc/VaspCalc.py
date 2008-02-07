import os
import re
import weakref
import numpy as num

from ASE import Atom,ListOfAtoms

import AbInitio.vasp.parsing.parser2
import AbInitio.vasp.potcar
from AbInitio.vasp.parsing.SystemPM import *

from AbInitio.vasp.vasp import VASP
#from AbInitio.AbiCalc.AbiCalc import PlaneWaveAbiCalc,PeriodicAbiCalc

class VaspCalc(PlaneWaveAbiCalc):
    """A wrapper class for the VASP calculator."""

    def __init__(self, unitCell=None, kpts=None, ekincutoff=None,
                 name='vasp', xc='pawpbe', vaspcmd='vasp'):
        self._kpts = kpts
        self._ekincutoff = ekincutoff
        PeriodicAbiCalc.__init__(self, unitCell)
        self._vasp = VASP(name=name, kpts=kpts, pw=ekincutoff, xc=xc, vaspcmd=vaspcmd)


    def getPotEnergy(self):
        """Returns the potential energy for the ionic configuration."""
        return self._vasp.GetPotentialEnergy()


    def getForces(self):
        """Returns the forces on the nuclei for the current atomic configuration."""
        return self._vasp.GetCartesianForces()


    def getStress(self):
        """Returns the stress tensor on the unit cell in reduced coordinates."""
        return self._vasp.GetStress(self,units='reduced')


    pass # enf class vaspCalc
