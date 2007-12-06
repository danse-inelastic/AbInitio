import os
import numpy as np
import scipy.linalg as la

try:
    from p4vasp.SystemPM import *
    from p4vasp.Structure import Structure
    import p4vasp.matrix as p4mat
except ImportError:
    print "P4Vasp could not be imported in Python."

from UnitCell import *
import CubicStandardStrainTensors as cubic


class StrainGenerator():
    """Component to generate strained unit cells."""

    def __init__(unitcell=None, tensor=None, strains=[0.01]):

        if unitcell is None: unitcell = UnitCell()
        if tensor is Note: tensor = cubic.C11

        self._uc = unitcell
        self._tensor = tensor
        self._strains = strains

        return

    def getUnitCell(self):
        return self._uc

    def setUnitCell(self, unitcell):
        self._uc = unitcell
        return

    def getStrainTensor(self):
        return self._tensor

    def setStrainTensor(self, tensor):
        self._tensor = tensor
        return

    def getStrains(self):
        return self._strains

    def setStrains(self,strains):
        self._strains = strains
        return

    def genStrainedPoscar(self,name="CXXstrain", struct=None, tensor=None):
        """Produces a sequence of strained P4vasp structures, \n\
        for strain values in the NumPy array strains."""

        # Parse the cell volume and basis vectors for the structure
        # with the P4vasp parser for the vasprun.xml 
        try:
            V0=struct.getCellVolume()
            print "Original volume=", V0
        except:
            print "Could not access structure."

        try:
            basis=struct.basis
        except:
            print "Could not access structure."
        
        b0=basis[0]
        b1=basis[1]
        b2=basis[2]

        # Generate the strained structures according to the strain sequence
        # and the strain tensor
        for x in strains:
            scale=p4mat.Matrix([[x, 0.0, 0.0],[0.0, x, 0.0],[0.0, 0.0, x]])
            try:
                tmp1=tensor.__rmul__(scale)
            except:
                print "Tensor is not a valid P4vasp matrix."
        
            tmp2=p4mat.Matrix(3,3)
            tmp2.identity()
            tmp=tmp1.__add__(tmp2)

            print tmp

            # Apply the tensor deformation to the unitcell vectors
            struct.basis[0]=tmp.__mul__(b0)
            print "b0=", struct.basis[0]
            struct.basis[1]=tmp.__mul__(b1)
            print "b1=", struct.basis[1]
            struct.basis[2]=tmp.__mul__(b2)
            print "b2=", struct.basis[2]

            # The transformation needs to be volume-conserving,
            # so we need to rescale the volume 
            V=struct.getCellVolume()
            print "New volume=", V
            y=float(np.power((V0/V), 1/3.))
            struct.scaleBasis(y,y,y)

            # Now write the strained structures to file
            # in VASP POSCAR format using the P4vasp parsing
            struct.write("pos_"+name+'_'+repr(x), newformat=0)
        
            # Reset to unstrained basis
            struct.basis[0]=b0
            struct.basis[1]=b1
            struct.basis[2]=b2

        continue # end of loop on strains
        return # end of genStrainedPoscar

    def getNormStrainUnitCell(self, strain=0.01, tensor=None):
        """Returns a unit cell object, obtained by applying the strain tensor,
        with specified amplitude, to the component unit cell.
        The unit cell is 'normalized' to conserve volume."""

        if tensor is None: tensor = self._tensor

        volume = self._uc.getVolume()

        straintensor = np.array(tensor) * (1.0 + strain) # array element-wise multiply
        cellvecs = np.array(self._uc.getCellVectors())
        newvecs = np.dot(straintensor, cellvectors) # matrix-like multiply on arrays

        # normalize the volume
        volume = self._uc.getVolume()
        newvol = abs(la.det(newvecs))
        cuberescale = np.power(volume / newvol, 1./3.)
        newvecs = newvecs * cuberescale #<< should check that this is what we want - O.D. 06/07

        atoms = self._uc.getAtoms()
        positions = self._uv.getPositions()
        uc = create_unitcell(newvecs, atoms, positions)
        
        return uc

    pass # end of class StrainGenerator
