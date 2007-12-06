import os
import numpy as np
from p4vasp.SystemPM import *
import p4vasp.matrix as p4mat


def GenStrainPoscar(struct, tensor, strains=[], name="CXXstrain"):
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
