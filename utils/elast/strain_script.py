#! /usr/bin/env python


__doc__ = """The purpose of this script is to integrate generation of strained structures and the calculation of energies with VASP, and then retrieve the energy-vs-strain for fitting."""

import os
import sys
import numpy as np
import scipy
import scipy.io
from CalcStrains import CalcStrains
from GenStrainPoscar import GenStrainPoscar
from FitCubicStrain import FitCubicStrain

try:
    from p4vasp.SystemPM import *
    from p4vasp.Structure import Structure
    import p4vasp.matrix as p4mat
except ImportError:
    print "P4Vasp could not be imported in Python."

def CalculateElasticConstant(straintype):
    """Calculate the elastic constant corresponding to straintype
    using ab-initio energy calculation (VASP).
    For a C11 deformation, one can use the deformation tensor:
    [[1  0  0]
     [0  0  0]
     [0  0  0]]
    
    For a C11-C12 deformation, one can use the deformation tensor:
    [[1  0  0]
     [0 -1  0]
     [0  0  0]]
    ( 8 symmetry operators)
    or:
    [[1  0  0]
     [0  1  0]
     [0  0 -2]]
    (16 symmetry operators)
    
    For a C44 deformation, one can use the deformation tensor:
    [[0  0  0]
     [0  0  1]
     [0  1  0]]
    (4 symmetry operators)
    or:
    [[0  1  1]
     [1  0  1]
     [1  1  0]]
    (12 symmetry operators)
    All these tensors are now listed in the module CubicStandardStrainTensors.
    """
    
    # load a structure from existing POSCAR file:
    struct = Structure("POSCAR")
    cellvolume = struct.getCellVolume()

    # define a sequence of strain amplitudes:
    strains=scipy.arange(-0.020, 0.022, 0.002).round(5)

    # For a C11 deformation, one can use the deformation tensor:
    # [[1  0  0]
    #  [0  0  0]
    #  [0  0  0]]
    #
    # For a C11-C12 deformation, one can use the deformation tensor:
    # [[1  0  0]
    #  [0 -1  0]
    #  [0  0  0]]
    # ( 8 symmetry operators)
    # or:
    # [[1  0  0]
    #  [0  1  0]
    #  [0  0 -2]]
    # (16 symmetry operators)
    #
    # For a C44 deformation, one can use the deformation tensor:
    # [[0  0  0]
    #  [0  0  1]
    #  [0  1  0]]
    # (4 symmetry operators)
    # or:
    # [[0  1  1]
    #  [1  0  1]
    #  [1  1  0]]
    # (12 symmetry operators)
    # all these tensors are now listed in the module CubicStandardStrainTensors

    import CubicStandardStrainTensors as cubic
    # straintype = 'C44_2'
    tens = cubic.__dict__[straintype]
    # tens = cubic.__dict__['C44_2']
    # tensor = p4mat.Matrix([[1.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]])
    tensor = p4mat.Matrix( tens.tolist() )

    # cmd='rsh -. n00 mpirun  -np 8 vasp'
    cmd = 'mr vasp'

    # Generate the strained POSCAR files
    GenStrainPoscar(struct, tensor, strains, straintype)

    # Execute the VASP command on the sequence of POSCAR files
    CalcStrains(strains, straintype, cmd)

    print "Finished calculating strains."

    # Collect the energies from the output vasprun.xml files.
    runs=[]
    energies=[]
    forces=[]

    for x in strains:
        runs.append(XMLSystemPM('vasprun_'+straintype+'_'+repr(x)+'.xml'))
        energies.append(runs[-1].FREE_ENERGY)
    pass

    for i in range(len(runs)):
        print runs[i].FINAL_STRUCTURE.getCellVolume()

    # for i in range(len(runs)):
    #    forces.append(runs[i].FORCES_SEQUENCE)
    pass

    evec=np.array(energies)
    data=scipy.zeros((len(evec),2),  dtype='f')

    # We build the data to outfile to file, which is used later for fitting
    # the strains are saved to file in Angstroems, not fractional coordinates
    # if the structure read from the POSCAR file is "direct", we need to apply scaling:
    # if struct.isDirect():
    #     scaling = struct.scaling[0]
    #    data[:,0] = strains * scaling
    # else the strains are already in units of length (assumed Angstroems):
    # elif struct.isCartesian():
    #  data[:,0] = strains

    data[:,0] = strains
    data[:,1] = evec
                
    datafile=open(straintype+"_data.dat", 'w')
    scipy.io.write_array(datafile, data)
    datafile.close()

    elastConst = FitCubicStrain(straintype+"_data.dat", straintype, cellvolume)
    print straintype, ": ", elastConst
    print straintype, " (GPa) : ", elastConst*160.22

    print "Finished all."

    pass # End of CalculateElasticConstant

if __name__ == "__main__":
	try:
		straintype = sys.argv[1]
		CalculateElasticConstant(straintype)
	except:
		print "Usage:", sys.argv[0], "strainstypestring (see CubicStandardStrainTensors)"

# End of file
