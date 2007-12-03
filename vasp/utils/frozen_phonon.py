#!/usr/bin/python

__doc__ = """The purpose of this script is to calculate energies for a series of supercells deformed according to symmetry vectors produced by FROZSL."""

import os
import sys
import numpy as np
import scipy
import scipy.io

try:
    from p4vasp.SystemPM import *
    from p4vasp.Structure import Structure
    import p4vasp.matrix as p4mat
except ImportError:
    print "P4Vasp could not be imported in Python."


def isNearlyZero(x, eps):
    """Evaluates wether some number is within eps of zero."""
    return (abs(x) < eps)

def vectorEqual(x,y, eps):
    """check wether two vectors (hashable containers of numbers) are equal."""
    if len(x) != len(y): return False
    res = True
    for (xi, yi) in zip(x,y):
        res = res & isNearlyZero(xi-yi, eps)
        # this could be made more efficient'
        # since they are unequal as soon as a pair of components are not equal
    return res

# name of distortion
name = "X3-_sym1"

# load a structure from existing POSCAR file:
struct = Structure("POSCAR_X3-")
print "Loaded POSCAR."
cellvolume = struct.getCellVolume()
print "Cell volume: ", cellvolume

# backup the original structure
#struct.write("POSCAR_X3-_ref", newformat=0)

# open output file:
outfile = open('frozen_phonon_'+name+'.dat', 'w')


# load the symmetry vector
try:
    symfile = open("symmetryvector.in", 'r')
except: IOError, "symmetryvector.in file not found"

# Load the symmetry vector (here for 2x2x2 BCC supercell):
symveclist = []
for line in symfile:
    symveclist.append([float(s) for s in line.split()])

symvec = np.array(symveclist)

# load the reference atomic positions as output from FROZSL
# (needed because FROZSL may have shuffled atomic positions order from POSCAR)
# These must be the cartesian atomic positions, as in the POSCAR file.
# (!) the 'dimensionless' coordinates in frozsl_init output are NOT fractional coords (!)
try:
    posfile = open("positionvector.in", 'r')
except: IOError, "positionvector.in file not found"

# Load the symmetry vector (here for 2x2x2 BCC supercell):
posveclist = []
for line in posfile:
    posveclist.append([float(s) for s in line.split()])

posvec = np.array(posveclist)


# define a sequence of strain amplitudes:
eps = 1e-6
strains = scipy.arange(-0.4, 0.42, 0.02).round(5)
    
#strains = np.arange(-0.1, 0.15, 0.05).round(3)
#strains = np.arange(0.0, 0.03, 0.01).round(3)
for i in range(len(strains)):
    if isNearlyZero(strains[i], eps): strains[i] = 0.0

print "Strains: ", strains

# Loop on all the strain amplitudes:

runflag = True

for strain in strains:
    print "Strain: ", strain
    # Now we need to add the symmetry vector * the strain.
    # For this, we need to identify which struct positions correspond
    # to which position in the posvec output by FROZSL,
    # and apply the corresponding symvec * strain.

    struct = Structure("POSCAR_X3-")
    
    for refposnum in range(len(struct.positions)):
        for posnum in range(len(posvec)):
            if vectorEqual(posvec[posnum], struct.positions[refposnum], 1e-4):
                print "Found the position to move."
                struct.positions[refposnum] = struct.positions[refposnum] + p4mat.Vector(symvec[posnum].tolist()) * strain
    print "New positions: "
    print struct.positions

    struct.write("POSCAR_"+name+"_u_"+repr(strain), newformat=0)

    if runflag:
        # run the calculation
        vaspcmd = 'mr vasp'

        os.system("cp POSCAR_"+name+"_u_"+repr(strain)+" POSCAR")
        os.system("more POSCAR")
        os.system(vaspcmd)
        os.system('cp vasprun.xml vasprun_'+name+'_u_'+repr(strain)+'.xml')
        os.system('cp OUTCAR OUTCAR_'+name+'_u_'+repr(strain))
        os.system('cp CHGCAR CHGCAR_'+name+'_u_'+repr(strain))
        currentrun = XMLSystemPM('vasprun.xml')
        currentenergy = currentrun.FREE_ENERGY
        os.system(" E=`tail -n 2 OSZICAR` ; echo " + ("%1.8g" % strain) + " $E  >> frozen_phonon.out")
        outfile.write(repr(strain)+' '+repr(currentenergy)+'\n')

outfile.close()


def writeXYZ(name, strains, outfilename):
    """Reads the output XML files from VASP computation of frozen_phonons
    and creates a XYZ format file for input into GULP."""

    from p4vasp.SystemPM import *
    from p4vasp.Structure import Structure

    outfile = open(outfilename, 'a')  # open the file in 'append mode'
    weight = 100.0 # this should be different weights for different configurations

    for x in strains:
        xmlfilestring = 'vasprun_'+name+'_u_'+repr(strain)+'.xml'
        poscarstring = 'POSCAR_'+name+'_u_'+repr(strain)
        xmlrun = XMLSystemPM(xmlfilestring)
        energy = xmlrun.FREE_ENERGY
        struct = xmlrun.FINAL_STRUCTURE
        # set the atomic positions to fractional coordinates:
        struct.setDirect()
        # write a description line
        outfile.write('# ' + name + ' ' + repr(strain) + '\n')
        # write the cell vectors
        outfile.write('vectors \n')
        for vec in struct.basis:
            outfile.write(str(vec) + '\n')
        # Write the atoms symbols and atomic positions:
        outfile.write('fractional ' + len(struct.positions) + '\n')
        for i in range(len(struct.positions)):
            atomstr = struct.getRecordForAtom(i)
            vec = str(struct.positions[i])
            outfile.write(atomstr + ' ' + str(vec) + '\n')

        outfile.write('# energy for this configuration: \n observable \n energy ev')
        outfile.write(energy + ' ' + weight + '\n')
        outfile.write('end \n')

    outfile.close()        

    pass # End of writeXYZ
        


pass # End of file frozen_phonon.py
    

    
