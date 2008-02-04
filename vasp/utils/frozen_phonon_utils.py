import os
import sys
import numpy as np
import scipy
import scipy.io

try:
    from vasp.parsing.SystemPM import XMLSystemPM
    from vasp.parsing.Structure import Structure
    #import vasp.parsing.matrix as p4mat
except ImportError:
    print "P4Vasp could not be imported in Python."

def writeXYZ(name, strains, outfilename):
    """Reads the output XML files from VASP computation of frozen_phonons
    and creates a XYZ format file for input into GULP."""
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
        
def parseConfigsAndEnergies(name, strains):
    """Parses the strained structures, from a list of VASP output XMl files,
    as a LsitOfAtoms list,
    and a list of corresponding energies."""
    from crystalIO.converters import p4vaspStruct2UnitCell,unitCell2ListOfAtom
    atomicConfigs = []
    energies = []
    for strain in strains:
        try:
            xmlfilestring = 'vasprun_'+name+'_u_'+repr(strain)+'.xml'
            xmlrun = XMLSystemPM(xmlfilestring)
        except:
            raise IOError, 'VASP XML output file '+xmlfilestring+' could not be opened.'        
        energy = xmlrun.FREE_ENERGY
        struct = xmlrun.FINAL_STRUCTURE
        # set the atomic positions to fractional coordinates:
        struct.setDirect()
        uc = p4vaspStruct2UnitCell(struct)
        loa = unitCell2ListOfAtom(uc)
        print "x = %s \n" %  strain
        print "Atomic configuration: \n"
        print loa
        print "Energy = %s" % energy
        atomicConfigs.append(loa)
        energies.append(energy)
    return atomicConfigs, energies

def writeConfigsAndEnergies(name, strains, pklfilename):
    """Parses the strained structures, from a list of VASP output XMl files,
    as a ListOfAtoms list,
    and a list of corresponding energies,
    and writes everything to a pickle file."""
    import pickle
    atomicConfigs, energies = parseConfigsAndEnergies(name, strains)
    pklfile = open(pklfilename, 'w')
    pickle.dump((atomicConfigs, energies), pklfile)
    pklfile.close()
    pass # End of writeUCsAndEnergies

