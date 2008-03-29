# Olivier Delaire

__doc__ = """A modules to implement a First-Principles Phonon Calculator with Phon."""

import os, sys
import shutil
from AbInitio.vasp.parsing import parser2
from AbInitio.vasp.parsing.Structure import Structure
from AbInitio.vasp.parsing.matrix import Vector
from AbInitio.AbiPhon.AbiPhonCalc import AbiPhonCalc
from pyphon._pyphon import phon

class PhonCalc(AbiPhonCalc):
    """A first-principles phonon calcor based on Phon."""

    def __init__(self, unitcell,
                 name='phon',
                 supersize=[1,1,1],
                 qgridsize=[10,10,10],
                 abiCalc=None,
                 dosmin=0.0, dosmax=50.0,
                 dosstep=0.2, dossmear=0.2,
                 temperature=300,
                 **miscargs):
        AbiPhonCalc.__init__(self,unitcell, supersize, abiCalc=abiCalc, qpts=None)
        self._name = name
        self._qgridsize = qgridsize
        self._dosmin = dosmin
        self._dosmax = dosmax
        self._dosstep = dosstep
        self._dossmear = dossmear
        self._temperature = temperature
        
        self._atomTypesNums = unitcell.getAtomTypeDenum()
        self._numtypes = len(self._atomTypesNums)
        self._masses = [a[1] for a in self._atomTypesNums]

        # uses a input parser that is derived from ordered dictionary,
        # and writes to file every time an entry is modified:
        self._inphon = parser2.INPUT2('INPHON')
        self._setupPhon()

        for keyarg in miscargs:
            self._inphon[keyarg] = miscargs[keyarg] 

        pass # end of __init__

    
    def _setupPhon(self):
        """helper function to setup Phon calculation."""

        # setup the inphon dictionary
        self._inphon = parser2.INPUT2('INPHON')
        self._inphon['NTYPES'] = self._numtypes
        self._inphon['MASS'] = self._masses
        self._inphon['NDIM'] = self._supersize
        self._inphon['DISP'] = int(1.0 / self._amplitude)
        self._inphon['QA'] = self._qgridsize[0]
        self._inphon['QB'] = self._qgridsize[1]
        self._inphon['QC'] = self._qgridsize[2]
        self._inphon['DOSIN'] = self._dosmin
        self._inphon['DOSEND'] = self._dosmax
        self._inphon['DOSSTEP'] = self._dosstep
        self._inphon['DOSSMEAR'] = self._dossmear
        self._inphon['TEMPERATURE'] = self._temperature
        
        # default values for the inphon dictionary
        self._inphon['RMAX'] = 15
        self._inphon['LFREE'] = '.TRUE.'
        self._inphon['LSYMM'] = '.TRUE.'
        self._inphon['NTI'] = 50
        self._inphon['IPRINT'] = 3       # required to get eigenvectors
        self._inphon['LGAMMA'] = '.FALSE.'
        self._inphon['LSUPER'] = '.TRUE.'

        return

    def _writeInphon(self):
        """Writes the INPHON input file for Phon."""
        self._inphon.write()
        return

    def name(self):
        """returns the name of the calculation."""
        return self._name

    def printInputs(self):
        """Prints out all the input arguments that are used to run Phon."""
        for key in self._inphon:
            print key, self._inphon[key]

    def setQgridSize(qgridsize=[10,10,10]):
        """Set the size of the (regular) grid of Q-points at which the phonons
        will be calculated."""
        self._qgridsize = qgridsize
        self._inphon['QA'] = self._qgridsize[0]
        self._inphon['QB'] = self._qgridsize[1]
        self._inphon['QC'] = self._qgridsize[2]


    def genSupercell(self, superdims=None):
        """Generates a supercell of dimensions given by superdims."""

        # This part should be replaced by a call to AbiPhonCalc method
        if superdims is None:
            superdims = self._supersize
        if superdims is None:
            raise ValueError, 'supercell should be integer multiple of unit cell.'
        from crystal.UnitCell import UnitCell
        from crystal.Atom import Atom
        # generate a supercell with multiplied lattice vectors:
        supercell = UnitCell(self._unitcell)
        cellvectors = self._unitcell.getCellVectors()
        supercellvectors = cellvectors * superdims
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
                        newsite = Site(newpos, site.getAtom())
                        supercell.addSite(newsite, '')
        self._supercell = supercell
        
        # prevent the generation of an even larger supercell:
        self._superCellReady = True

        # Phon-specific actions:
        self._inphon['LSUPER'] = '.FALSE.'
        self._inphon.Write()
        # write the SPOSCAR file
        from crystal.crystalIO.converters import unitCell2P4vaspStruct
        superstruct = unitCell2P4vaspStruct(self._supercell)
        superstruct.write(f='SPOSCAR', newformat=0)

        return
   
    def genPhonSupercell(self, supersize=None):
        """Generate supercell and displacements by launching Phon."""
        if self._superCellReady:
            return
        else:
            from crystal.crystalIO.converters import unitCell2P4vaspStruct, p4vaspStruct2UnitCell
            from AbInitio.vasp.parsing.Structure import Structure
            
            if supersize == None:
                supersize = self._supersize
            self._supersize = supersize
            self._inphon['NDIM'] = self._supersize
            self._inphon['LSUPER'] = '.TRUE.'
            # write the INPHON file:
            self._inphon.Write()

            # First, we need to write the unitcell to a POSCAR for Phon():
            
            struct = unitCell2P4vaspStruct(self._unitcell)
            try:
                struct.write(f='POSCAR', newformat=0)
            except:
                raise IOError, 'Could not write the structure POSCAR to file.'
            
            # run Phon to build the supercell and calculate displacements:
            #phon()
            

            # now copy the SPOSCAR over to POSCAR
            try:
                shutil.copy('POSCAR', 'POSCAR_eq')
            except:
                pass 
            try:
                shutil.copy('SPOSCAR', 'SPOSCAR_eq')
                shutil.copy('SPOSCAR', 'POSCAR')
            except:
                raise IOError, 'Error while copying SPOSCAR.'

            # load the supercell created by Phon (copied to POSCAR)
            struct = Structure('POSCAR')
            uctmp = p4vaspStruct2UnitCell(struct)
            self._supercell = uctmp   # this may be a problem if VASP atom types are not passed correctly...
            
            # prevent the generation of an even larger supercell:
            self._superCellReady = True
            self._inphon['LSUPER'] = '.FALSE.'
            self._inphon.Write()
            return
        

    def _makePosFiles(self):
        """Generates the Poscar files for all atomic displacements."""
        try:
            struct = Structure("POSCAR")
        except:
            raise IOError, "POSCAR file not found: make sure supercell (SPOSCAR) was generated."

        try:
            dispfile=open("DISP",'r')
        except:
            raise IOError, "DISP file not found, make sure supercell (SPOSCAR) was generated."
        # this section creates a 'runhf' file,
        # that can be used to lauch the batch vasp job
        runfile=open("runhf", 'w')
        runfile.write("#!/bin/bash \n")
        runfile.write("VASP='vasp' \n\n")
        # end section

        disp=[]
        pos=[]
        pos.append(struct)
        pos[len(pos)-1].write('pos'+repr(len(pos)-1))

        for line in dispfile:       
            print line
            str=eval(line.strip(' \\ \n'))
            liststr = str.split()
            print liststr
            listnum = []
            for item in liststr:
                listnum.append(eval(item))
            disp.append(listnum)
            atomindex = listnum[0]
            atomdisp = listnum[-3:]
            pos.append(struct)
            print "Length of positions file list:", len(pos) 
            pos[len(pos)-1].positions[atomindex-1]=pos[len(pos)-1].positions[atomindex-1]+Vector(atomdisp)
            print "Writing pos file ", len(pos)-1
            pos[len(pos)-1].write('pos'+repr(len(pos)-1), newformat=0)
            pos[len(pos)-1].positions[atomindex-1]=pos[len(pos)-1].positions[atomindex-1]-Vector(atomdisp)
            # this section creates a 'runhf' file,
            # that can be used to lauch the batch vasp job
            runfile.write("cp pos"+repr(len(pos)-1)+" POSCAR \n")
            runfile.write("$VASP >>out.vasp 2>>err.vasp & \n")
            runfile.write("wait \n")
            runfile.write("cp OUTCAR out_"+repr(len(pos)-1)+" \n\n")
            runfile.write("cp CHGCAR CHGCAR_"+repr(len(pos)-1)+" \n\n")
            runfile.write("cp vasprun.xml vasprun.xml_"+repr(len(pos)-1)+" \n\n")
            # end section

        dispfile.close()
        runfile.close()
        os.system('chmod +x runhf')
        return # end of _makePosFiles()

    def calcForces(self):
        """This calls the abinitio engine to evaluate the forces on the atoms,
        for all the distorted supercells."""
        # For now, just run the 'runhf' script that calls VASP on all structures.
        os.popen('./runhf')
        os.system('cp SPOSCAR_eq SPOSCAR')
        print "Finished computing the forces for all displacements."
        return

    def _gatherForces(self):
        """This is a helper function to collect all the forces on the atoms,
        for every distorted supercell."""
        from AbInitio.vasp.parsing.SystemPM import XMLSystemPM
        # we get back the undistorted base cell into POSCAR
        os.system('cp SPOSCAR_eq POSCAR')
        struct = Structure("POSCAR")
        dispfile=open("DISP",'r')
        displist=[]
        forcefile=open("FORCES", 'w')

        for line in dispfile:
            str=eval(line.strip(' \\ \n'))
            displist.append(str)

        dispfile.close()
        print "Number of displacements:", len(displist)
        print "Displacements:"
        
        for disp in displist:
            print disp.strip()

        npos=len(displist)
        forcefile.write(repr(npos)+"\n")
        posnum=range(npos)

        for i in posnum:
            run=XMLSystemPM('vasprun.xml_'+repr(i+1))
            forces=run.FORCES_SEQUENCE
            forcefile.write(displist[i])
            forcefile.write("\n")
            for atomforce in forces[0]:
                fx = "%8f" % atomforce[0]
                fy = "%8f" % atomforce[1]
                fz = "%8f" % atomforce[2]
                forcefile.write(fx+' '+fy+' '+fz+'\n')     

        forcefile.close()
        return
        
    def calcPhonons(self, phonexe='phon'):
        """Generate the Phonons (DOS and/or dispersions)."""
        import shutil
        # gather forces
        # this currently assumes that VASP has written output files...
        self._gatherForces()
        
        # Call the Phon exec again to run BvK part of calculation.
        self._makesupercell = False        
        self._inphon['LSUPER'] = '.FALSE.'
        try:
            shutil.copyfile('SPOSCAR_eq', 'SPOSCAR')
            #os.system('cp SPOSCAR_eq POSCAR')
        except:
            raise IOError, 'SPOSCAR_eq not found.'
        print "Calling the Phon executable."
        outfile = open('phon.out',"w")
        outfile.write(''.join(os.popen(phonexe).readlines()))
        #phon()
        #os.system('phon > phon.out')
        print "Done running Phon."
        return
