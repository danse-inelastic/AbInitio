import os, sys
import vasp.parsing.parser2
from vasp.vasp import VASP
from vasp.parsing.Structure import Structure
from vasp.parsing import parser2
from pyphon._pyphon import phon

__doc__ = """A module to implement a Python Phonon calculator,
based on the PHON Fortran code."""

vaspcmd = 'vasp'

class OldPhonCalc:
    """A Python Phonon calculator,
    based on the PHON Fortran code."""
    def __init__(self,
                 name='phoncalc',
                 cellsize=[1,1,1],
                 amplitude=0.01,
                 kgrid=[20,20,20],
                 dosmin=0.0,
                 dosmax=20.0,
                 dosstep=0.1,
                 dossmear=0.1,
                 temperature=300
                 ):
        """
        * 'name' is the nameof the calculation
        * 'amplitude' is the amplitude of the atomic displacement(s),
        in fraction of the unitcell length along displacement direction.
        * 'cellsize' is the size of the supercell, in numbers (integer) of
        replications of base unit cell.
        * 'kgrid' is the order of the Monkhorst-Pack grid for phonon k-points
        * 'dosmin' is the energy minimum for the phonon DOS histogram
        * 'dosmax' is the energy maximum for the phonon DOS histogram
        * 'dosstep' is the step size for the phonon DOS histogram
        * 'dossmear' is the energy broadening of the phonon DOS
        * 'temperature' is the temperature for calculating phonon free energy
        """
        
        self._name = name
        self._cellsize = list(cellsize)
        self._amplitude = amplitude
        self._kgrid = kgrid
        self._dosmin = dosmin
        self._dosmax = dosmax
        self._dosstep = dosstep
        self._dossmear = dossmear
        self._temperature = temperature
        self._makesupercell = True

        # number of types of ions and masses need to be passed in...:
        self._ntypes = 2
        self._masses = [26.98, 55.85]

        # setup the inphon dictionary
        self.inphon = parser2.INPUT2('INPHON')
        self.inphon['NTYPES'] = self._ntypes
        self.inphon['MASS'] = self._masses
        self.inphon['NDIM'] = cellsize
        self.inphon['DISP'] = int(1.0 / self._amplitude)
        self.inphon['QA'] = kgrid[0]
        self.inphon['QB'] = kgrid[1]
        self.inphon['QC'] = kgrid[2]
        self.inphon['DOSIN'] = dosmin
        self.inphon['DOSEND'] = dosmax
        self.inphon['DOSSTEP'] = dosstep
        self.inphon['DOSSMEAR'] = dossmear
        self.inphon['TEMPERATURE'] = temperature
        
        # default values for the inphon dictionary
        self.inphon['RMAX'] = 15
        self.inphon['LFREE'] = '.TRUE.'
        self.inphon['LSYMM'] = '.TRUE.'
        self.inphon['NTI'] = 50
        self.inphon['IPRINT'] = 3       # required to get eigenvectors
        self.inphon['LGAMMA'] = '.FALSE.'
        self.inphon['LSUPER'] = '.TRUE.'

        # write the INPHON file
        self.inphon.Write()
        return

    def generateSupercell(self,cellsize=None):
        """Generate supercell and displacements by launching the Phon exe."""

        if self._makesupercell == False:
            return
        else:
            if cellsize == None:
                cellsize = self._cellsize
            self._cellsize = cellsize
            self.inphon['NDIM'] = self._cellsize
            self.inphon['LSUPER'] = '.TRUE.'
            # write the INPHON file:
            self.inphon.Write()
            # run Phon to build the supercell and calculate displacements:
            # This is assuming that the POSCAR exists...
            phon()

            # now copy the SPOSCAR over to POSCAR
            # (this will only work on linux...)
            #os.system('cat %s >> POTCAR' % potcarpath)
            os.system('cp POSCAR POSCAR_eq')
            os.system('cp SPOSCAR SPOSCAR_eq')
            os.system('cp SPOSCAR POSCAR')
            
            # prevent the generation of an even larger supercell:
            self._makesupercell = False
            self.inphon['LSUPER'] = '.FALSE.'
            self.inphon.Write()
            return

    def _makePosFiles(self):
        """Generates the Poscar files for all atomic displacements."""
        from vasp.parsing.matrix import Vector
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
    

    def calculateForces(self):
        """Calls the computation engine to calculate the forces,
        for every distorted supercell."""
        # This will call VASP for now
        self._makePosFiles()
        # For now, just run the 'runhf' script that calls VASP on all structures.
        os.popen('./runhf')
        os.system('cp SPOSCAR_eq SPOSCAR')
        print "Finished computing the forces for all displacements."
        return

    def _gatherForces(self):
        """This is a helper function to collect all the forces on the atoms,
        for every distorted supercell."""
        from vasp.parsing.SystemPM import XMLSystemPM
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
        self.inphon['LSUPER'] = '.FALSE.'
        try:
            shutil.copyfile('SPOSCAR_eq', 'SPOSCAR')
            #os.system('cp SPOSCAR_eq POSCAR')
        except:
            raise IOError, 'SPOSCAR_eq not found.'
        print "Calling the Phon executable."
        open('phon.out',"w").write(''.join(os.popen(phonexe).readlines()))
        #phon()
        #os.system('phon > phon.out')
        print "Done running Phon."
        return

# we also need to implement parsers to retrieve the DOS, dispersions, polarizations, etc...
