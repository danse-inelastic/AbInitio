#!/usr/bin/env python
"""
A python module to use the ASE environment to run vasp calculations.
Adapted from bindings written by John Kitchin.

Olivier Delaire
"""

import os
import re
import weakref
import numpy as num

from ASE import Atom,ListOfAtoms

import parsing.parser2
import potcar
from parsing.SystemPM import *

class VASP:

    def __init__(self,
                 name='vasp',
                 kpts=(1,1,1),
                 pw=300,
                 xc='pawpbe',
                 vaspcmd='vasp'):
        """
        * 'name' is used in comments in the input files and the directory
        that the results are stored in.

        * 'pw' is the ENCUT parameter in the INCAR
        
        * 'kpts' is a tuple of integers that describe an a*b*c Monkhorst Pack grid

        * 'xc' is a string for the exchange-correlation functional
        """

        self.VASPPOTCARPATH = os.environ['VASPPOTCARPATH']
        self.vaspcmd = vaspcmd
        self.name = name
        self.kpts = tuple(kpts)

        #INPUT is a class object that is a dictionary to store all the
        # input variables in from the parser module
        self.incar = parser2.INPUT2('INCAR')

        self.incar['SYSTEM'] = name
        self.incar['ENCUT'] = pw

        self.incar['ISMEAR'] = 2
        self.incar['SIGMA'] = 0.2
        self.incar['PREC'] = 'Normal'
        self.incar['LREAL'] = '.False.'

        self.xc = xc
            
        #flag to be used to signify if a calculation should be run
        # True means no calculation is necessary
        # False means something has been changed that requires a recalculate
        # I assume in the beginning the calculator is not ready
        self.ready = False

    def SetName(self,name):
        self.name = name

    def _SetListOfAtoms(self,atoms):
        """Internal function used by the calculator to reference the atoms.
        to get the atoms, you must use a line like:
        self.atoms()
        note the parentheses.
        """
        self.atoms = weakref.ref(atoms)
        self.UpdateAtomsInfo()

    def UpdateAtomsInfo(self):
        """
        store info on atoms to check later if anything has changed since
        the last update
        """
        atoms = self.atoms()
        self.atoms_pos = atoms.GetCartesianPositions()
        self.atoms_uc  = atoms.GetUnitCell()
        self.atoms_atnum = atoms.GetAtomicNumbers()

    def CheckAtoms(self):
        """
        determine if the atoms have changed since they were assigned
        this calculator.  if they have, then a calculation is required
        to be current
        """
        from numpy import allclose
        
        atoms = self.GetListOfAtoms()

        #check that all atom positions are the same or not
        if not allclose(atoms.GetCartesianPositions(),self.atoms_pos,
                        rtol=1.0e-05,
                        atol=1e-08):
            print 'Cartesian positions have changed!'
            self.ready = False
        #check that the unit cell has not changed
        if not allclose(atoms.GetUnitCell(), self.atoms_uc,
                        rtol=1.0e-05,
                        atol=1e-08):
            print 'Unit cell has changed!'
            self.ready = False

        #check if an atom has been added or removed or if its identity has changed
        #this should probably raise an exception because it means the potcar must
        #be rewritten.But that is done in self.Calculate
        if len(atoms) != len(self.atoms_atnum):
            #an atom has been added or removed
            print 'The number of atoms has changed!'
            self.ready = False
        else:
            if not allclose(atoms.GetAtomicNumbers(),self.atoms_atnum):
                print 'the identity of an atom has changed!'
                self.ready = False

        return self.ready
        
    def SetPlanewaveCutoff(self,pw):
        """
        sets the ENCUT variable in the INCAR file
        """
        
        self.incar['ENCUT'] = pw

        self.ready = False

    def GetListOfAtoms(self):
        """Return the ListOfAtoms."""
        return self.atoms()

    def GetPotentialEnergy(self):
        """
        reads the energy from OUTCAR and returns it.

        first determines if a calculation is necessary.
        """
        self.Calculate()

        s=XMLSystemPM("vasprun.xml")
        energy = s.FREE_ENERGY
        return energy
    

    def GetCartesianForces(self):
        """
        return the cartesian forces on the atoms
        """
        self.Calculate()
        s=XMLSystemPM("vasprun.xml")
        f = array(s.FORCES_SEQUENCE[-1])

        return f


    def GetStress(self,units='reduced'):
        """
        Return the stress on the unit cell in eV/reduced length
        units can be 'reduced' or 'kB' kB = kBar
        """
        self.Calculate()

        if os.path.exists('OUTCAR'):
            f = open('OUTCAR','r')
            regexp = re.compile('FORCE on cell =-STRESS in cart. coord.  units \(eV/reduce length\):')
            lines = f.readlines()
            f.close()

            stressind = []
            for i,line in enumerate(lines):
                if regexp.search(line):
                    stressind.append(i)

            stressi = stressind[-1]
            
            #the stress line starts 12 or 13 lines after the
            #pattern is matched depending on the units chosen
            if units == 'reduced':
                stressi = stressi + 12
            else:
                stressi = stressi + 13
                
            fields = lines[stressi].split()
            stress = num.array([float(x) for x in fields[1:]])
            x = stress[0]
            y = stress[1]
            z = stress[2]
            xy = stress[3]
            yz = stress[4]
            zx = stress[5]

            stress_tensor =  num.array([[x, xy, zx],
                                        [xy, y , yz],
                                        [zx, yz, z]])
        else:
            stress_tensor = None

        return stress_tensor
        
    def Save(self,directoryname,force=False):
        """
        Save important files in a directory called directoryname.
        You should be able to rerun the job from all those files. 
        Do not copy over files unless force is set to True.
        """
        import shutil
        
        if not os.path.isdir(directoryname):
            os.mkdir(directoryname)

        files = ['POSCAR','POTCAR','INCAR','KPOINTS','OUTCAR','vasprun.xml']

        for f in files:
            if (os.path.exists(os.path.join(directoryname,f))
                and force is False):
                raise Exception, '%s exists, I do not want to overwrite it' % os.path.join(directoryname,f)

            shutil.copy(f,directoryname)
            
        return 1
    

    def GetMTime(self,filename):
        """
        Get the last time the file was modified.
        Usually for comparison with another file to see which one is older.
        """
        if os.path.exists(filename):
            mtime = os.stat(filename).st_mtime
        else:
            mtime = None

        return mtime

    def makePotcarFile(self): 
        """Creates the POTCAR file from the ASE list of atoms and XC type."""

        atoms = self.atoms()        
        f = open('POTCAR','w').close()
        for atom in atoms:
            symbol = atom.GetChemicalSymbol()

            if self.xc in potcar.xcs :
                print "Performing a %s-type calculation." % self.xc
                potcardir = potcar.getPath(self.xc,symbol) 
            else:
                raise ValueError, 'unrecognized XC type.'

            potcarpath = '%s/%s/POTCAR' % (self.VASPPOTCARPATH,potcardir,)
            os.system('cat %s >> POTCAR' % potcarpath)
        return

    def makePoscarFile(self):
        """ Create the POSCAR file from the ASE list of atoms.
        This should only be done when necessary. not every time.
        It is necessary only when the atoms have moved or the unit cell
        has changed. this should be explicitly checked for before
        updating the POSCAR file."""
        
        atoms = self.atoms()
        f = open('POSCAR','w')
        f.write('%s  Auto-generated by vasp.py\n' % self.name)
        ### I do not believe in scaling factors. I use the real unit cell
        f.write('1.0\n')
        uc = atoms.GetUnitCell()
        f.write('%f %f %f\n' % tuple(uc[0]))
        f.write('%f %f %f\n' % tuple(uc[1]))
        f.write('%f %f %f\n' % tuple(uc[2]))
        natoms = ''
        for atom in atoms:
            natoms += '1 '
        f.write('%s\n' % natoms) #number of each species type
        f.write('direct\n')

        sp = atoms.GetCartesianPositions()
        for pos in sp:
            f.write('%f %f %f\n' % tuple(pos))

        f.close()
        return

    def makeKpointsFile(self):
        """ Create the KPOINTS file."""
        f = open('KPOINTS','w')
        f.write('Auto-generated kpoint file\n')
        f.write('0                Automatic generation of k-points\n')
        f.write('Monkhorst        M use Monkhorst Pack\n')
        f.write('%i %i %i         grid\n' % self.kpts)
        f.write('.0 .0 .0         shift\n')
        f.close()
        return

        
    def Calculate(self):
        """
        1. figure out if a calculation needs to be run
        2. write out relevant files
        3. call vasp

        it is kind of tricky to figure out when a calculation
        should be run. some easy cases are whenever something is changed
        like atomic positions, calculation parameters, etc... or whenever
        the output file does not exist.
        """
            

        #check if the atoms have changed since the last update
        self.CheckAtoms()

        # is OUTCAR newer than POSCAR?
        if os.path.exists('OUTCAR'):
            t1 = self.GetMTime('OUTCAR')
            t2 = self.GetMTime('POSCAR')

            if t1 > t2 and self.ready:
                #outcar is newer, and a calculation is not required
                #Probably this needs to be changed in case the user wants to
                #re-run with a better accuracy
                print 'OUTCAR is newer, probably no calculation required'
                return
            

        #if you get here, a calculation should be run
        print 'running a calculation'

        self.makePotcarFile()
        self.makePoscarFile()

        ### create INCAR. This should be redundant, but I do it anyway.
        self.incar.Write()

        ### create the KPOINTS file:
        self.makeKpointsFile()
            
        ### Now run vasp
        status = os.system(self.vaspcmd)

        self.UpdateAtomsInfo()
        self.ready = True
        return status


    def ReadAtoms(name='.'):
        """Static method to read in the atoms."""

        f = open('POSCAR','r')
        lines = f.readlines()
        f.close()

        comment = lines[0]

        uc = []
        uc.append([float(x) for x in lines[2].split()])
        uc.append([float(x) for x in lines[3].split()])
        uc.append([float(x) for x in lines[4].split()])

        scalefactor = float(lines[1].strip())
        if scalefactor < 0:
            #that means this is the volume of the cell
            #vol = determinant of the uc
            vol0 = abs(num.determinant(uc))
            uc = abs(scalefactor)/vol0 * uc
        else:
            uc = scalefactor * num.array(uc)

        atomcounts = [int(x) for x in lines[5].split()]
        natoms = 0
        for count in atomcounts:
            natoms += count
            
        if lines[6][0] in ['s','S']:
            #selective dynamics were chosen, and positions start on line 7
            coordsys = lines[7][0]
            poscounter = 8
        else:
            coordsys = lines[6][0]
            poscounter = 7

        positions = []
        for i in range(natoms):
            pos = num.array([float(x) for x in lines[poscounter + i].split()])
            if coordsys[0] in ['C','c','K','k']:
                #cartesian coordinates, do nothing
                pass
            else:
                #direct coordinates. calculate cartesian coords
                pos = pos[0]*uc[0] + pos[1]*uc[1] + pos[2]*uc[2]

            positions.append(pos)

        positions = num.array(positions)

        #now get the identities from the POTCAR file.
        f = open('POTCAR','r')
        lines = f.readlines()
        f.close()
        #the start of each psp is either 'US symbol' 'PAW_GGA symbol
        #comment' or 'PAW_PBE symbol comment'
        tag,symbol = lines[0].split()
        regexp = re.compile('^\s+%s' % tag)

        psps = []
        for line in lines:
            if regexp.search(line):
                psps.append(line)

        #print atomcounts, psps
        if len(atomcounts) != len(psps):
            raise Exception,'number of atom counts in POSCAR does not equal # of psps in POTCAR'

        from ASE import Atom,ListOfAtoms
        atoms = ListOfAtoms([])
        poscounter = 0
        for i,count in enumerate(atomcounts):
            tag,symbol = psps[i].split()
            for j in range(count):
                atoms.append(Atom(symbol,position=positions[poscounter]))
                poscounter += 1

        atoms.SetUnitCell(uc,fix=True)

        if os.path.exists('OUTCAR'):
            f = open('OUTCAR','r')
            
            regexp = re.compile('TOTAL-FORCE \(eV/Angst\)')
            lines = f.readlines()
            f.close()

            forcei = None
            for i,line in enumerate(lines):
                if regexp.search(line):
                    forcei = i+2 #linenumber that forces start on
                    break
            if forcei is None:
                raise Exception,'forcei is none, no forces found!'
            
            forces = []
            for i in range(len(atoms)):
                posforce = [float(x) for x in lines[forcei + i].split()]
                forces.append(posforce[3:])
        else:
            forces = [None for atom in atoms]

        for atom, force in zip(atoms, forces):
            #print force
            atom._SetCartesianForce(force)
            
        ### more code must be added to read in the calculator.
        calc = VASP()
        calc.incar = parser2.INPUT2('INCAR')
        return atoms

    ReadAtoms = staticmethod(ReadAtoms)   
            

def test1(pw=300, kpts=(6,6,6), vaspcmd='mr vasp'):
    """Performs a test single-point energy calculation with vasp module for B2 FeAl."""

    print "Test single-point energy calculation with vasp module for B2 FeAl."
    
    v = VASP(pw=pw,
             kpts=kpts,
             name='FeAl',
             vaspcmd=vaspcmd)
    
    from ASE import ListOfAtoms, Atom

    atoms = ListOfAtoms([Atom('Fe', [0.0,0.0,0.0]),
                         Atom('Al', [0.5,0.5,0.5])])

    uc = [[2.87, 0.00, 0.00],
          [0.00, 2.87, 0.00],
          [0.00, 0.00, 2.87]]
    
    atoms.SetUnitCell(uc)
    atoms.SetCalculator(v)

    print "Potential energy: ", v.GetPotentialEnergy()
    return

def test2(pw=300, kpts=(6,6,6), vaspcmd='mr vasp'):
    """Performs a test stress-strain calculation with vasp module for fcc Pt."""

    v = VASP(pw=pw,
             kpts=kpts,
             name='fccPt',
             vaspcmd=vaspcmd)

    from ASE import ListOfAtoms, Atom

    atoms = ListOfAtoms([Atom('Pt',[0,0,0]),
                         Atom('Pt',[0.5,0.5,0]),
                         Atom('Pt',[0.5,0,0.5]),
                         Atom('Pt',[0,0.5,0.5])])

    uc = [[4.05, 0, 0],
          [0.0, 4.05, 0],
          [0.0, 0, 4.05]]

    atoms.SetUnitCell(uc)

    atoms.SetCalculator(v)


    import ASE.Utilities.GeometricTransforms as ASEgeotrans
    
    initvol = atoms.GetUnitCellVolume()

    vols, energies = [],[]

    #for f in [0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15]:
    for f in [0.9, 0.95, 1.0, 1.05, 1.1]:
        ASEgeotrans.SetUnitCellVolume(atoms,f*initvol)

        #v.SetName('%1.2f_eos' % f)
        #print v.GetStress()
        vols.append(atoms.GetUnitCellVolume())
        energies.append(atoms.GetPotentialEnergy())


    import pylab as pl

    pl.plot(vols,energies,'ko ')
    pl.show()

    import ASE.Utilities.EquationOfState as ASEeos

    eos = ASEeos.EquationOfState('Murnaghan',vols,energies)
    # print the minimum volume, energy, bulk modulus and pressure
    print eos         
    g = eos.GetPlot()
    eos.SavePlot('murn.png') #save the figure as a png
    return

def test3(pw=129, kpts=(4,4,4), vaspcmd='mr vasp', xc='uspplda'):

    v = VASP(pw=pw,
             kpts=kpts,
             name='fccAl',
             vaspcmd=vaspcmd)

    from ASE import ListOfAtoms, Atom
    
    atoms = ListOfAtoms([Atom('Al',[0,0,0]),
                         Atom('Al',[0.5,0.5,0]),
                         Atom('Al',[0.5,0,0.5]),
                         Atom('Al',[0,0.5,0.5])])

    uc = [[4.05, 0, 0],
          [0.0, 4.05, 0],
          [0.0, 0, 4.05]]

    atoms.SetUnitCell(uc)
    atoms.SetCalculator(v)


    import ASE.Utilities.GeometricTransforms as ASEgeotrans
    
    initvol = atoms.GetUnitCellVolume()

    vols, energies = [],[]

    for f in [0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15]:
        ASEgeotrans.SetUnitCellVolume(atoms,f*initvol)
        v.SetName('%1.2f_eos' % f)
        print v.GetStress()
        vols.append(atoms.GetUnitCellVolume())
        energies.append(atoms.GetPotentialEnergy())

    import pylab as pl

    pl.plot(vols,energies,'ko ')
    pl.show()

    import ASE.Utilities.EquationOfState as ASEeos

    eos = ASEeos.EquationOfState('Murnaghan',vols,energies)
    # print the minimum volume, energy, bulk modulus and pressure
    print eos         
    g = eos.GetPlot()
    eos.SavePlot('murn.png') #save the figure as a png
    return

if __name__ == '__main__':

    test1(pw=264, kpts=(6,6,6))
#    test2(pw=300, kpts=(6,6,6))
