import os
import numpy as np
from p4vasp.SystemPM import *
import p4vasp.matrix as p4mat


def CalcStrains(strains, name, vaspcmd):
    """Calls vasp to calculate a sequence of structures,
    corresponding to POSCAR files derived from name,
    appended with '_x', for x in strains."""

    # Call the vasp executable (vaspcmd) through a os.system call
    # for every structure "name" found, appended with x in strains

    for x in strains:
        print "strain=", x
        os.system('cp pos_'+name+'_'+repr(x)+' POSCAR')
        os.system(vaspcmd)
        os.system('cp vasprun.xml vasprun_'+name+'_'+repr(x)+'.xml')
        os.system('cp OUTCAR OUTCAR_'+name+'_'+repr(x))


    
