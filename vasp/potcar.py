#!/usr/bin/env python


'''
This file contains the directory names where the POTCAR files for each
element are stored.

You need to define an environement variable called VASPPATH that contains
the path to the Potcar_Files: e.g. in csh:

setenv VASPPATH $HOME/myPython/VASP/Vasp/Potcar_Files

In a script, if you want to change the potcar file you would


from VASP.potcar import *

gga['C'] = 'selected_C_potcar_dir'


Eventually this file should contain dictionaries for each XC available, e.g.:
USPP LDA
USPP GGA
PAW  GGA

etc...
and the right dictionary will get loaded depending on the XC chosen in the
input files.


modified: Olivier Delaire 05/17/07
tried to implement access to different XC POTCARs and uniform method
to retrieve POTCAR from XC-type and atomic symbol.
Still need to distinguish between different POTCARs according to number of valence electrons (eg: Fe / Fe_pv, Al / Al_h)
'''


xcs = {'uspplda' : 'pot',
      'usppgga' : 'pot_GGA',
      'pawlda' : 'potpaw',
      'pawgga' : 'potpaw_GGA',
      'pawpbe' : 'potpaw_PBE',
      }

uspplda = {'Al' : 'Al'}

usppgga = {'H':'H_200eV',
           'Pt':'Pt',
           'C':'C',
           'O':'O',
           'Sn':'Sn',
           'Al':'Al',
           'Fe':'Fe'}

pawlda = {}

pawgga = {'Al' : 'Al',
          'V' : 'V',
          'Fe' : 'Fe',
          'Pt' : 'Pt',
          'Si' : 'Si'}

pawpbe = {'Al' : 'Al',
          'V' : 'V',
          'Fe' : 'Fe',
          'Pt' : 'Pt',
          'Si' : 'Si'}

def getPath(xc, symbol):
    """Returns a path string for the location of the POTCAR file."""

    if xc not in xcs:
        raise ValueError, "Unrecognized XC type."
    else:
        mapper = eval(xc)
        filename = mapper.get(symbol)
        if not filename: filename = symbol # default
        dir = xcs[xc]
        return os.path.join( dir, filename)

import os
