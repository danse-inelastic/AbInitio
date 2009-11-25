#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Nikolay Markovskiy
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import string
import qe_io_dict as io_dict
import numpy
from baseoutput import BaseOutput


class Output(BaseOutput):

    def __init__(self):
        BaseOutput.__init__(self)
        self.parsers = {
                    'single phonon' : self.getSinglePhonon
                }

    def getSinglePhonon(self, setting):
        """Obtain a list of phonon modes from output generated by dynmat.x"""
        return [self.dynmat_modes( setting.dynmatOutput )[1]]
        dynmatOut = read_file(setting.dynmatOutput)
        keyStart = find_key_from_marker_string(dynmatOut, '#', 'mode')
        modeNum = 1
        key = keyStart + 1
        words = string.split(dynmatOut[key])
        modes = []
        while words[0]	 == str(modeNum):
            modes.append( float(words[1]))
            key = key + 1
            modeNum = modeNum + 1
            words = string.split(dynmatOut[key])
        return [(modes, 'cm-1')]
    
    def dynmat_modes(self, fname):
        matdynDict = io_dict.read_file(fname)
        qKeys = io_dict.find_all_keys_from_string(matdynDict, 'q =')
        qP = []
        for i in qKeys:
            qP.append( [ float(qi) for qi  in string.split( matdynDict[ i ] )[2:5] ] )
        qPoints = numpy.array(qP)


    # find number of atoms  per unit cell and dimensionality
    # get frequencies keys for the last q-point:
        fKeys = io_dict.find_all_keys_from_string_afterkey( matdynDict, qKeys[-1], 'omega')
    #  get sample displacement vector for the first frequency of the last q point and find its dimensionality
    #  here omegaShift = 1 means no space between displacements and frequencies (2 means 1 extra line ...)
        omegaShift = 1
        nDim =  len( matdynDict[ fKeys[0] + omegaShift  ].split()[1:-1] )/2
    #  get number of atoms in unit cell
        nAtom =  fKeys[1] - fKeys[0] - omegaShift

    # qShift = 2 means 1 exra line between q line and frequencies line
        qShift = 2

    #  create numpy array in accordance with idf format specification

    #    Pol = [ [ [ [ ] ] ] ]
        Pol = []
        Omegas = []
        for i in qKeys:
    #        Mode = [ [ [ ] ] ]
            Mode = []
            Omega = []
            for iOmega in range( i + qShift, nDim*nAtom*(nAtom + omegaShift) +  i + qShift, nAtom+omegaShift):
            # get omegas in THz:

    #            print float( matdynDict[ iOmega].split('=')[1].split()[0] )
                Omega.append( float( matdynDict[ iOmega].split('=')[2].split()[0] ) )
    #	    Atom = [ [ ] ]
                Atom = []
                for iAtom in range( iOmega + omegaShift, iOmega + omegaShift + nAtom ):
                    vecStr = matdynDict[ iAtom ].split()[1:-1]
    #		        vec = [  ]
                    vec = [ float(v1) + 1.0j*float(v2) for v1,v2 in zip( vecStr[:-1:2], vecStr[1::2] ) ]
                    Atom.append(vec)
                Mode.append( Atom )
            Pol.append( Mode )
            Omegas.append( Omega )
        npOmega = numpy.array(Omegas)
        npPol = numpy.array(Pol)
    #    THz2meV = 4.1357 # meV * s
    #    output Omega in cm-1
        return (npPol, None), (npOmega, 'cm-1'), (qPoints, None)

if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 18, 2009 7:28:06 PM$"