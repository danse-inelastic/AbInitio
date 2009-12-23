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

import mpi4py
import os

inputs = ['pwscfInput', 'd3Input', 'dosInput', 'dynmatInput', 'matdynInput', 'phInput', 'q2rInput']
outputs = [pwscfOutput, d3Output, d3fildyn, dosOutput, fldos, dynmatOutput, matdynOutput, matdynModes, matdynFreqs, matdynfldos, phOutput, phFildyn, q2rOutput ]



if __name__ == "__main__":
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Dec 14, 2009 2:47:12 PM$"
