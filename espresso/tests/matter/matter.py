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

import matter
from qecalc.pwcalc import PWCalc

configString = """
[pw.x]
pwInput: scf.in
"""

def testMatter():

    pwcalc = PWCalc( configString = configString )
    pwcalc.pw.syncSetting()

    

if __name__ == "__main__":
    testMatter()

__author__="Nikolay Markovskiy"
__date__ ="$Jan 13, 2010 11:46:07 AM$"
