#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# Performs self-consistent calculation and parses PW output file

# Requires diffpy.Structure
from qecalc.qecalc import QECalc
qe  = QECalc("multi-config.ini")
qe.pwscfLauncher()
energy  = qe.getTotalEnergy()[0]
force   = qe.getForces()[0]
print "Example: (using QECalc)"
print "Total energy = %f" % energy
print "Force        = %s" % force



__date__ = "$Oct 15, 2009 3:20:12 PM$"


