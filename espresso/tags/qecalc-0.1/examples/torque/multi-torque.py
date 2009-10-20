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

# Example: Finding total energy and force for Ni
# Performs self-consistent calculation using torque and parses PW output file
# Make sure to create 'temp' directory
#   Tested on octopod

# Requires diffpy.Structure

from qecalc.qecalc import QECalc
qe  = QECalc("multi-config.ini")
qe.pwscfLauncher()
energy  = qe.getTotalEnergy()[0]
force   = qe.getForces()[0]
print "Example: (using QECalc)"
print "Total energy = %f" % energy
print "Force        = %s" % force



__date__ = "$Oct 14, 2009 5:19:05 PM$"

