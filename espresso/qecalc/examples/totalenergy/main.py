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

from property import Property

p  = Property("config.ini")
energy  = p.getTotalEnergy()[0]
force   = p.getForces()[0]
print "Total energy = %f" % energy
print "Force        = %s" % force



__date__ = "$Oct 14, 2009 5:19:05 PM$"


