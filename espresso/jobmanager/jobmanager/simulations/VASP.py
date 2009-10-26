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

"""
Runs simulation for VASP
"""

from Simulator import Simulator

class VASP(Simulator):

    def render(self, computation, db=None, dds=None):
        self.dds = dds

        self._files = []

        return self._files

__date__ = "$Oct 25, 2009 7:11:37 PM$"


