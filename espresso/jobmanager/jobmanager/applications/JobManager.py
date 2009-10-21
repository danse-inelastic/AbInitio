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

from pyre.applications.Application import Application

class JobManager(Application):
    
    class Inventory(Application.Inventory):
        import pyre.inventory

    def __init__(self, name):
        Application.__init__(self, name, facility='something')

    def _configure(self): pass

    def _init(self): pass

__date__ = "$Oct 21, 2009 7:11:53 AM$"


