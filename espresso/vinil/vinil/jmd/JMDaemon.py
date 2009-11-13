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
from pyre.applications.Daemon import Daemon


class JMDaemon(Application, Daemon):


    class Inventory(Application.Inventory):
        import pyre.inventory

        # ?
        depositories = pyre.inventory.list("depositories")
        depositories.meta['tip'] = """extra depositories for my harnessed components"""

    def __init__(self, name):
        Application.__init__(self, name, facility='daemon')
        Daemon.__init__(self)
        return

    def main(self, *args, **kwds):
        print "Hello world!"

    def _configure(self):
        super(JMDaemon, self)._configure()

        self.depositories = self.inventory.depositories
        return


    def _init(self):
        super(JMDaemon, self)._init()

        curator = self.getCurator()
        curator.addDepositories(*self.depositories)

        return
    
__date__ = "$Nov 12, 2009 12:35:26 PM$"


