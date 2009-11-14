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

class EchoDaemon(Application, Daemon):

    class Inventory(Application.Inventory):
        import pyre.inventory
        home = pyre.inventory.str("home", default="/tmp")   # ?


    def main(self, *args, **kwds):
        pass

    def __init__(self, name):
        Application.__init__(self, name, facility='daemon')
        Stager.__init__(self)


    def _configure(self):
        super(Daemon, self)._configure()

        import os
        self.home = os.path.abspath(self.inventory.home)


if __name__ == "__main__":
    echo    = EchoDaemon()
    echo.run()


#    def _init(self):
#        super(Daemon, self)._init()
#
#        curator = self.getCurator()
#        curator.addDepositories(*self.depositories)
#
#        return
#        clientName = pyre.inventory.str("client-name")
#        serverInfo = pyre.inventory.str("server-info", default="server-info")



__date__ = "$Nov 13, 2009 4:28:05 PM$"


