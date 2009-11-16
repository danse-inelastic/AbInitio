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
from EchoService import EchoService

class EchoDaemon(Application, Daemon):

    class Inventory(Application.Inventory):
        import pyre.inventory
        home = pyre.inventory.str("home", default="/tmp")   # ?


    def main(self, *args, **kwds):
        
        service  = self.configureComponent(EchoService, 'service') #    EchoService()
        service.serve()
        
        # harness service?

#        service.init()
#        service.eventLoop()


    def __init__(self, name=None):
        Application.__init__(self, name, facility='daemon')
        Daemon.__init__(self)

    def _defaults(self):
        print "defaults"
        Application._defaults(self)


#    def _init(self):
#        curator = self.getCurator()
#        curator.addDepositories(*self.home)


    def _configure(self):
        print "configure"
        super(EchoDaemon, self)._configure()

        import os
        self.home = os.path.abspath(self.inventory.home)


if __name__ == "__main__":
    echo    = EchoDaemon('echo')
    echo.run(spawn=False)


#    def _init(self):
#        super(Daemon, self)._init()
#
#        curator = self.getCurator()
#        curator.addDepositories(*self.depositories)
#
#        return
#        clientName = pyre.inventory.str("client-name")
#        serverInfo = pyre.inventory.str("server-info", default="server-info")

#    def _getPrivateDepositoryLocations(self):
#        return ['.']



__date__ = "$Nov 13, 2009 4:28:05 PM$"


