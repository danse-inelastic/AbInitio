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

# Should I inherit from pyre.services.Session instead?

from pyre.services.Session import Session as base

class EchoSession(base):

    class Inventory(base.Inventory):
        import pyre.inventory

        from pyre.services import pickler
        marshaller = pyre.inventory.facility("marshaller", factory=pickler)


    def __init__(self, name):
        super(EchoSession, self).__init__(name)
        self.client = None # set by our parent
        return


    def _configure(self):
        super(EchoSession, self)._configure()
        self.marshaller = self.inventory.marshaller
        return



#    from pyre.services.RequestError import RequestError


#    def _init(self):
#        from pyre.ipc.Socket import Socket
#
#        # connect to the server
#        try:
#            self._connect()
#        except Socket.ConnectionError, error:
#            raise self.RequestError(str(error))
#
#        # start the session
#        address = (self.host, self.port)
#        self.session = self.client.newSession(
#            self._connection, address, self.marshaller
#            )
#
#        return
__date__ = "$Nov 13, 2009 4:29:55 PM$"


