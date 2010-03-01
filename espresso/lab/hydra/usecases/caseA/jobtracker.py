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

TT_PORT = 8020

from twisted.protocols.basic import LineReceiver
from twisted.internet.protocol import ClientFactory
from twisted.internet import reactor

class JobTracker(LineReceiver):
    def connectionMade(self):
        print "Got new client!"
        self.factory.clients.append(self)


    def connectionLost(self, reason):
        print "Lost a client!"
        self.factory.clients.remove(self)

    def lineReceived(self, line):
        print "received", repr(line)
        if len(self.factory.clients) == 1:
            try:
                factory = JobClientFactory()
                reactor.connectTCP('localhost', TT_PORT, factory)
            except:
                pass


        for c in self.factory.clients:
            c.message(line)

    def message(self, message):
        self.transport.write(message + ' JobTracker ')


class JobClientFactory(ClientFactory):
    protocol = JobTracker

    def clientConnectionFailed(self, connector, reason):
        print 'connection failed:', reason.getErrorMessage()
        reactor.stop()

    def clientConnectionLost(self, connector, reason):
        print 'connection lost:', reason.getErrorMessage()
        reactor.stop()


#        try:
#            factory = JobClientFactory()
#            reactor.connectTCP('localhost', TT_PORT, factory)
#        except:
#            self.factory.clients.remove(self)



__date__ = "$Feb 28, 2010 9:44:36 PM$"


