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

JT_PORT = 8021

from twisted.protocols import basic

class JobTracker(basic.LineReceiver):
    def connectionMade(self):
        print "Got new client!"
        self.factory.clients.append(self)

    def connectionLost(self, reason):
        print "Lost a client!"
        self.factory.clients.remove(self)

    def lineReceived(self, line):
        print "received", repr(line)
        for c in self.factory.clients:
            c.message(line)

    def message(self, message):
        self.transport.write(message + '\n')


class JobClientFactory(ClientFactory):
    protocol = JobTracker

    def clientConnectionFailed(self, connector, reason):
        print 'connection failed:', reason.getErrorMessage()
        reactor.stop()

    def clientConnectionLost(self, connector, reason):
        print 'connection lost:', reason.getErrorMessage()
        reactor.stop()



from twisted.internet import protocol
from twisted.application import service, internet

factory = protocol.ServerFactory()
factory.protocol = JobTracker
factory.clients = []

application = service.Application("jobtracker")
internet.TCPServer(JT_PORT, factory).setServiceParent(application)


__date__ = "$Feb 28, 2010 9:44:36 PM$"


