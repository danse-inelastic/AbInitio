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

from twisted.protocols import basic

class TaskTracker(basic.LineReceiver):
    def connectionMade(self):
        print "Got new client!"
        self.factory.clients.append(self)

    def connectionLost(self, reason):
        print "Lost a client!"
        self.factory.clients.remove(self)

    def lineReceived(self, line):
        print "received", repr(line)
        if len(self.factory.clients) > 0:
            c   = self.factory.clients[0]
            c.message(line)

    def message(self, message):
        self.transport.write(message + ' TaskTracker \n')


from twisted.internet import protocol
from twisted.application import service, internet

factory = protocol.ServerFactory()
factory.protocol = TaskTracker
factory.clients = []

application = service.Application("tasktracker")
internet.TCPServer(TT_PORT, factory).setServiceParent(application)

__date__ = "$Feb 28, 2010 9:44:22 PM$"


