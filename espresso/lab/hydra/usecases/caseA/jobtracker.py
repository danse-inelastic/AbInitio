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

JOB_MESSAGE = " -> JobTracker"

from ports import TT_PORT
from twisted.protocols.basic import LineReceiver
from twisted.internet.protocol import ClientFactory
from twisted.internet import reactor

"""
Starts as a daemon. Receives the line from the client, attaches job message and
sends to the TaskTracker, receives line from the TaskTracker and sends back to
the client

Notes:
    1. The mechanism of including client to the JobTracker is borrowed from
    twisted.web.proxy.Proxy which is indirect subclass of twisted.protocols.basic.LineReceiver
"""


class JobTracker(LineReceiver):

    def connectionMade(self):
        print "JobTracker: got new client!"
        
        self.clientFactory     = JobClientFactory()
        self.clientFactory.client   = None
        self.clientFactory.proxy    = self      # Important line!
        reactor.connectTCP('localhost', TT_PORT, self.clientFactory)


    def connectionLost(self, reason):
        print "JobTracker: lost a client!"
        self.factory.client = None


    def lineReceived(self, line):
        print "JobTracker: received", repr(line)
        message     = line + JOB_MESSAGE

        client      = self.clientFactory.client
        client.sendLine(message)


class JobClient(LineReceiver):

    def connectionMade(self):
        print "JobClient: got new client!"
        self.factory.client   = self    # Save the client


    def lineReceived(self, line):
        print "JobClient: received", repr(line)
        proxy   = self.factory.proxy
        proxy.sendLine(line + JOB_MESSAGE)


class JobClientFactory(ClientFactory):
    protocol       = JobClient

    def clientConnectionFailed(self, connector, reason):
        print 'connection failed:', reason.getErrorMessage()
        reactor.stop()

    def clientConnectionLost(self, connector, reason):
        print 'connection lost:', reason.getErrorMessage()
        reactor.stop()



__date__ = "$Feb 28, 2010 9:44:36 PM$"


