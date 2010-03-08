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
from twisted.internet.protocol import ClientFactory #, ClientCreator
from twisted.internet import reactor

# The mechanism of including client to the JobTracker is borrowed from twisted.web.proxy.Proxy
# which is indirect subclass of twisted.protocols.basic.LineReceiver

class JobTracker(LineReceiver):

    def connectionMade(self):
        print "JobTracker: got new client!"
        
        self.clientFactory     = JobClientFactory()    #self)
        self.clientFactory.client   = None
        reactor.connectTCP('localhost', TT_PORT, self.clientFactory)
        #self.client = self.fc.protocol() #self.fc.client

        #print "self.client " + self.client
        #self.factory.clients.append(self)


    def connectionLost(self, reason):
        print "JobTracker: lost a client!"
        #self.factory.clients.remove(self)


    def lineReceived(self, line):
        print "JobTracker: received", repr(line)
        message     = line + ' JobTracker '

        print dir(self.clientFactory.client)
        #print dir(self.clientFactory)
        #client      = self.clientFactory.protocol
        #print dir(client.transport)
#        client.makeConnection(client)
        #print client.transport
        #print dir(client)
        client.send(message)


#        if len(self.factory.clients) == 1:
#            try:
#                self.transport.write(message + 'try')
#
##                factory = JobClientFactory(self)
##                reactor.connectTCP('localhost', TT_PORT, factory)
#            except:
#                self.transport.write("TaskTracker not found!")
#
#        else:
#        self.transport.write(message)

#        for c in self.factory.clients:
#            c.message(line)
#
#    def message(self, message):
#        self.transport.write(message + ' JobTracker ')


class JobClient(LineReceiver):

    def connectionMade(self):
        print "JobClient: got new client!"
        print dir(self.factory)
        
#        LineReceiver.connectionMade(self)
        #self.factory.client  = self

#        self.sendLine("Mur-mur")

    def lineReceived(self, line):
        print "JobClient: received", repr(line)


    def send(self, line):
        print "send"
        #print self.transport
        self.sendLine(line)
        
#        try:
#            factory = JobClientFactory()
#            reactor.connectTCP('localhost', TT_PORT, factory)
#        except:
#            self.factory.clients.remove(self)

class JobClientFactory(ClientFactory):
    protocol       = JobClient

#    def __init__(self, proxy):
#        #print JobClient
#        self.protocol       = JobClient
#        #self.protocol.proxy = proxy

    def clientConnectionFailed(self, connector, reason):
        print 'connection failed:', reason.getErrorMessage()
        reactor.stop()

    def clientConnectionLost(self, connector, reason):
        print 'connection lost:', reason.getErrorMessage()
        reactor.stop()

#    def buildProtocol(self, addr):
#        #protocol       = JobClient
#        return self.protocol()


#        client = ClientCreator(reactor, JobClient)
#        client.connectTCP('localhost', TT_PORT)



__date__ = "$Feb 28, 2010 9:44:36 PM$"


