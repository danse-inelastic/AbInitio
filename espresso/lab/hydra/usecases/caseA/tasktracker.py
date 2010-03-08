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

from twisted.protocols import basic

class TaskTracker(basic.LineReceiver):
    def connectionMade(self):
        print "TaskTracker: got new client!"


    def connectionLost(self, reason):
        print "TaskTracker: lost a client!"


    def lineReceived(self, line):
        print "received", repr(line)
        self.transport.write(line + ' TaskTracker ')


#        if len(self.factory.clients) > 0:
#            c   = self.factory.clients[0]
#            c.message(line)
#
#    def message(self, message):
        
        #self.factory.clients.remove(self)
        
        #self.factory.clients.append(self)

#

__date__ = "$Feb 28, 2010 9:44:22 PM$"


