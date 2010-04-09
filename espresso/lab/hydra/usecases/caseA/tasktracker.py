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

TASK_MESSAGE = " -> TaskTracker"
from twisted.protocols import basic

"""
Starts as a daemon. Receives a line and sends back the attached task message
"""

class TaskTracker(basic.LineReceiver):
    def connectionMade(self):
        print "TaskTracker: got new client!"
        self.factory.client = self  # Save the client


    def connectionLost(self, reason):
        print "TaskTracker: lost a client!"
        self.factory.client = None


    def lineReceived(self, line):
        print "TaskTracker: received", repr(line)
        client  = self.factory.client
        self.sendLine(line + TASK_MESSAGE)


__date__ = "$Feb 28, 2010 9:44:22 PM$"


