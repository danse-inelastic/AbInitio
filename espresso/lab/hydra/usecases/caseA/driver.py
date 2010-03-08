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
TT_PORT = 8020

from twisted.internet import protocol
from twisted.application import service, internet
from jobtracker import JobTracker
from tasktracker import TaskTracker

application = service.Application("BrokenPhone")

# JobTracker
ttFactory = protocol.Factory()
ttFactory.protocol = JobTracker
#ttFactory.clients = []

ttService   = internet.TCPServer(JT_PORT, ttFactory)
ttService.setName("JobTracker")
ttService.setServiceParent(application)

# TaskTracker
ttFactory = protocol.ServerFactory()
ttFactory.protocol = TaskTracker
#ttFactory.clients = []

ttService   = internet.TCPServer(TT_PORT, ttFactory)
ttService.setName("TaskTracker")
ttService.setServiceParent(application)


__date__ = "$Mar 1, 2010 12:57:11 AM$"


