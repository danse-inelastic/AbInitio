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

# Twisted specific
from twisted.internet import protocol
from twisted.application import service, internet

# CaseA specific
from ports import JT_PORT, TT_PORT
from jobtracker import JobTracker
from tasktracker import TaskTracker

application = service.Application("ProxyEcho")

# JobTracker
jtFactory = protocol.Factory()
jtFactory.protocol  = JobTracker

jtService   = internet.TCPServer(JT_PORT, jtFactory)
jtService.setName("JobTracker")
jtService.setServiceParent(application)

# TaskTracker
ttFactory = protocol.ServerFactory()
ttFactory.protocol = TaskTracker

ttService   = internet.TCPServer(TT_PORT, ttFactory)
ttService.setName("TaskTracker")
ttService.setServiceParent(application)


__date__ = "$Mar 1, 2010 12:57:11 AM$"


