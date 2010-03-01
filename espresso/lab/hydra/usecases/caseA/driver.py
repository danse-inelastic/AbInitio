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

from twisted.internet import protocol
from twisted.application import service, internet
from jobtracker import JobTracker

jtFactory = protocol.ServerFactory()
jtFactory.protocol = JobTracker
jtFactory.clients = []

application = service.Application("BadPhone")

jtService   = internet.TCPServer(JT_PORT, jtFactory)
jtService.setName("JobTracker")
jtService.setServiceParent(application)


__date__ = "$Mar 1, 2010 12:57:11 AM$"


