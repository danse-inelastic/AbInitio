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

from twisted.application import service
from twisted.python.logfile import LogFile
from twisted.python.log import ILogObserver, FileLogObserver
from cassandra.applications.Worker import Worker

basedir     = r'/home/dexity/danse-workspace/AbInitio/espresso/lab/hydra/usecases/cassandra/config/worker'
configfile  = r'worker.cfg'

application = service.Application('worker')
logfile     = LogFile.fromFullPath("twistd.log")
application.setComponent(ILogObserver, FileLogObserver(logfile).emit)

Worker(basedir, configfile).setServiceParent(application)

__date__ = "$Apr 6, 2010 10:57:59 AM$"


