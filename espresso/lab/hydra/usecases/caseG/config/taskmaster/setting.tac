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
from buildbot.master import BuildMaster

basedir = r'/home/dexity/danse-workspace/AbInitio/espresso/lab/hydra/usecases/caseG/config/taskmaster/'
configfile = r'taskmaster.cfg'
rotateLength = 1000000
maxRotatedFiles = None

application = service.Application('taskmaster')
try:
    from twisted.python.logfile import LogFile
    from twisted.python.log import ILogObserver, FileLogObserver
    logfile = LogFile.fromFullPath("twistd.log",
                                    rotateLength    = rotateLength,
                                    maxRotatedFiles = maxRotatedFiles)
    application.setComponent(ILogObserver, FileLogObserver(logfile).emit)
except ImportError:
    # probably not yet twisted 8.2.0 and beyond, can't set log yet
    pass

BuildMaster(basedir, configfile).setServiceParent(application)


__date__ = "$Apr 6, 2010 10:57:59 AM$"


