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

"""
There should be just one binary that launches daemons

"""

import sys
from twisted.python import usage


def stop(so):
    print "Stop "


class MakerBase(usage.Options):
    opt_h = usage.Options.opt_help

    def parseArgs(self, *args):
        if len(args) > 0:
            self['basedir'] = args[0]
        else:
            self['basedir'] = None
        if len(args) > 1:
            raise usage.UsageError("I wasn't expecting so many arguments")


    def postOptions(self):
        if self['basedir'] is None:
            raise usage.UsageError("<basedir> parameter is required")
        #self['basedir'] = os.path.abspath(self['basedir'])


class StartOptions(MakerBase):
    def getSynopsis(self):
        return "Usage:    cassandra start <basedir>"


class StopOptions(MakerBase):
    def getSynopsis(self):
        return "Usage:    cassandra stop <basedir>"


class Options(usage.Options):
    synopsis = "Usage:    cassandra <command> [command options]"
    
    subCommands = [
                    ['start', None, StartOptions, "Start taskmaster or worker daemons"],
                    ['stop', None, StopOptions, "Stop taskmaster or worker daemons"],
                  ]

    def postOptions(self):
        if not hasattr(self, 'subOptions'):
            raise usage.UsageError("Must specify a command")


def run():
    config = Options()
    try:
        config.parseOptions()
    except usage.error, e:
        print "%s:  %s" % (sys.argv[0], e)
        print
        c = getattr(config, 'subOptions', config)
        print str(c)
        sys.exit(1)

    command = config.subCommand     # passed command!
    so = config.subOptions

    if command == "start":
        from scripts.startup import start
        start(so)
    elif command == "stop":
        stop(so)


__date__ = "$Mar 8, 2010 5:19:13 PM$"


