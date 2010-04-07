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

import os
import sys
import signal
import time
from twisted.python import usage
from cassandra.scripts.startup import start

def stop(config, signame="TERM"):
    basedir = config['basedir']
    os.chdir(basedir)
    try:
        f = open("twistd.pid", "rt")
    except:
        print "Daemon is not running"
        
    pid     = int(f.read().strip())
    signum  = getattr(signal, "SIG"+signame)
    timer   = 0
    os.kill(pid, signum)
    time.sleep(0.1)
    while timer < 10:
        # Poll once per second until twistd.pid goes away, up to 10 seconds
        try:
            os.kill(pid, 0)
        except OSError:
            print "Cassandra process %d is dead" % pid
            return      # When pid goes away
        timer += 1
        time.sleep(1)   # Sleep for 1 sec
    print "Never saw process go away"


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
        start(so)
    elif command == "stop":
        stop(so)


__date__ = "$Mar 8, 2010 5:19:13 PM$"


