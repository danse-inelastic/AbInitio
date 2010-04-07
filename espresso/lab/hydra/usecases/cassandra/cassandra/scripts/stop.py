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

import os
import signal
import time

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


__date__ = "$Apr 6, 2010 6:14:49 PM$"


