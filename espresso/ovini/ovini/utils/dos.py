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

from subprocess import Popen,  PIPE

def runPW_DOS(infile=None, config=None):
    # If infile is None use database
    # else parse infile
    if infile is not None:
        f = open(infile)
        buf = f.read() # buf is 'string'
        f.close()

    if config is not None:
        buf = config

    proc = Popen("dos.x", stdin=PIPE, stdout=PIPE,  stderr=PIPE,  shell="/bin/bash")
    (stdout,   stderr) = proc.communicate(buf) # stdout is 'string'
    return (stdout,   stderr)

__date__ = "$Jul 30, 2009 12:08:03 PM$"


