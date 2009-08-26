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

# For now it's just a function
def runPWSimulation(outfile, infile=None):
    # If infile is None use database
    # else parse infile
    f = open(infile)
    buf = f.read() # buf is 'string'
    f.close()

    proc = Popen("pw.x", stdin=PIPE, stdout=PIPE,  stderr=PIPE,  shell="/bin/bash")
    (stdout,   stderr) = proc.communicate(buf) # stdout is 'string'
    out = open(outfile,  "w")
    out.write(stdout)
    out.close()

__date__ = "$Jul 30, 2009 12:07:54 PM$"


