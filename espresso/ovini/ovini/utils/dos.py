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

def run_pw_dos(infile):
    f = open(infile)
    buf = f.read() # buf is 'string'
    f.close()
    proc = Popen("bin/dos.x", stdin=PIPE, stdout=PIPE,  stderr=PIPE,  shell="/bin/bash")
    (stdout,   stderr) = proc.communicate(buf) # stdout is 'string'

__date__ = "$Jul 30, 2009 12:08:03 PM$"


