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

import re

text = """
 &control
    pseudo_dir = '/home/markovsk/projects/pslib/espresso/mgb2/'
    outdir='temp/'
   /
&electrons
conv_thr = 1.0d-12
diago_full_acc=.TRUE.
/

"""

line = "blah '/home/markovsk/projects/pslib/espresso/mgb2/'"

NAME            = '([a-zA-Z_]*)[^/]'    # Extracts namelist name ()
SPACES          = '[ \t]*'              # Spaces and tabs
NAMELIST        = """%s&%s%s(([^/]*|(['"][^'"]*['"])?)*)/""" % (SPACES, SPACES, NAME)        # Namelist block # [^/]
PATH            = """([^/]*(['"][^'"]*['"]))*"""    #  # ([^/]*(['"][^'"]*['"]))*

def testSlash(text):
    p   = re.compile(NAMELIST)  # PATH)  #
    m   = p.match(text)
    if m is not None:
        print m.group(0)

if __name__ == "__main__":
    testSlash(text)

__date__ = "$Oct 20, 2009 11:04:00 AM$"


