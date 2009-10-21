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

# The namelist block is marked by

line = "blah '/home/markovsk/projects/pslib/espresso/mgb2/' foo 'temp/'/"

NAME            = '([a-zA-Z_]*)[^/]'    # Extracts namelist name ()
SPACES          = '[ \t]*'              # Spaces and tabs
NAMELIST        = """%s&%s%s(([^/]*(['"][^'"]*['"])?)*)/""" % (SPACES, SPACES, NAME)        # Namelist block # [^/]
DIR             = """['"]([^/'"]*/[^/'"]*)*['"]"""
PATH            = """%s&%s%s([^&]*)/""" % (SPACES, SPACES, NAME)   #

def testSlash(text):
    p   = re.compile(PATH)  # NAMELIST)  #
    m   = p.findall(text)
    if m is not None:
        for n in m:
            print n[1]

if __name__ == "__main__":
    testSlash(text)

__date__ = "$Oct 20, 2009 11:04:00 AM$"

# """%s&%s%s(([^/]*%s[^/]*)*)[^'"]*/"""
# ([^/]*(['"][^'"]*['"]))*

