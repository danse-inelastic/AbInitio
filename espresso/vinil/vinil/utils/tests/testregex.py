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

textA   = """hi ! b?#@(
blah ! antoher comment""" #"&car ! some comment"

textSimple = """ &system
    ibrav=2,
    celldm(1) =6.65,
/
"""

def parser(text, pattern):
    p = re.compile(pattern)
    s = p.sub('', text)
    #p.sub()
    #m = p.match(text)
    #print m.span()
    #g = m.group(1)
    #g.sub('', )
    print s


if __name__ == "__main__":
    p   = '!.*'
    parser(textA, p)

__date__ = "$Oct 8, 2009 11:52:34 AM$"


