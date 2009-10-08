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
Parsing steps:
1. Remove comments!
2. Remove empty new lines, 
"""

import re

textA   = """hi ! b?#@(
blah ! another comment"""


textSimple = """ 	&system ! blah
    ibrav=2, ! another ! comment
    celldm(1) =6.65,
/
! comment
"""

textB = """ & system # ! comment
    ibrav=2, ! fignja vsja_kaja!
    celldm(1) =6.65, ! tozhe fignja
/ blah!
"""

comment     = '!.*'
namelist    = '([a-zA-Z_]*).*'    # Extracts namelist name ()
spaces      = '[ \t]*'          # Spaces and tabs
newline     = '[\n\r]*'         # New line ()

def parser(s):
    p   = re.compile(comment)
    s1  = re.sub(p, '', s)      # Remove comments
    patstr  = '%s&%s%s([\s\w,.=()]*)/' % (spaces, spaces, namelist)
    p2  = re.compile(patstr)
    m   = p2.match(s1)

    s2  = m.group(2)            # Namelist content
    s3  = re.sub(newline, '', s2)   # Remove empty new lines


    print s3


if __name__ == "__main__":
    parser(textB)

__date__ = "$Oct 8, 2009 11:52:34 AM$"


