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
Generates python dictionary of QE pseudo-potentials

Dictionary of Pseudo Potentials generated from espresso_pp.tar.gz:
http://www.quantumespresso.org/pseudo.php . Requires file 'pseudo.list' with
list of pseudo-potentials as in example:

[pseudo.list]
Ag.pbe-d-rrkjus.info
Ag.pbe-d-rrkjus.UPF
Ag.pz-d-rrkjus.UPF
Al.blyp-n-van_ak.info
...

"""

from orderedDict import OrderedDict

def testPseudo():
    from pseudotest import PSEUDO
    return PSEUDO


OUTPUT  = "pseudo.py"
INPUT   = "pseudo.list"
SPACE   = "                "
HEADER  = """#!/usr/bin/env python
from orderedDict import OrderedDict

PSEUDO = OrderedDict()
"""

def toString(pseudo):
    s = HEADER
    for name in pseudo.keys():
        s   += 'PSEUDO["%s"] = (\n' % name
        for item in pseudo[name]:
            s   += '%s "%s",\n' % (SPACE, item)
        s   += '%s)\n\n' % SPACE

    return s

def generate():
    inp    = open("pseudo.list")
    #Generate structure
    pseudo = OrderedDict()
    out    = open("pseudo.py", "w")
    #out.write(toString(pseudo))
    out.write(toString(testPseudo()))
    out.close()

if __name__ == "__main__":
    generate()

__date__ = "$Jan 19, 2010 6:21:45 AM$"


