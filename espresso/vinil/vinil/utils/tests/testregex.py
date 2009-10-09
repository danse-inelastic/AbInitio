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

& systemB # ! comment
    ibrav=2, ! fignja vsja_kaja!
    celldm(1) =6.65, ! tozhe fignja
/ blah!
"""

textC = """
 &control
    calculation='scf'
    restart_mode='from_scratch',
    tprnfor = .true.
    prefix='ni',
    pseudo_dir = '',
    outdir=''
 /
 &system
    ibrav=2,
    celldm(1) =6.65,
    nat=  1,
    ntyp= 1,
    nspin=2,
    starting_magnetization(1)=0.5,
    degauss=0.02,
    smearing='gauss',
    occupations='smearing',
    ecutwfc =27.0
    ecutrho =300.0
 /
 &electrons
    conv_thr =  1.0d-8
    mixing_beta = 0.7
 /
"""


COMMENT     = '!.*'                 # Comment
NAME        = '([a-zA-Z_]*)[^/]'    # Extracts namelist name ()
SPACES      = '[ \t]*'              # Spaces and tabs
NO_SPACES   = '[^\s]*'              # No spaces
NEWLINE     = '[\n\r]*'             # New line ()
PARAMTER    = '[\w,()]+'            # Parameter characters (space is not allowed)
VALUE       = '[^\s,]+'              # Parameter's value (numerate all possible characters)
EXPRESSION  = '(%s%s=%s%s)' % (PARAMTER, SPACES, SPACES, VALUE) # Parameter's expression
NAMELIST    = """%s&%s%s([^/]*)/""" % (SPACES, SPACES, NAME) # Namelist block

def parser(text):
    nl  = parseNamelists(text)
    parseCards(text)
    print nl

def parseNamelists(text):
    namelists = []
    p   = re.compile(COMMENT)
    s1  = re.sub(p, '', text)       # Remove comments
    p2  = re.compile(NAMELIST)
    matches   = p2.findall(s1) #match(s1)
    for m in matches:
        name    = m[0]
        params  = parseParams(m[1])
        namelists.append((name, params))

    return namelists

def parseParams(text):
    params  = []
    p   = re.compile(EXPRESSION)    # Match expression
    matches = p.findall(text)
    for m in matches:
        pl  = getParams(m) # Parameters list
        params.append(pl)

    return params

def getParams(text):
    """ Takes string like 'a = 2' and returns tuple ('a', 2) """
    s = text.split('=')
    for i in range(len(s)):
        s[i] = s[i].strip()

    param   = s[0]
    val     = s[1]
    # Assume that there are two values only: (variable, value) pair
    assert len(s) == 2
    
    return (param, val)

def parseCards(text):
    #print "Parsing Cards"
    return

if __name__ == "__main__":
    parser(textC)

__date__ = "$Oct 8, 2009 11:52:34 AM$"


