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

 
ATOMIC_SPECIES 
 Ni  26.98  Ni.pbe-nd-rrkjus.UPF
 
ATOMIC_POSITIONS
 Ni 0.00 0.00 0.00
K_POINTS AUTOMATIC
4 4 4 1 1 1
blah

"""

CARD_NAMES      = ["atomic_species", "atomic_positions", "k_points"]
NAMELIST_NAMES  = [ "system", "electrons"] # "control",

COMMENT     = '!.*'                 # Comment
NAME        = '([a-zA-Z_]*)[^/]'    # Extracts namelist name ()
SPACES      = '[ \t]*'              # Spaces and tabs
NO_SPACES   = '[^\s]*'              # No spaces
NEWLINE     = '[\n\r]*'             # New line ()
PARAMTER    = '[\w,()]+'            # Parameter characters (space is not allowed)
VALUE       = '[^\s,]+'             # Parameter's value (numerate all possible characters)
EXPRESSION  = '(%s%s=%s%s)' % (PARAMTER, SPACES, SPACES, VALUE)     # Parameter's expression
NAMELIST    = """%s&%s%s([^/]*)/""" % (SPACES, SPACES, NAME)        # Namelist block
CARD        = '(%s[\w]+)' % (SPACES)
EMPTY_LINE  = r'^\s*'                # Empty line

from vinil.utils.orderedDict import OrderedDict

def parser(text):
    (namelists, cardsText)  = parseNamelists(text)
    cards                   = parseCards(cardsText)
    print (namelists, cards)

def parseNamelists(text):
    namelists = OrderedDict()
    p   = re.compile(COMMENT)
    s1  = re.sub(p, '', text)       # Remove comments
    p2  = re.compile(NAMELIST)
    matches     = p2.findall(s1)      # Finds all namelist blocks
    for m in matches:
        name    = m[0]
        if name in NAMELIST_NAMES:
            params  = parseParams(m[1])     # Parse parameters from a namelist block
            namelists[name.lower()] = params
        
    noneparsed = re.sub(p2, '', s1)
    return (namelists, noneparsed)

def parseParams(text):
    params  = []
    p   = re.compile(EXPRESSION)    # Match expression
    matches = p.findall(text)
    for m in matches:
        pl  = getParams(m)          # Parameters list
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
    rawlist = []
    #cards   = []
    p   = re.compile(EMPTY_LINE)
    s   = text.split('\n')
    for line in s:
        line    = line.strip()
        if line != '':
            rawlist.append(line)

    cards   = getCards(rawlist)
    return cards

def getCards(rawlist):
    cards       = OrderedDict()
    cardName    = None
    for l in rawlist:
        firstPart   = l.split()[0].lower()
        if firstPart in CARD_NAMES:
            cardName    = firstPart
            cards[cardName]    = []
        elif cardName is not None:
            cards[cardName].append(l)

    return cards

if __name__ == "__main__":
    parser(textC)

__date__ = "$Oct 8, 2009 11:52:34 AM$"


