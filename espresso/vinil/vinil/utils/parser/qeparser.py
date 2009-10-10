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

textPW = """
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
NAMELIST_NAMES  = [ "control", "system", "electrons"] 

# Regular expressions
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

import re
from vinil.utils.orderedDict import OrderedDict
from namelist import Namelist
from card import Card

class QEParser:
    def __init__(self, filename=None, configText=None, type='pw'):
        self.namelists  = OrderedDict()
        self.cards      = OrderedDict()
        self.filename   = filename
        self.configText = configText

    def parse(self):
        text = self.configText
        self._parseNamelists(text)
        self._parseCards(text)
        return (self.namelists, self.cards)

    def _parseNamelists(self, text):
        namelists  = OrderedDict()
        p   = re.compile(COMMENT)
        s1  = re.sub(p, '', text)           # Remove comments
        p2  = re.compile(NAMELIST)
        matches     = p2.findall(s1)        # Finds all namelist blocks
        for m in matches:
            name    = m[0].lower()
            if name in NAMELIST_NAMES:
                params  = self._parseParams(m[1])     # Parse parameters from a namelist block
                namelists[name.lower()] = params

        self._convertNamelists(namelists)

    # Converts from dictionary to Namelist
    def _convertNamelists(self, namelists):
        for name in namelists.keys():
            nl      = Namelist(name)
            for p in namelists[name]:
                nl.add(p[0], p[1])
                
            self.namelists[name] = nl

#        for n in self.namelists.keys():
#            print self.namelists[n].toString()

    # Parses parameters
    def _parseParams(self, text):
        params  = []
        p   = re.compile(EXPRESSION)        # Match expression
        matches = p.findall(text)
        for m in matches:
            pl  = self._getParams(m)         # Parameters list
            params.append(pl)

        return params

    def _getParams(self, text):
        """ Takes string like 'a = 2' and returns tuple ('a', 2) """
        s = text.split('=')
        for i in range(len(s)):
            s[i] = s[i].strip()

        param   = s[0]
        val     = s[1]
        # Assume that there are two values only: (variable, value) pair
        assert len(s) == 2

        return (param, val)

    def _parseCards(self, text):
        p   = re.compile(COMMENT)
        s1  = re.sub(p, '', text)       # Remove comments
        p2  = re.compile(NAMELIST)
        s2  = re.sub(p2, '', s1)
        rawlist = []

        p   = re.compile(EMPTY_LINE)
        s   = s2.split('\n')
        for line in s:
            line    = line.strip()
            if line != '':
                rawlist.append(line)

        self._convertCards(self._getCards(rawlist))

    def _getCards(self, rawlist):
        cards       = OrderedDict()
        cardName    = None
        for l in rawlist:
            firstPart   = l.split()[0].lower()
            if firstPart in CARD_NAMES:
                cardName    = l.lower()
                cards[cardName]    = []
            elif cardName is not None:
                cards[cardName].append(l)

        return cards

    def _convertCards(self, cards):
        for cname in cards.keys():
            c   = Card(cname)
            for l in cards[cname]:
                c.addLine(l)

            self.cards[cname]    = c

        for c in self.cards.keys():
            print self.cards[c].toString()



if __name__ == "__main__":
    qeparser    = QEParser(configText = textPW)
    qeparser.parse()

__date__ = "$Oct 9, 2009 4:34:28 PM$"


