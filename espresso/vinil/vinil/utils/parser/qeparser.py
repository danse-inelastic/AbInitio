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

class QEParser:
    def __init__(self, filename=None, configText=None, type='inputpw'):
        pass

    def parser(self, text):
        namelists = self.parseNamelists(text)
        cards     = self.parseCards(text)
        print (namelists, cards)

    def parseNamelists(self, text):
        namelists = OrderedDict()
        p   = re.compile(COMMENT)
        s1  = re.sub(p, '', text)       # Remove comments
        p2  = re.compile(NAMELIST)
        matches     = p2.findall(s1)      # Finds all namelist blocks
        for m in matches:
            name    = m[0]
            if name in NAMELIST_NAMES:
                params  = self.parseParams(m[1])     # Parse parameters from a namelist block
                namelists[name.lower()] = params

        #noneparsed = re.sub(p2, '', s1)
        return namelists #(namelists, noneparsed)

    def parseParams(self, text):
        params  = []
        p   = re.compile(EXPRESSION)    # Match expression
        matches = p.findall(text)
        for m in matches:
            pl  = getParams(m)          # Parameters list
            params.append(pl)

        return params

    def getParams(self, text):
        """ Takes string like 'a = 2' and returns tuple ('a', 2) """
        s = text.split('=')
        for i in range(len(s)):
            s[i] = s[i].strip()

        param   = s[0]
        val     = s[1]
        # Assume that there are two values only: (variable, value) pair
        assert len(s) == 2

        return (param, val)

    def parseCards(self, text):
        p   = re.compile(COMMENT)
        s1  = re.sub(p, '', text)       # Remove comments
        p2  = re.compile(NAMELIST)
        s2  = re.sub(p2, '', s1)
        rawlist = []
        #cards   = []
        p   = re.compile(EMPTY_LINE)
        s   = s2.split('\n')
        for line in s:
            line    = line.strip()
            if line != '':
                rawlist.append(line)

        cards   = self.getCards(rawlist)
        return cards

    def getCards(self, rawlist):
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
    pass

__date__ = "$Oct 9, 2009 4:34:28 PM$"


