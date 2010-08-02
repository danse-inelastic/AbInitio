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

SPACES      = '[ \t]*'              # Spaces and tabs
OPEN_BRACKET    = '[({]?'           # Open bracket
CLOSE_BRACKET   = '[)}]?'           # Close bracket
CARD        = '(%s[\w]+)%s%s(%s[\w]*%s)%s' % (SPACES, SPACES, OPEN_BRACKET, SPACES, SPACES, CLOSE_BRACKET)

cardRef   = ('atomic_positions', 'k_points')

import re

def testCards(rawlist):
    cards   = {}
    cardName    = None
    for l in rawlist:
        p   = re.compile(CARD)
        m   = p.match(l)
        firstPart   = m.group(1).lower()
        if firstPart in cardRef:
            secondPart  = m.group(2).strip().lower()    # Catch argument of the card
            cardName    = firstPart
            cards[cardName]    = {}
            if (secondPart != ''):
                cards[cardName]["args"] = secondPart
            else:
                cards[cardName]["args"] = None
                
            cards[cardName]["values"]   = []
        elif cardName is not None:
            cards[cardName]["values"].append(l)

    print cards

if __name__ == "__main__":
    rawlist = ['ATOMIC_POSITIONS (alat)', 'Ni 0.00 0.00 0.00', 
               'K_POINTS AUTOMATIC', '4 4 4 1 1 1', 'blah']

    cardA   = 'ATOMIC_POSITIONS (alat)'
    line    = 'Ni 0.00 0.00 0.00'

    testCards(rawlist)
    
    # Output: {'atomic_positions': {'args': 'alat', 'values': ['Ni 0.00 0.00 0.00']}, 'k_points': {'args': 'automatic', 'values': ['4 4 4 1 1 1', 'blah']}}

__date__ = "$Oct 19, 2009 2:54:35 PM$"


