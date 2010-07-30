# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2010  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

"""
Filter - class for filtering card and namelist parameters in configuration input
"""

from namelist import Namelist
from card import Card

CARD_KEYS       = ("name", "lines", "arg")
NAMELIST_KEYS   = ("name", "params")

# Auxiliary functions
ifelse  = lambda a,b,c: (b,c)[not a]

class Filter(object):

    def __init__(self, name = None):
        """
            name:  (str) -- Name of the filter
        """
        self._name      = name
        self._namelists  = []    # Filtered namelists
        self._cards      = []    # Filtered cards


    def name(self):
        "Returns name of the filter"
        return self._name


    def setParam(self, namelist, param, value):
        """
        Set parameter of the namelist. If the parameter or namelist do not exist,
        it creates them.

            namelist: (str) -- Name of the namelist
            param: (str) -- Name of the parameter
            value: (str) -- Valur of the parameter
        """


    def removeParam(self, namelist, param):
        """
        Removes parameter from the namelist. If the parameter or namelist do
        not exist it will just ignore.

            namelist: (str) -- Name of the namelist
            param: (str) -- Name of the parameter
        """
        pass


    def setNamelist(self, namelist):
        """
        Sets namelist as dictionaty:

            namelist: (dict) -- Dictionary of namelist

        Format:
        {"name":      <namelist name>,
         "params":    {param1:    value1,
                       param2:    value2,
                       ...}}

        Example:
        {"name":      "control",
         "params":    {calculation:   'scf',
                       restart_mode:  'from_scratch',
                       tprnfor:       .true.}}
        """
        pass


    def removeNamelist(self, name):
        pass


    def namelists(self):
        """
        Returns namelists
        """
        return self._namelists


    def setCard(self, card):
        """
        Set card as dictionary

            card: (dict) -- Dictionary of card

        Format:
        {"name":      <card name>,
         "lines":    (line1,
                      line2,
                       ...),
         "arg":       <argument>}

        Example:
        {"name":      "k_points",
         "lines":     ("4 4 4 1 1 1",),
         "arg":       "automatic"}
        """
        if not card:    # Ignore empty card 
            return

        if not type(card) == dict:
            raise TypeError("Parameter 'card' must be dictionary")

        for k in card.keys():
            if not k in CARD_KEYS:
                raise KeyError("Invalid key: %s" % k)

        if not card.has_key("name"):
            raise KeyError("No key: 'name'")
        
        c       = Card( card["name"],
                        arg = ifelse(card.has_key("arg"), card["arg"], None))
        if card.has_key("lines"):
            c.setLines(card["lines"])
            
        self._cards.append(c)


    def removeCard(self, name):
        """
        Removes card
        """


    def cards(self):
        """
        Returns filter cards
        """
        return self._cards
    

    def apply(self, input):
        """
        Applies filter to the input

            input: (object: QEInput) -- Input object
        """
        pass

__date__ = "$Jul 28, 2010 2:35:31 PM$"


