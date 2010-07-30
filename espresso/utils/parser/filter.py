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
CARD_REQ        = "name"
NAMELIST_KEYS   = ("name", "params")
NAMELIST_REQ    = "name"

# Auxiliary functions
ifelse  = lambda a,b,c: (b,c)[not a]

class Filter(object):

    def __init__(self, name = None):
        """
            name:  (str) -- Name of the filter
        """
        self._name      = name
        self._fnamelists  = []    # Filtered namelists
        self._fcards      = []    # Filtered cards


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
        return self._fnamelists


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

        self._checkDictFormat(card, CARD_KEYS, CARD_REQ)
        
        c       = Card( card["name"],
                        arg = ifelse(card.has_key("arg"), card["arg"], None))
        if card.has_key("lines"):
            c.setLines(card["lines"])
            
        self._fcards.append(c)


    def removeCard(self, name):
        """
        Removes card
        """
        if not name:    # No name, just ignore
            return

        for c in self._fcards:
            if c.name() == name:  # Remove the first card which matches name
                self._fcards.remove(c)


    def cards(self):
        """
        Returns filter cards
        """
        return self._fcards
    

    def apply(self, input):
        """
        Applies filter to the input

            input: (object: QEInput) -- Input object
        """
        pass


    def _checkDictFormat(self, item, keys=None, reqkeys=None ):
        """
        Checks dictionary format

            item: (dict) -- Dictionary which format is being checked
            keys: (list) -- Keys allowed in the dictionary
            reqkeys: (str) -- Required keys in the dictionary, single string for now
        """
        if not type(item) == dict:
            raise TypeError("Parameter '%s' must be dictionary" % str(item))

        for k in item.keys():
            if not k in keys:
                raise KeyError("Invalid key: %s" % k)

        if not item.has_key(reqkeys):
            raise KeyError("No key: '%s'", reqkeys)


__date__ = "$Jul 28, 2010 2:35:31 PM$"


