#!/usr/bin/env python

"""
QEInput - parses, modifies and creates Quantum Espresso (QE) configuration files

Features:
    - Parses existing configuration file
    - Adds, edits and removes parameters from/to a namelist or a card
    - Creates new configuration file

Parameters of QE configuration files are defined in the Doc directory of the
source code. You can download it from:

    http://www.quantum-espresso.org/download.php


Main Use Cases:
    - Parse existing configuration file, modify parameters in namelists or cards
    and save updated configuration to the same file.
    - Create new configuration file from scratch

Stability Issues:
    - Card starts with card title on a separate line and values between card titles.

Example:

ATOMIC_SPECIES
 Al  26.9815 Al.blyp-n-van_ak.UPF

ATOMIC_POSITIONS (crystal)
 Al      0.00000000  0.00000000  0.00000000


Implementation Issues:
    - Saves both Namelists and Cards titles in capital
    - Refactoring?  Introduce class relation: Namelist(Block), Card(Block)
"""

from orderedDict import OrderedDict
from namelist import Namelist
from card import Card
from qeparser import QEParser

"""
Supported types of the configuration file:

type =
    'pw'               - (default)
    'ph'               -
    'pp'               -
    'bands'            -
    'cp'               - 
    'cppp'             -
    'd3'               -
    'dos'              -
    'dynmat'           -
    'initial_state'    -
    'gipaw'            -
    'd1'               -
    'matdyn'           -
    'projwfc'          -
    'pwcond'           -
    'q2r'              -
"""

class QEInput(object):

    def __init__(self, filename=None, config=None, type='pw'):
        """
        Initializes QEInput by passing either filename or config (not both)
        parameters
        
        Parameters:
            filename:   str
                Absolute or relative filename of file to be parsed
            config:     str
                Configuration text to be parsed
            type:       str
                Type of the simulation
        
        """

        self.header         = None
        self.filename       = filename
        self.config         = config
        self.parser         = QEParser(filename, config, type)
        self.type           = type
        self.namelists      = OrderedDict()
        self.cards          = OrderedDict()
        self.attach         = None          # Specific for 'matdyn', 'dynmat', etc.
        self.namelistRef    = None
        self.cardRef        = None
        self.qe             = [self.header, self.namelists, self.cards, self.attach]
        self.filters        = []


    def parse(self):
        """
        Parses the configuration file and stores the values in qe dictionary
        """
        (self.header, self.namelists, self.cards, self.attach) = self.parser.parse()


    def createNamelist(self, name):
        """
        Creates namelist and adds to self.namelists. Return namelist if it exists.

        Parameters:
            name:       str
                Name of the new namelist
        """        
        # If namelist with the name exists, then return it
        if self.namelistExists(name):
            return  self.namelists[name]

        # Otherwise create a new namelist
        name    = name.lower()
        nl  = Namelist(name)
        self.namelists[name] = nl
        return nl


    def removeNamelist(self, name):
        """
        Remove namelist if it exists or ignore otherwise

        Parameters:
            name:       str
                Name of the existing namelist
        """
        name    = name.lower()
        try:
            del(self.namelists[name])
        except KeyError:    # parameter is not present
            return


    def addNamelist(self, namelist):
        """
        Adds namelist to self.namelists

        Parameters:
            namelist:   object
                Namelist object
        """
        if not namelist:    # No namelist, just ignore it!
            return

        self.namelists[namelist.name()] = namelist


    def namelist(self, name):
        """
        Returns namelist specified by name if exists or create a new one

        Parameters:
            name:       str
                Name of the namelist
        """
        if self.namelistExists(name):   # If exists, return namelist
            return self.namelists[name]

        return self.createNamelist(name) # otherwise create a new one


    def namelistExists(self, name):
        """
        Checks if namelist specified by name (lowered before) exists

        Parameters:
            name:       str
                Name of the namelist        
        """
        return self._exists(name, self.namelists.keys())


    def createCard(self, name):
        """
        Creates card and adds to self.cards
        
        Parameters:
            name:       str
                Name of the card        
        """
        card    = Card(name)
        self.cards[name] = card
        return card


    def addCard(self, card):
        """
        Adds card to self.cards
        
        Parameters:
            card:       object
                Card object
        """
        self.cards[card.name()] = card


    def removeCard(self, name):
        """
        Remove card if it exists or ignore otherwise

        Parameters:
            name:       str
                Name of the existing card
        """
        name    = name.lower()
        try:
            del(self.cards[name])
        except KeyError:    # parameter is not present
            return


    def attach(self):
        "Returns attachment"
        return self.attach


    def addAttach(self, text):
        """
        Sets attachment to text. If attachment is not None it still will be
        overwritten

        Parameters:
            text:       str
                Attachment text, usually is appended to the end of the
                configuration file
        """
        self.attach = text


    def removeAttach(self):
        "Sets attachment to None"
        self.attach = None


    def card(self, name):
        """
        Returns card specified by name if exists or create a new one
        """
        if self.cardExists(name):        # If exists, return card
            return self.cards[name]

        return self.createCard(name)


    def cardExists(self, name):
        "Checks if card specified by name exists"
        return self._exists(name, self.cards.keys())


    def getObject(self, name, dict):
        """
        Returns object specified by name
        
        Parameters:
            name:   str
                Name that identifies an object through name() method
            dict:   dict
                Dictionary that stores objects as {name: object}
        """
        if not name or not dict:
            return None

        for n in dict.values():
            if n.name() == name:
                return dict[name]

        return None


    def save(self, filename=None):
        """
        Saves the QEInput structure to the configuration file

        Parameters:
            filename:       str
                File name where the configuration is stored to
        """
        default = "config.out"

        if filename is None:
            if self.filename is not None:
                filename = self.filename
            else:
                filename = default

        f = open(filename, "w")
        f.write(self.toString())
        f.close()


    def type(self):
        "Returns type of the configuration file"
        return self.type


    def applyFilter(self, filter):
        """
        Applies filter to the QEInput

        Parameters:
            filter:     object (Filter)
                Filter that applies changes to QEInput
        """
        if not filter:  # No filter, just ignore it
            return None

        filter.apply(self)  # Apply filter to QEInput
        self.filters.append(filter)


    def toString(self):
        """
        Dumps QEInput structure to string
        """
        (self.namelistRef, self.cardRef)    = self.parser.setReferences()
        s = ''
        if self.header:             # Add header
            s   += "%s\n" % self.header

        for name in self.namelistRef:   # Add namelists
            nl  = self.getObject(name, self.namelists)
            if nl is not None:
                s   += nl.toString()

        for name in self.cardRef:   # Add cards
            c  = self.getObject(name, self.cards)
            if c is not None:
                s   += c.toString()

        if self.attach:             # Add attribute (e.g. for type='matdyn')
            s   += self.attach

        return s


    def structure(self):
        """
        Returns basic atomic structure information as list of tuples.
        Specific for 'pw' type

        Example: [('Ni', '52.98', 'Ni.pbe-nd-rrkjus.UPF'), (...)]
        """
        # Extract structure not from pw type input. Hard to do it from other types
        if self.type != "pw":
            return None

        list    = []        # list of structure
        card    = self.card("atomic_species")

        for l in card.lines():     # Should have format: "<Label> <Mass> <Pseudo-Potential>"
            l   = l.strip()
            if l == "":     # Empty line
                continue
            list.append(l.split())

        return list


    def _exists(self, name, list):
        "Checks if lowered name is in the list"
        if not name:
            return False
        
        name    = name.lower()
        if name in list:
            return True

        return False        


def _import(package):
    return __import__(package, globals(), locals(), [''], -1)




