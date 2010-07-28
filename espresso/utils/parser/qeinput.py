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

#Type of the configuration file can be:
#type =
#    'pw'               - (default)
#    'ph'               -
#    'pp'               -
#    'bands'            -
#    'cp'               - 
#    'cppp'             -
#    'd3'               -
#    'dos'              -
#    'dynmat'           -
#    'initial_state'    -
#    'gipaw'            -
#    'd1'               -
#    'matdyn'           -
#    'projwfc'          -
#    'pwcond'           -
#    'q2r'              -

class QEInput(object):

    # Either filename or config (not both) can be specified
    def __init__(self, filename=None, config=None, type='pw'):
        self.header     = None
        self.filename   = filename
        self.config     = config
        self.parser     = QEParser(filename, config, type)
        self.type       = type
        self.namelists  = OrderedDict()
        self.cards      = OrderedDict()
        self.attach     = None          # Specific for 'matdyn', 'dynmat'
        self.namelistRef    = None
        self.cardRef        = None
        self.qe         = [self.header, self.namelists, self.cards, self.attach]

    def parse(self):
        """ Parses the configuration file and stores the values in qe dictionary """
        (self.header, self.namelists, self.cards, self.attach) = self.parser.parse()


    def createNamelist(self, name):
        """Creates namelist and adds to QEInput. """
        nl  = Namelist(name)
        self.namelists[name] = nl
        return nl


    def addNamelist(self, namelist):
        """Adds namelist. """
        self.namelists[namelist.name()] = namelist


    def removeNamelist(self, name):
        try:
            del(self.namelists[name])
        except KeyError:    # parameter is not present
            return


    def namelist(self, name):
        "Returns namelist specified by name if exists or create a new one"
        if self.namelistExists(name):   # If exists, return namelist
            return self.namelists[name]

        return self.createNamelist(name) # otherwise create a new one


    def namelistExists(self, name):
        "Checks if namelist specified by name exists"
        return self._exists(name, self.namelists.keys())


    def createCard(self, name):
        "Creates card and adds to QEInput. "
        card    = Card(name)
        self.cards[name] = card
        return card


    def addCard(self, card):
        "Adds card"
        self.cards[card.name()] = card


    def removeCard(self, name):
        try:
            del(self.cards[name])
        except KeyError:    # parameter is not present
            return


    def attach(self):
        "Returns attachment"
        return self.attach


    def addAttach(self, text):
        """
        Sets attachment to some string.
        If attachment is not None it still will be overwritten
        """
        self.attach = text


    def removeAttach(self):
        "Sets attachment to None"
        self.attach = None


    def card(self, name):
        "Returns card specified by name if exists or create a new one"
        if self.cardExists(name):        # If exists, return card
            return self.cards[name]

        return self.createCard(name)


    def cardExists(self, name):
        "Checks if card specified by name exists"
        return self._exists(name, self.cards.keys())


    def getObject(self, name, dict):
        """Returns object that corresponds to 'name'"""
        for n in dict.values():
            if n.name() == name:
                return dict[name]

        return None


    def save(self, filename=None):
        """ Saves the QEInput to the configuration file"""
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
        return self.type


    def structure(self):
        """Returns basic structure information as list tuples
        Example: [('Ni', '52.98', 'Ni.pbe-nd-rrkjus.UPF'), (...)]
        """
        # Hard to extract structure not from pw type input
        # TODO: Should also have "atomic_species" card
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


    def toString(self):
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


    def _exists(self, name, list):
        if name in list:
            return True

        return False
        


def _import(package):
    return __import__(package, globals(), locals(), [''], -1)




