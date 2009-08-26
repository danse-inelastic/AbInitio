#!/usr/bin/env python

# Parses configuration file (for pw.x) and stores it in dictionary.
# It also allows to dump the dictionary to create the configuration file
# The main use case is to parse the existing configuration file, change 
# some values and save it back to the configuration file.

# Input parameters are defined in INPUT_PW.html 

"""
Stability issues:
- Parsing goes line by line. 
- Namelist starts with '&' and ends with '/' on a separate line
- Card starts with card title on a separate line and values between card titles.
- Prints both Namelists and Cards in capital 
- Would like to use ordered dictionary?
- Refactoring?  Introduce class relation: Namelist(Block), Card(Block)
"""

import configPW

ref = {'control': configPW.namelist_control}

namelistsPW = ('control',
             'system',
             'electrons',
             'ions',
             'cell',
             'phonon')

cardsPW = ('atomic_species', 
           'atomic_positions',
           'k_points automatic',
           'cell_parameters',
           'climbing_images',
           'constraints',
           'collective_vars',
           'occupations')

import sys

class Namelist():
    """Namelist class that corresponds to Namelist in QE config file"""

    def __init__(self, name):
        # Verifies if the namelist is valid
        try:
            if name not in namelistsPW:
                raise NameError('Not valid namelist')
        except NameError:
            raise

        self.__name = name.lower() # keeps lower name
        self.params = {}

    def name(self):
        return self.__name

    def setName(self, name):
        self.__name = name.lower()

    def param(self, param):
        """Returns value of parameter 'param'"""
        if self.__paramExists(param):
            return self.params[param]

    def addParam(self, param, val):
        # Add verification?
        param = param.lower()
        self.params[param]  = val

    def editParam(self, param, val):
        """Edits parameter. If it doesn't exist, it just ignores it """
        if self.__paramExists(param):
            self.params[param] = val

    def removeParam(self, param):
        """Deletes parameter"""
        if self.__paramExists(param):
            del(self.params[param])

    def __paramExists(self, param):
        try:
            param = param.lower()
            self.params[param]
            return True
        except KeyError:    # parameter is not present
            return False 

    def toString(self, indent="    ", br="\n"):
        # Dump namelist
        s = '&%s%s' % (self.name().upper(), br)

        for p in self.params.keys():
            s += '%s%s = %s%s' % (indent, p, self.params[p], br)

        s += "/%s" % br
        return s

class Card():
    """Card class that corresponds to Card in QE config file"""
    # May be add some more convenience methods?

    def __init__(self, name):
        # Verifies if the card is valid
        try:
            if name not in cardsPW:
                raise NameError('Not valid card')
        except NameError:
            raise

        self.__name     = name.lower() # keeps lower name
        self.__lines    = []

    def name(self):
        return self.__name

    def setName(self, name):
        self.__name = name.lower()

    def line(self, num):
        """Returns value of parameter 'param'"""
        self.__checkRange(num)
        return self.__lines[num]

    def addLine(self, line):
        self.__lines.append(line)

    def editLines(self, lines):
        """Replaces lines by new 'lines' (list) """
        self.__lines    = lines

    def removeLine(self, num):
        self.__checkRange(num)
        self.lines.pop(num)

    def removeLines(self):
        self.__lines = []

    def __checkRange(self, num):
        assert num > 0
        assert len(self.__lines) > num

    def toString(self, indent=" ", br="\n"):
        # Dump card
        s = '%s%s' % (self.name().upper(), br)

        for l in self.__lines:
            s += '%s%s%s' % (indent, l, br)

        return s

class QEConfig(object):
    """Quantum Espresso configuration class. It can:
    - Parse existing configuration file
    - Add, Edit or Remove parameters from/to namelist or card
    """

    # Either filename or configstr (not both) can be specified
    def __init__(self, filename=None, configstr=None):
        self.filename   = filename
        self.configstr  = configstr
        self.namelists  = {}
        self.cards      = {}
        self.qe         = [self.namelists, self.cards]

    def createNamelist(self, name):
        """Creates namelist and adds to QEConfig. """
        nl  = Namelist(name)
        self.namelists[name] = nl

    def addNamelist(self, namelist):
        """Adds namelist. """
        self.namelists[namelist.name()] = namelist

    def removeNamelist(self, name):
        try:
            del(self.namelists[name])
        except KeyError:    # parameter is not present
            return

    def namelist(self, name):
        # Do I need editNamelist()?
        try:
            return self.namelists[name]
        except KeyError:    # parameter is not present
            raise

    def createCard(self, name):
        """Creates card and adds to QEConfig. """
        self.cards[name] = Card(name)

    def addCard(self, card):
        """Adds card. """
        self.cards[card.name()] = card

    def removeCard(self, name):
        try:
            del(self.cards[name])
        except KeyError:    # parameter is not present
            return

    def card(self, name):
        # Do I need editNamelist()?
        try:
            return self.cards[name]
        except KeyError:    # parameter is not present
            raise

    def toString(self):
        s = ''
        for nl in self.namelists.values():
            s += nl.toString()

        for c in self.cards.values():
            s += c.toString()
        return s

    def save(self, filename=None):
        """ Saves the QEConfig to the configuration file"""
        default = "config.out"

        if filename is None:
            if self.filename is not None:
                filename = self.filename
            else:
                filename = default

        f = open(filename, "w")
        f.write(self.toString())
        f.close()

    def parse(self):
        """ Parses the configuration file and stores the values in qe dictionary """

        if self.filename is not None:
            try:
                f = open(self.filename)
            except IOError:
                print "I/O error"
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise

            lines       = f.readlines()         # Returns list of lines.
            f.close()

        elif self.configstr is not None:
            lines       = self.configstr.splitlines()

        else:
            print "Error: You haven't specify any config file"

        lines       = self.__clearLines(lines)
        marks       = self.__getMarks(lines)

        for m in marks:
            block   = self.__addBlock(m[0], lines[m[1]:m[2]])


    def __clearLines(self, lines):
        """ Strips lines from white spaces, commas and empty lines"""

        cl = []     # Lines without white spaces and empty lines
        for l in lines:
            l = l.strip().strip(',') # Remove both lead and trailing whitespace, including '\n' and comma
            if l == '':
                continue

            cl.append(l)

        return cl

    def __addBlock(self, type, slice):
        """ Adds block (namelist or card to ) """

        # Improve calls?
        if type == 'namelist':
            self.__addNamelist(slice)
        elif type == 'card':
            self.__addCard(slice)

        return

    def __addNamelist(self, slice):
        """Adds namelist based on slice """
        name    = slice[0].strip('&')
        nl      = Namelist(name)

        for s in slice[1:]:
            p   = self.getParam(s)
            nl.addParam(p[0], p[1])

        self.namelists[name] = nl

    def __addCard(self, slice):
        """Adds card"""
        name    = slice[0].lower()
        c = Card(name)

        for s in slice[1:]:
            c.addLine(s)

        self.cards[name]    = c

    def __getMarks(self, lines):
        # TODO: Cumbersome method, rewrite it
        """
        Determines start and end of namelist and card blocks: [type, start, end]
        E.g ['namelist', 0, 7] for CONTROL namelist
        Iterate over number of lines. Empty lines are included
        Not tested very well yet
        """
        blocklist   = []
        isNamelist  = False
        isCard      = False
        size        = len(lines)

        for i in range(size):
            l = lines[i]
            # We suppose that namelists and card do not intersect
            # Namelist part

            # Namelist end
            if l[0] == '/' and isNamelist:
                isNamelist  = False
                block.append(i)
                blocklist.append(block)

            # Namelist start
            if l[0] == '&' and not isNamelist:
                name = l[1:].lower()

                if not name in namelistsPW:
                    continue             # namelist is not recognizable

                block       = []
                isNamelist  = True
                block.append('namelist')
                block.append(i)

            # Card part
            line    = l.lower()
            # Card end
            if line in cardsPW and isCard:
                #print "End: %s, line: %d" % (line, i-1)
                isCard  = False
                block.append(i)
                blocklist.append(block)

            if i == size-1 and isCard:
                isCard  = False
                block.append(i+1)
                blocklist.append(block)

            # Card start
            if line in cardsPW and not isCard:
                #print "Start: %s, line: %d" % (line, i)
                block   = []
                isCard  = True
                block.append('card')
                block.append(i)

        return blocklist

        # Example return: [['namelist', 0, 7], ['namelist', 8, 20]]

    def getParam(self, s):
        """ Takes string like 'a = 2' and returns tuple ('a', 2) """

        ss = s.split('=')
        for i in range(len(ss)):
            ss[i] = ss[i].strip()

        val = ss[1]

        # Assume that there are two values only: (variable, value) pair
        assert len(ss) == 2

        return (ss[0], val)


def testCreateConfig():
    print "Testing creation of config file"
    qe  = QEConfig()
    nl  = Namelist('control')
    nl.addParam('title', "'Ni'")
    nl.addParam('restart_mode', "'from_scratch'")
    print "Adding parameters to namelist:\n%s" % nl.toString()
    nl.editParam('title', "'Fe'")
    qe.addNamelist(nl)
    print "Adding namelist to QEConfig:\n%s" % qe.toString()

    c = Card('atomic_species')
    c.addLine('Ni  26.98  Ni.pbe-nd-rrkjus.UPF')
    print "Adding line to card:\n%s" % c.toString()
    qe.addCard(c)
    print "Adding card to QEConsig:\n%s" % qe.toString()
    qe.save()

def testParseConfig():
    print "Testing parsing config file"
    qe  = QEConfig("vini/ni.scf.in")
    qe.parse()
    print qe.toString()
    nl  = qe.namelist('control')
    nl.addParam('title', 'Ni')
    nl.removeParam('restart_mode')
    qe.removeCard('atomic_species')
    nl.editParam('calculation', "'nscf'")
    c = qe.card('atomic_positions')
    c.editLines(['Say Hi! :)'])
    print qe.toString()
    qe.save("ni.scf.in.mod")

if __name__ == "__main__":
    testCreateConfig()
    testParseConfig()
