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
Card - class that corresponds to Card in QE config file
"""

# TODO: Add some more convenience methods?

from block import Block

class Card(Block):

    def __init__(self, name, arg = None):
        self._arg      = arg
        self._name     = name.lower() # keeps lower name
        self._lines    = []


#    def name(self):
#        "Returns name of the card"
#        return self._name
#
#
#    def setName(self, name):
#        """
#        Sets name
#        """
#        self._name = name.lower()


    def arg(self):
        return self._arg


    def setArg(self, arg):
        if arg is None:
            self._arg   = None
            return

        self._arg    = arg.lower()


    def line(self, num):
        """Returns value of parameter 'param'"""
        self._checkRange(num)
        return self._lines[num]


    def lines(self):
        "Returns list of lines"
        return self._lines


    def addLine(self, line):
        self._lines.append(line)


    def editLines(self, lines):
        """Replaces lines by new 'lines' (list) """
        self._lines    = lines


    def removeLine(self, num):
        self._checkRange(num)
        self.lines.pop(num)


    def removeLines(self):
        self._lines = []


    def __checkRange(self, num):
        assert num > 0
        assert len(self._lines) > num


    def toString(self, indent=1, br="\n"):
        """
        Dumps card as a sting

        Parameters:
            indent:     int
                Number of spaces in indent for parameters
            br:         str
                Separator between parameters
        """
        ind  = ""
        for i in range(indent):    # Populate indent
            ind += " "

        s = self.name().upper()
        if self._arg is not None:
            s += ' (%s)%s' % (self._arg, br)
        else:
            s += br

        for l in self._lines:
            s += '%s%s%s' % (ind, l.strip().strip('\n').strip(br), br)

        return s

__date__ = "$Aug 27, 2009 7:34:32 AM$"




