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

import configPW
from ovini.utils.orderedDict import OrderedDict

class Card():
    """Card class that corresponds to Card in QE config file"""
    # May be add some more convenience methods?

    def __init__(self, name):
        # Verifies if the card is valid
        try:
            if name not in configPW.cards:
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

__date__ = "$Aug 27, 2009 7:34:32 AM$"


