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
from vinil.utils.orderedDict import OrderedDict

class Namelist():
    """Namelist class that corresponds to Namelist in QE config file"""

    def __init__(self, name):
        # Verifies if the namelist is valid
        try:
            if name not in configPW.namelists:
                raise NameError('Not valid namelist')
        except NameError:
            raise

        self.__name = name.lower() # keeps lower name
        self.params = OrderedDict() # Replace dictionary by ordered dictionry

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
        # Should I use space?
        s = '&%s%s' % (self.name().upper(), br)

        for p in self.params.keys():
            s += '%s%s = %s,%s' % (indent, p, self.params[p], br)

        s += "/%s" % br
        return s

__date__ = "$Aug 27, 2009 7:30:39 AM$"


