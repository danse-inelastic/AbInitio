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
Namelist class that corresponds to Namelist in QE configuration file
"""

from orderedDict import OrderedDict
from block import Block

class Namelist(Block):
    
    def __init__(self, name):
        """
        Parameters:
            name:       str
                Name of the namelist in lower case. Example: "control"
        """
        self._name = name.lower()  # keeps lower name
        self.params = OrderedDict() # Replace dictionary by ordered dictionry


    def get(self, param, quotes = True):
        """
        Returns paramater value. If no parameter exists, return None.
        When quotes=True quotes are not added to parameter's value.

        Parameters:
            param:      str
                Parameter of the namelist
            quotes:     bool
                True - if add quotes '' to parameters value, False - otherwise
                
        Note: replaces param()
        """
        if not self._paramExists(param):
            return None

        param   = param.lower()
        if quotes:
            return self.params[param]
        
        return self._unquote(self.params[param])


    def set(self, param, val, quotes = False):
        """
        Sets existing parameter to the specified value.
        If no parameter exists, create one
        
        Parameters:
            param:      str
                Parameter name
            val:        str
                Parameter value
            quotes:     bool
                Add quotes to the value or not
        """
        param   = param.lower()
        if quotes:
            val     = self._quote(val)

        self.params[param] = val


    def remove(self, param):
        """
        Deletes parameter

        Parameters:
            param:      str
                Name of the parameter
        """
        if self._paramExists(param):
            del(self.params[param])


    def exists(self,param):
        """
        Checks if parameter exists in the namelist

        Parameters:
            param:      str
                Name of the parameter
        """
        return self._paramExists(param)


    def _quote(self, val):
        """
        Quotes value with "'" quote mark

        Parameters:
            val:        str
                Value to be quoted
        """
        return "'" + val.strip('"').strip("'") + "'"
    
    
    def _unquote(self, val):
        """
        Removes quotes "'" (unquotes) on both sides of the string

        Parameters:
            val:        str
                Value to be unquoted
        """
        return val.strip('"').strip("'")


    def toString(self, indent = 4, br = "\n"):
        """
        Dumps namelist as a sting
        
        Parameters:
            indent:     int
                Number of spaces in indent for parameters
            br:         str
                Separator between parameters
        """
        ind  = ""
        for i in range(indent):    # Populate indent
            ind += " "

        s = '&%s%s' % (self.name().upper(), br)

        for p in self.params.keys():
            s += '%s%s = %s,%s' % (ind, p, self.params[p], br)

        s += "/%s" % br 
        return s


    def _paramExists(self, param):
        """
        Checks if parameter exists in self.params

        Parameters:
            param:      str
                Name of the parameter
        """
        try:
            param = param.lower()
            self.params[param]
            return True
        except KeyError:    # parameter is not present
            return False


    # Depricated methods:
    # Depricated: Use get() instead
    def param(self, param, quotes = True):
        """Returns value of parameter 'param'"""
        return self.get(param, quotes)
    

    # Depricated: Use set() instead!
    def add(self, param, val, quotes = False):
        "Adds parameter to the namelist"
        self.set(param, val, quotes)
        

__date__ = "$Aug 27, 2009 7:30:39 AM$"





