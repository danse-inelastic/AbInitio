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

class Namelist:
    
    def __init__(self, name):
        """
        Parameters:
            name:       str
                Name of the namelist in lower case. Example: "control"
        """
        self.__name = name.lower()  # keeps lower name
        self.params = OrderedDict() # Replace dictionary by ordered dictionry


    def name(self):
        "Return name of the namelist"
        return self.__name


    def setName(self, name):
        "Set name in lower case"
        self.__name = name.lower()


    def get(self, param, quotes = True):
        """
        Returns paramater value. If not exists, return None

        Parameters:
            param:      str
                Parameter of the namelist
            quotes:     bool
                True - if add quotes '' to parameters value, False - otherwise
                
        Note: replaces param()
        """
        if not self.__paramExists(param):
            return None
        
        if quotes:
            return self.params[param]
        
        return self._unquote(self.params[param])

#        return self.params[param]



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
#        if not self.__paramExists(param):

        param   = param.lower()
        if quotes:
            val     = self._quote(val)

        self.params[param] = val


    def remove(self, param):
        """Deletes parameter"""
        if self.__paramExists(param):
            del(self.params[param])


    def exists(self,param):
        return self.__paramExists(param)


    def _quote(self, val):
        return "'" + val.strip('"').strip("'") + "'"
    
    
    def _unquote(self, val):
        return val.strip('"').strip("'")


    def toString(self, indent="    ", br="\n"):
        # Dump namelist
        # Should I use space?
        s = '&%s%s' % (self.name().upper(), br)

        for p in self.params.keys():
            s += '%s%s = %s,%s' % (indent, p, self.params[p], br)

        s += "/%s" % br 
        return s


    def __paramExists(self, param):
        try:
            param = param.lower()
            self.params[param]
            return True
        except KeyError:    # parameter is not present
            return False


    # Depricated methods:
    # Depricated
    def param(self, param, quotes = True):
        """Returns value of parameter 'param'"""
        return self.get(param, quotes)
    
#        if self.__paramExists(param):
#            if quotes:
#                return self.params[param]
#            else:
#                return self._unquote(self.params[param])
#            return self.params[param]
#
#        return None

    # Depricated: Use set() instead!
    def add(self, param, val, quotes = False):
        "Adds parameter to the namelist"
        self.set(param, val, quotes)
        
        #param = param.lower()   # Should be lowered?
#        if quotes:
#            val     = self._quote(val)
#
#        self.params[param]  = val

__date__ = "$Aug 27, 2009 7:30:39 AM$"





