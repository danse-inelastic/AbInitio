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

class Filter(object):

    def __init__(self, name = None):
        """
        Parameters:
            name:       str
                Name of the filter
        """
        self._name  = name


    def addNamelist(self, namelist, param):
    def addCard

    def apply(self, input):
        """
        Applies filter to the input

        Parameters:
            input:      object (QEInput)
                Input object 
        """
        pass

__date__ = "$Jul 28, 2010 2:35:31 PM$"


