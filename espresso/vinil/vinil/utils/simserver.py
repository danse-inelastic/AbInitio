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

from luban.content.Splitter import Splitter
from luban.content.Paragraph import Paragraph
from luban.content import load
from luban.content.Link import Link
import ConfigParser

class SimServer:

    def __init__(self, director):
        self._director  = director
        self._clerk     = director.clerk
        self._type      = "settings"

    def getServer(self, id):      # simulation id
        settings  = self._clerk.getConfigurations(where="simulationId='%s' AND type='%s'" % (id, self._type))

        if self._serverIsSet(settings):
            text    = Link(label=self._label(settings), Class="action-link",
                        onclick=load(actor="espresso/settings-view",
                        routine="link", id=id))
        else:
            text    = Paragraph(text="None")

        return text

    def _servreIsSet(settings):
        """Checks if server is set"""
        return False
        #if not settings



    def _label(self, settings):
        """Returns filename"""
        if settings:
            return settings[0].filename

        return "Add"


    def _getActor(self, settings):
        """Returns proper actor depending if 'input' exists"""
        if settings:   # View
            return "espresso/settings-view"

        return "espresso/settings-add" # Create New


if __name__ == "__main__":
    chain   = SimParams(None)


__date__ = "$Nov 11, 2009 1:21:52 PM$"


