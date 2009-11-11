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


class SimParams:

    def __init__(self, director):
        self._director  = director
#        self._simtype   = type
#        self._simlist   = self._getSimlist(type)


    def getLink(self, id):
        #inputs  = self._director.clerk.getConfigurations(where="simulationId='%s' AND type='settings'" % id)

        link = Link(label="settings.conf", Class="action-link", onclick=load(actor="espresso/server-view", routine="link", id=id))
#        link    = Link(label=self._inputText(orderedInputs[i]), Class="action-link",
#                             onclick=load(actor=self._getActor(orderedInputs[i]),
#                                          routine="link", id=id, type=self._simlist[i],
#                                          configid=self._getId(self._simlist[i], inputs)))
#        orderedInputs   = self._orderInput(self._simlist, inputs)
#
#        splitter    = Splitter(orientation='horizontal')
#        listsize    = len(self._simlist)
#
#        for i in range(listsize):
#            section     = splitter.section()
#            section.add(Paragraph(text=self._simlist[i]))   # Simulation type
#            section.add(Link(label=self._inputText(orderedInputs[i]), Class="action-link",
#                             onclick=load(actor=self._getActor(orderedInputs[i]),
#                                          routine="link", id=id, type=self._simlist[i],
#                                          configid=self._getId(self._simlist[i], inputs)))   # Passes config type (not id)
#                        )
#
#            if i != listsize - 1:   # No arrow for last config
#                sep     = splitter.section()        # Separator
#                sep.add(Paragraph(text=" ----> "))

        return link


    def _inputText(self, input):
        """Returns"""
        if input:
            return input.filename

        return "Add"


    def _getActor(self, input):
        """Returns proper actor depending if 'input' exists"""
        if input:   # View
            return "espresso/input-view"

        return "espresso/input-add" # Create New


#    def _orderInput(self, simlist, inputs):
#        """Orders input according to simlist (E.g. simlist = ("PW", "PH") )"""
#        newinputs   = []
#
#        for name in simlist:
#            newinputs.append(self._configObject(name, inputs))
#
#        return newinputs
#
#
#    def _configObject(self, type, inputs):
#        """Returns object if simulation type exists or None otherwise"""
#        for sim in inputs:
#            if sim.type == type:
#                return sim
#
#        return None


#    def _getSimlist(self, type):
#        if type in SIMCHAINS:
#            return SIMCHAINS[type]
#
#        return ()


    def _getId(self, type, inputs):
        config  = self._configObject(type, inputs)
        if config:
            return config.id

        return ""

if __name__ == "__main__":
    chain   = SimParams(None)


__date__ = "$Nov 10, 2009 5:52:02 PM$"


