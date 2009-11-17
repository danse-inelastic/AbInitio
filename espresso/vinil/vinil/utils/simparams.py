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
        self._type      = "settings"

    def getLink(self, id):      # simulation
        settings  = self._director.clerk.getConfigurations(where="simulationId='%s' AND type='%s'" % (id, self._type))

        link = Link(label=self._label(settings), Class="action-link",
                    onclick=load(actor      = self._getActor(settings),
                                 routine    = "link",
                                 id         = id,
                                 configid   = self._configid(settings))
                    )

        return link


    def _label(self, settings):
        """Returns filename"""
        if settings:
            return settings[0].filename

        return "Add"


    def _configid(self, settings):
        """Returns config id"""
        if settings:
            return settings[0].id

        return None


    def _getActor(self, settings):
        """Returns proper actor depending if 'input' exists"""
        if settings:   # View
            return "espresso/settings-view"

        return "espresso/settings-add" # Create New


if __name__ == "__main__":
    chain   = SimParams(None)


__date__ = "$Nov 10, 2009 5:52:02 PM$"



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
#        self._simtype   = type
#        self._simlist   = self._getSimlist(type)


#    def _getSimlist(self, type):
#        if type in SIMCHAINS:
#            return SIMCHAINS[type]
#
#        return ()

#    def _getId(self, type, inputs):
#        config  = self._configObject(type, inputs)
#        if config:
#            return config.id
#
#        return ""



