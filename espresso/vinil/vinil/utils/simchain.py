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

from vinil.utils.const import SIMULATIONS

from luban.content.Splitter import Splitter
from luban.content.Paragraph import Paragraph
from luban.content import load
from luban.content.Link import Link


class SimChain:

    def __init__(self, director, type=None):
        self._director  = director
        self._simtype   = type
        self._simlist   = self._getSimlist(type)


    def chain(self, id):
        inputs          = self._director.clerk.getConfigurations(where="simulationId='%s'" % id)
        orderedInputs   = self._orderInput(self._simlist, inputs)

        splitter    = Splitter(orientation='horizontal')
        listsize    = len(self._simlist)

        for i in range(listsize):
            section     = splitter.section()
            section.add(Paragraph(text=self._simlist[i]))
            section.add(Link(label=self._inputText(orderedInputs[i]), Class="action-link",
                             onclick=load(actor="espresso-input", routine="link",
                             id=id, type=self._simlist[i]))
                        )

            if i != listsize - 1:   # No arrow for last config
                sep     = splitter.section()        # Separator
                sep.add(Paragraph(text=" ----> "))

        return splitter


    def _inputText(self, input):
        if input:
            return input.filename

        return "Add"


    def _orderInput(self, simlist, inputs):
        """E.g. simlist = ("PW", "PH")"""
        newinputs   = []

        for name in simlist:
            newinputs.append(self._simObject(name, inputs))

        return newinputs


    def _simObject(self, name, inputs):
        """Returns object if simulation type exists or None otherwise"""
        for sim in inputs:
            if sim.type == name:
                return sim

        return None


    def _getSimlist(self, type):
        if type in SIMULATIONS:
            return SIMULATIONS[type]

        return ()


if __name__ == "__main__":
    chain   = SimChain(None, None)


__date__ = "$Nov 9, 2009 10:50:54 AM$"


# *************** DEAD CODE ********************
#        one     = splitter.section()
#        one.add(Paragraph(text="PW "))
#        one.add(Link(label=filename, Class="action-link", onclick=load(actor="espresso-set-config", routine="link", id=id)))
#        sep     = splitter.section()        # Separator
#        sep.add(Paragraph(text=" -> "))
#        two     = splitter.section()
#        two.add(Paragraph(text="DOS"))
#        two.add(Link(label="Add"))


