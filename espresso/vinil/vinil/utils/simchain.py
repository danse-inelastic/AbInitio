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


class SimChain():

    def __init__(self, director):
        self._director  = director

    def chain(self, id):
        config  = self._director.clerk.getConfigurations(where="simulationId='%s'" % id)
        if config:
            filename    = config[0].filename
        else:
            filename    = "Add"

        splitter    = Splitter(orientation='horizontal')
        one     = splitter.section()
        one.add(Paragraph(text="PW "))
        one.add(Link(label=filename, Class="action-link", onclick=load(actor="espresso-set-config", routine="link", id=id)))
        sep     = splitter.section()
        sep.add(Paragraph(text=" ----> "))
        two     = splitter.section()
        two.add(Paragraph(text="DOS"))
        two.add(Link(label="Add"))

        return splitter
    
__date__ = "$Nov 9, 2009 10:50:54 AM$"


