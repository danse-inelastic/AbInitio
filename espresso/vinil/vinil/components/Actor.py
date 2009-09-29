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

from luban.content import select
from luban.content.Paragraph import Paragraph
from luban.content.Document import Document

from opal.components.Actor import Actor as base
class Actor(base):
    
    def default(self, director):
        page = director.retrieveVisual('template')

        page.skeleton.path.add(self.path())
        page.maindoc.add(self.content(director))

        return page

    def link(self, director, document=None):
        actions = []
        actions.append(select(id='path-content').replaceContent(self.path()))
        actions.append(select(id='maindoc').replaceContent(self.content(director)))

        return actions

    def content(self, director):
        document = Document(title='Not Implemented')

        text = "Sorry, the page is not implemented!"
        p = Paragraph(text=text)
        document.add(p)

        return document

    # set self.list before using path!
    def path(self):
        from vinil.utils.pathbuilder import PathBuilder
        pb = PathBuilder(self.list)
        return pb.buildPath()

    def setPathList(self, list):
        self.list = list
    
    def __init__(self, *args, **kwds):
        super(Actor, self).__init__(*args, **kwds)
        self.list = [("Home","greet","link"),
                     ("Not Implemented", None, None)]
        return

__date__ = "$Sep 28, 2009 4:20:04 PM$"
