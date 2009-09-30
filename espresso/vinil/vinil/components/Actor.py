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

        page.skeleton.path.add(self.path(self.pathlist))
        page.maindoc.add(self.content(director))

        return page

    def link(self, director, document=None):
        return self.getActions(self.content(director), self.pathlist)

    def content(self, director):
        document = Document(title='Not Implemented')

        text = "Sorry, the page is not implemented!"
        p = Paragraph(text=text)
        document.add(p)

        return document

    # set self.list before using path!
    def path(self, pathlist):
        from vinil.utils.pathbuilder import PathBuilder
        pb = PathBuilder(pathlist)
        return pb.buildPath()

    def setPathList(self, pathlist):
        self.pathlist = pathlist

    def getActions(self, document, pathlist):
        actions = []
        actions.append(select(id='path-content').replaceContent(self.path(pathlist)))
        actions.append(select(id='maindoc').replaceContent(document))

        return actions

    def __init__(self, *args, **kwds):
        super(Actor, self).__init__(*args, **kwds)
        self.list = [("Home","greet","link"),
                     ("Not Implemented", None, None)]
        return

__date__ = "$Sep 28, 2009 4:20:04 PM$"
