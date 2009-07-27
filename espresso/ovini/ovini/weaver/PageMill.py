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

from opal.weaver.PageMill import PageMill as base

class PageMill(base):

    def onPage(self, page):

        head = page._head
        body = page._body

        #header = body.pageHeader()

        """
        from opal.content.Literal import Literal
        l = Literal()
        l.text = ['<h2><a href="/">Ovini</a> (Opal VNF Mini)</h2>']
        header.add(l)
        
        footer = body.pageFooter()
        """


        #print head, body, page

        """
        # render javascripts
        includes, scripts = self.javascriptWeaver.render(page)

        # add includes to head
        for inc in includes:
            head.script(src=inc)

        # add scripts to head
        for script in scripts:
            s = head.script()
            s.script = script
        """

        return base.onPage(self, page)

    """
    def onBody(self, body):
        from opal.content.Literal import Literal
        l = Literal()
        l.text = ['<h2><a href="/">Ovini</a> (Opal VNF Mini)</h2>']
        body.add(l)

        return base.onBody(self, body)
    """

    def __init__(self, configurations):
        base.__init__(self)

        from BodyMill import BodyMill
        self.bodyMill = BodyMill(self.tagger)
        #print self.bodyMill


        """
        from JavaScriptWeaver import JavaScriptWeaver
        self.javascriptWeaver = JavaScriptWeaver(configurations)

        tagger = self.bodyMill.tagger
        from StructuralMill import StructuralMill
        self.bodyMill.structuralMill = StructuralMill(tagger, configurations)
        self.bodyMill.structuralMill.master = self
        """

        return

# version
__id__ = "$Id$"
__date__ = "$Jul 20, 2009 10:20:19 AM$"


