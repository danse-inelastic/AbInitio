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

from luban.content.Splitter import SplitSection
from luban.content.Paragraph import Paragraph
from luban.content.Document import Document
from luban.content.HtmlDocument import HtmlDocument
from luban.content import load


# Excepts double list of the form:
# ((actorBack, routineBack),(actorNext, routineNext))

class Stepper:
    def __init__(self, linklist):
        self.checkRange(linklist)
        self.linklist = linklist

    def checkRange(self, linklist):
        size = len(linklist)
        if size != 2:
            raise IndexError

    def next(self, link):
        if len(link) != 0 and link[0] is not None:

            return HtmlDocument(text='<img src="images/icons/next.png"/>',
                            Class = "stepper-item",
                            onclick=load(actor=link[0], routine=self.routine(link[1]), id=link[2]))

        return Paragraph(Class="stepper-space-next")

    def back(self, link):
        if len(link) != 0 and link[0] is not None:

            return HtmlDocument(text='<img src="images/icons/back.png"/>',
                            Class = "stepper-item",
                            onclick=load(actor=link[0], routine=self.routine(link[1]), id=link[2]))

        return Paragraph(Class="stepper-space-back")


    def routine(self, name):
        if name is not None:
            return name

        return 'default'

    def getStepper(self):
        s_stepper   = SplitSection()
        d = Document(Class="stepper")
        d.add(self.next(self.linklist[1]))
        d.add(Paragraph(Class="stepper-space"))
        d.add(self.back(self.linklist[0]))

        s_stepper.add(d)

        return s_stepper

    def toString(self):
        s = ""
        if len(self.linklist[0]) != 0:
            s += "(%s, %s, %s) <== " % (self.linklist[0][0], self.linklist[0][1], self.linklist[0][2])

        if len(self.linklist[1]) != 0:
            s += " ==> (%s, %s, %s)" % (self.linklist[1][0], self.linklist[1][1], self.linklist[1][2])

        print s

def testStepper():
    linklist = (("espresso", "link", "5"), ("espresso-material", "link", "6"))
    stepper     = Stepper(linklist)
    stepper.toString()


if __name__ == "__main__":
    testStepper()

__date__ = "$Sep 30, 2009 3:15:16 PM$"


