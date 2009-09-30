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

defaultSteps = ["Set Material",
                "Set Simulation Parameters",
                "Review Configuration",
                "Add to Jobs Queue"]

# SimulationSteps display the steps user should follow to run simulation
# in Quantum Espresso (make general?)
# Steps are counted from 1 (current=1 is the first step)

class SimulationSteps:
    def __init__(self, stepslist=None, current=1):
        self.setStepsList(stepslist)
        self.checkRange(current, self.stepslist)
        self.current = current


    def addStepItem(self):
        pass

    def setStepsList(self, stepslist):        
        if stepslist is None:
            self.stepslist = defaultSteps
        else:
            self.stepslist  = stepslist


    def getSteps(self):
        s_steps     = SplitSection()

        for i in range(len(self.stepslist)):
            d = Document(Class="step-item")
            num     = i+1
            classes = self.setClasses(num, self.current)
            d.add(Paragraph(text="%s" % num, Class=classes[0] ))
            d.add(Paragraph(text=self.stepslist[i], Class=classes[1]))
            s_steps.add(d)

        return s_steps


    # Return tumple of two class names (number, text)
    def setClasses(self, index, current):
        if index != current:
            return ("step-number-disabled", "step-text-disabled")
        
        return ("step-number", "step-text")

    def checkRange(self, index, list):
        if not self.withinRange(index-1, list):
            raise IndexError

    def withinRange(self, index, list):
        if index < len(list) and index >= 0:
            return True
        
        return False

    def toString(self):
        s = ""
        for i in range(len(self.stepslist)):
            num = i+1
            if num != self.current:
                s += "%s. " % num
            else:
                s += "[%s]. " % num
            s += "%s  " % self.stepslist[i]

        print s


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
                            onclick=load(actor=link[0], routine=self.routine(link[1])))

        return Paragraph(Class="stepper-space-next")

    def back(self, link):
        if len(link) != 0 and link[0] is not None:

            return HtmlDocument(text='<img src="images/icons/back.png"/>',
                            Class = "stepper-item",
                            onclick=load(actor=link[0], routine=self.routine(link[1])))

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
            s += "(%s, %s) <== " % (self.linklist[0][0], self.linklist[0][1])

        if len(self.linklist[1]) != 0:
            s += " ==> (%s, %s)" % (self.linklist[1][0], self.linklist[1][1])

        print s


def testSimulationSteps():
    simsteps    = SimulationSteps(current=3)
    simsteps.toString()

def testStepper():
    linklist = (("espresso", "link"), ("espresso-material", "link"))
    stepper     = Stepper(linklist)
    stepper.toString()


if __name__ == "__main__":
    testSimulationSteps()
    testStepper()
    

__date__ = "$Sep 30, 2009 11:07:27 AM$"


