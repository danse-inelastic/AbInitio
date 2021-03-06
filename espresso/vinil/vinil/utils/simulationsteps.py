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

from vinil.utils.const import STEPS

# SimulationSteps display the steps user should follow to run simulation
# in Quantum Espresso (make general?)
# Steps are counted from 1 (current=1 is the first step)

class SimulationSteps:
    def __init__(self, stepslist=None, current=1):
        self.setStepsList(stepslist)
        self.checkRange(current, self.stepslist)
        self.current = current

    def setStepsList(self, stepslist):        
        if stepslist is None:
            self.stepslist = STEPS
        else:
            self.stepslist  = stepslist

    def getSteps(self):
        s_steps     = SplitSection()

        for i in range(len(self.stepslist)):
            if i == 0:
                itemclass   = "step-item-first"
            else:
                itemclass   = "step-item"
                
            d = Document(Class=itemclass)
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

    def addStepItem(self):
        pass


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



def testSimulationSteps():
    simsteps = SimulationSteps(current=3)
    simsteps.toString()


if __name__ == "__main__":
    testSimulationSteps()
    

__date__ = "$Sep 30, 2009 11:07:27 AM$"


