#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                      Nikolay Markovskiy, Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

class Calc: pass

class Task: pass

class VASPCalc(Calc): pass

class QECalc(Calc):
    tasks
    settings

    def run(self): pass

    def submit(self, tasks): pass

class QEPhon(QECalc):
    self.tasks = [PWTask(), PHTask(), Q2RTask(), MATDYNTask()]

    def getPWTask(self):
        return self.tasks[0]

    def getOutput(self): pass

    def run(self):
        for t in self.tasks:
            t.launch()

class QETask(Task):
    self.inputParser
    self.outputParse
    self.qeconfig
    self.parameters

    def __init__(self): pass

    def parseInput(self, filename):
        self.qeconfig = self._parse(self.inputParser)
    
    def parseOutput(self, params):
        self.parameters = self._parse(self.ouputParser)

    def _parse(self, parser): pass
    
    def energy(self):
        self.parameters['energy']

    def launch(self):
        self.parseInput("filename")
        self.run()
        self.parseOutput(params)

class PWTask(QETask):
    self.qeconfig   = QEConfig(type="pw")
    def __init__(self, input):
        self.energy;



class PHTask(QETask): pass

class Q2RTask(QETask): pass


# Use cases:
# - Find energy
def calculateEnergy():
    qe  = QEPhon()
    qe.run()
    e1  = qe.getPWTask().energy()
    return e1

def test():
    qe  = QEPhon()
    pwtask  = PWTask("filename2")
    qe.setTask(task=pwtask)
    qe.run()
    pw  = qe.getPWTask()
    e1  = pw.energy()
    pw.set("ecut", "27")
    qe.run()
    e2  = pw.energy()


    e3  = qe.getPHTask()

    

__date__ = "$Oct 16, 2009 6:09:58 PM$"


