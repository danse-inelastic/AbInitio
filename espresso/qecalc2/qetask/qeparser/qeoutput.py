#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# QEcalc              by DANSE Inelastic group
#                     Brent Fultz
#                     California Institute of Technology
#                     (C) 2009  All Rights Reserved
#
# File coded by:      Nikolay Markovskiy
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class QEOutput(object):
    def __init__(self, filename, type='pw'):
        self.filename   = filename
        self.type       = type
        self.parsers    = None
        outModule = __import__("outputs." + self.type, globals(), \
                                locals(), ['Output'], -1)
        self.output = outModule.Output()

    def listParsers(self):
        parserList = []
        for parserName in self.output.parsers:
            parserList.append(parserName)
        return parserList

    def parse(self, parserList = 'all', fname = None):
        if parserList == 'all':
            parserList = self.listParsers()
        if fname == None:
            filename = self.filename
        else:
            filename = fname
        properties = {}
        for parserName in parserList:
            properties[parserName] = self.output.parse(parserName, filename)
        return properties


def test():
    qeOut = QEOutput('scf.out', 'pw')
    print qeOut.listParsers()
    properties = qeOut.parse()
    print properties
    properties = qeOut.parse(['total energy', 'stress'])
    print properties

if __name__ == "__main__":
    test()
    print "Hello World";

__author__="Nikolay Markovskiy"
__date__ ="$Oct 18, 2009 7:51:00 PM$"