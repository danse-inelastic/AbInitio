# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2010  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
"""
Unit tests for the following classes (modules):
    Namelist, Card, QEParser, QEInput
"""

import unittest
import fixtures
from namelist import Namelist
from card import Card
from qeparser import QEParser
from qeinput import QEInput

class QEParserTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_namelist_name(self):
        nl  = Namelist("Control")
        self.assertEqual(nl.name(), "control")
        nl.setName("SyStem")
        self.assertEqual(nl.name(), "system")


    def test_namelist_set_get(self):
        nl  = Namelist("control")
        self.assertEqual(nl.get("title"), None)
        
        nl.set("title", "hello")
        self.assertEqual(nl.get("title"), "hello")
        self.assertEqual(nl.get("title", quotes=True), "hello") # should not add quotes
        
        nl.set("Title", "'hello'")
        self.assertEqual(nl.get("titlE", quotes=False), "hello")
        self.assertEqual(nl.get("title"), "'hello'")


    def test_namelist_remove_exists(self):
        nl  = Namelist("control")
        nl.set("title", "hello")
        self.assertEqual(nl.exists("title"), True)
        self.assertEqual(nl.exists("Title"), True)

        nl.remove("title")
        self.assertEqual(nl.exists("title"), False)


    def test_namelist_tostring(self):
        nl  = Namelist("control")
        nl.set("title", "hello")
        self.assertEqual(nl.toString(), fixtures.assertNL)
        self.assertEqual(nl.toString(indent=3), fixtures.assertNL_space_3)


if __name__ == '__main__':
    unittest.main()
    
# card
# qeparse
# qeinput


#
#def testMatdyn():
#    parser    = QEParser(configText = textMatdyn, type="matdyn")
#    parser.parse()
#    print parser.toString()
#
#
#def testDynmat():
#    parser    = QEParser(configText = textDynmat, type="dynmat")
#    parser.parse()
#    print parser.toString()
#
#
#def testFile():
#    parser    = QEParser(filename = "../tests/ni.scf.in")
#    parser.parse()
#    print parser.toString()
#
#def testCards():
#    parser    = QEParser(configText = textCards)
#    parser.parse()
#    print parser.toString()
#
#
#def testComma():
#    parser          = QEParser(configText = textComma, type="matdyn")
#    parser.parse()
#    print parser.toString()
#
#
#def testHeader():
#    parser          = QEParser(configText = textHeader, type="ph")
#    parser.parse()
#    print parser.toString()
#
#def testMgB2():
#    parser          = QEParser(filename = "../tests/ph.mgb2.in", type="ph")
#    parser.parse()
#    print parser.toString()
#
#
#if __name__ == "__main__":
#    #testMatdyn()
#    #testDynmat()
#    #testFile()
#    #testCards()
#    #testComma()
#    #testHeader()
#    testMgB2()
#
#
## !!! FINISH
#
## Tests
#def testCreateConfig():
#    print "Testing creation of config file"
#    qe  = QEInput()
#    nl  = Namelist('control')
#    nl.add('title', "'Ni'")
#    nl.add('restart_mode', "'from_scratch'")
#    print "Adding parameters to namelist:\n%s" % nl.toString()
#    nl.set('title', "'Fe'")
#    qe.addNamelist(nl)
#    print "Adding namelist to QEInput:\n%s" % qe.toString()
#
#    c = Card('atomic_species')
#    c.addLine('Ni  26.98  Ni.pbe-nd-rrkjus.UPF')
#    print "Adding line to card:\n%s" % c.toString()
#    qe.addCard(c)
#    print "Adding card to QEInput:\n%s" % qe.toString()
#    #qe.save()
#
#
#def testParseConfig():
#    print "Testing parsing config file"
#    qe  = QEInput("../tests/ni.scf.in")
#    qe.parse()
#    print qe.toString()
#    nl  = qe.namelist('control')
#    nl.add('title', 'Ni')
#    nl.remove('restart_mode')
#    qe.removeCard('atomic_species')
#    nl.set('calculation', "'nscf'")
#    c = qe.card('atomic_positions')
#    c.editLines(['Say Hi! :)'])
#    print qe.toString()
#    #qe.save("../tests/ni.scf.in.mod")
#
#def testAttach():
#    qe  = QEInput("../tests/si.ph.in", type="ph")
#    qe.parse()
#    qe.save("../tests/si.ph.in.mod")
#
#
#if __name__ == "__main__":
#    #testCreateConfig()
#    #testParseConfig()
#    testAttach()

__date__ = "$Jul 28, 2010 2:32:23 PM$"


