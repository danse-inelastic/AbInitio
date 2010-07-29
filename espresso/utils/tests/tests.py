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

# Before to run tests, please add "../parser" to PYTHON_PATH environmental variable

from namelist import Namelist
from card import Card
from qeparser import QEParser
from qeinput import QEInput

import unittest

class QEParserTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_namelist_name(self):
        nl  = Namelist("control")
        self.assertEqual(nl.name(), "control")


if __name__ == '__main__':
    unittest.main()
    
# namelist
# card
# qeparse
# qeinput



#textPW = """
# &control
#    calculation='scf'
#    restart_mode='from_scratch',
#    tprnfor = .true.
#    prefix='ni',
#    pseudo_dir = '',
#    outdir=''
# /
# &system
#    ibrav=2,
#    celldm(1) =6.65,
#    nat=  1,
#    ntyp= 1,
#    nspin=2,
#    starting_magnetization(1)=0.5,
#    degauss=0.02,
#    smearing='gauss',
#    occupations='smearing',
#    ecutwfc =27.0
#    ecutrho =300.0
# /
# &electrons
#    conv_thr =  1.0d-8
#    mixing_beta = 0.7
# /
#
#
#ATOMIC_SPECIES
# Ni  26.98  Ni.pbe-nd-rrkjus.UPF
#
#ATOMIC_POSITIONS
# Ni 0.00 0.00 0.00
#K_POINTS AUTOMATIC
#4 4 4 1 1 1
#blah
#
#"""
#
#textCards = """
#ATOMIC_POSITIONS
# Ni 0.00 0.00 0.00
#K_POINTS AUTOMATIC
#4 4 4 1 1 1
#blah
#
#"""
#
#textProblem = """
#CELL_PARAMETERS
#   0.993162743  -0.000000000   0.000000000
#  -0.496581371   0.860104165  -0.000000000
#  -0.000000000  -0.000000000   4.345938530
#ATOMIC_POSITIONS
# Ni 0.00 0.00 0.00
#K_POINTS AUTOMATIC
#4 4 4 1 1 1
#blah
#
#"""
#
#textMatdyn = """
# &input
#    asr='crystal',
#    amass(1)=24.305, amass(2)=11.000,
#    flfrc='mgb2666.fc'
# /
#176
#0.000000    0.000000    0.456392    0.000000
#0.000000    0.000000    0.447264    0.009128
#0.000000    0.000000    0.438137    0.018256
#0.000000    0.000000    0.429009    0.027384
#0.000000    0.000000    0.419881    0.036511
#"""
#
#textDynmat = """
#&input  fildyn='mgb2.dynG', asr='simple',
#        q(1)=0.0, q(2)=0.0, q(3)=0.0 /
#"""
#
#textPh  = """
# &inputph
#  tr2_ph=1.0d-10,
#  amass(1)=24.305,
#  amass(2)=11.000,
#  prefix='mgb2',
#  outdir='/scratch/markovsk/mgb2'
#  fildyn='mgb2.dynG',
# /
#
#"""
#
#textHeader  = """
#&INPUTPH
#   tr2_ph = 1.0d-12,
#   prefix = 'si',
#   epsil = .false.,
#   trans = .true.,
#   zue = .false.,
#   outdir = '/scratch/si',
#   amass(1) = 28.0855,
#   fildyn = 'si.dyn_G',
#   fildrho = 'si.drho_G',
#/
#0.0 0.0 0.0
#"""
#
## This is not a problem text (just add spaces between commas)
#textComma   = """&input
#   asr='crystal',  dos=.true.
#   amass(1)=26.982538, amass(2)=11.000,
#   flfrc='mgalb4666.fc', fldos='mgalb4.666.phdos', nk1=28,nk2=28,nk3=28
#/
#"""
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


