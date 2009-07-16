from configParser import *
from configPW import *
from MatterBase import MatterBase

def testCreateConfig():
    print "Testing creation of config file"
    qe  = QEConfig()
    nl  = Namelist('control')
    nl.addParam('title', "'Ni'")
    nl.addParam('restart_mode', "'from_scratch'")
    print "Adding parameters to namelist:\n%s" % nl.toString()
    nl.editParam('title', "'Fe'")
    qe.addNamelist(nl)
    print "Adding namelist to QEConfig:\n%s" % qe.toString()

    c = Card('atomic_species')
    c.addLine('Ni  26.98  Ni.pbe-nd-rrkjus.UPF')
    print "Adding line to card:\n%s" % c.toString()
    qe.addCard(c)
    print "Adding card to QEConsig:\n%s" % qe.toString()
    qe.save()

if __name__ == "__main__":
    testCreateConfig()


