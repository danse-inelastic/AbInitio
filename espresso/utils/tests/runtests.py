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

# !!! FINISH

# Tests
def testCreateConfig():
    print "Testing creation of config file"
    qe  = QEInput()
    nl  = Namelist('control')
    nl.add('title', "'Ni'")
    nl.add('restart_mode', "'from_scratch'")
    print "Adding parameters to namelist:\n%s" % nl.toString()
    nl.set('title', "'Fe'")
    qe.addNamelist(nl)
    print "Adding namelist to QEInput:\n%s" % qe.toString()

    c = Card('atomic_species')
    c.addLine('Ni  26.98  Ni.pbe-nd-rrkjus.UPF')
    print "Adding line to card:\n%s" % c.toString()
    qe.addCard(c)
    print "Adding card to QEInput:\n%s" % qe.toString()
    #qe.save()


def testParseConfig():
    print "Testing parsing config file"
    qe  = QEInput("../tests/ni.scf.in")
    qe.parse()
    print qe.toString()
    nl  = qe.namelist('control')
    nl.add('title', 'Ni')
    nl.remove('restart_mode')
    qe.removeCard('atomic_species')
    nl.set('calculation', "'nscf'")
    c = qe.card('atomic_positions')
    c.editLines(['Say Hi! :)'])
    print qe.toString()
    #qe.save("../tests/ni.scf.in.mod")

def testAttach():
    qe  = QEInput("../tests/si.ph.in", type="ph")
    qe.parse()
    qe.save("../tests/si.ph.in.mod")


if __name__ == "__main__":
    #testCreateConfig()
    #testParseConfig()
    testAttach()

__date__ = "$Jul 28, 2010 2:32:23 PM$"


