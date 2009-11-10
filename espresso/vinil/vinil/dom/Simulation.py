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

"""
Simulation  - table that contains simulation data
"""
#from vinil.utils.const import SIMULATIONS

from vinil.components.DBTable import DBTable

class Simulation(DBTable):

    name = "simulation"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    sname = pyre.db.varchar(name="sname", length=128, default='')
    sname.meta['tip'] = ""

    package = pyre.db.varchar(name="package", length=128, default='')
    package.meta['tip'] = ""

    type = pyre.db.varchar(name="type", length=128, default='')
    type.meta['tip'] = ""

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = ""

    formula = pyre.db.varchar(name="formula", length=32, default='')
    formula.meta['tip'] = ""

    timeCreated = pyre.db.varchar(name="timeCreated", length=16, default='')
    timeCreated.meta['tip'] = "timeCreated"

    timeModified = pyre.db.varchar(name="timeModified", length=16, default='')
    timeModified.meta['tip'] = "timeModified"

    isFavorite = pyre.db.boolean(name="isFavorite", default=True)   #?
    isFavorite.meta['tip'] = ""

    isExample = pyre.db.boolean(name="isExample", default=False)
    isExample.meta['tip'] = ""



# Default records
defaults    = ({"id": 1, "sname": 'MgB2_SP', "package": 'Quantum Espresso',
                "type": 'Single-Phonon', "description": 'Single-Phonon simualtion',
                "formula": 'MgB2'},
                {"id": 2, "sname": 'MgB2_E', "package": 'Quantum Espresso',
                "type": 'Total Energy', "description": 'Electron simualtion',
                "formula": 'MgB2'},
                {"id": 3, "sname": 'MgB2_MP', "package": 'Quantum Espresso',
                "type": 'Multi-Phonon DOS', "description": 'Multy-Phonon simualtion',
                "formula": 'MgB2'},
                {"id": 4, "sname": 'Ni_Energy', "package": 'Quantum Espresso',
                "type": 'Total Energy', "description": 'Total Energy simualtion',
                "formula": 'Ni', "isFavorite": False, "isExample": True},
                {"id": 5, "sname": 'Ni_E_DOS', "package": 'Quantum Espresso',
                "type": 'Electron DOS', "description": 'Electron DOS simualtion',
                "formula": 'Ni', "isFavorite": False, "isExample": True},
                {"id": 6, "sname": 'Ni_Ph_DOS', "package": 'Quantum Espresso',
                "type": 'Multi-Phonon DOS', "description": 'Multy-Phonon DOS simualtion',
                "formula": 'Ni', "isFavorite": False, "isExample": True})

# Init tables
def inittable(clerk):
    for params in defaults:
        r   = Simulation()
        r.setClerk(clerk)
        r.createRecord(params)


def test():
    for e in examples:
        s = ""
        for v in e:
            s += "%s " % v
        print s

if __name__ == "__main__":
    test()



__date__ = "$Oct 5, 2009 8:53:44 AM$"


