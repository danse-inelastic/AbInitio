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

from pyre.db.Table import Table

examples = (
            (1, 'MgB2_SP', 'Quantum Espresso', 'Single-Phonon', 'Single-Phonon simualtion', 'MgB2', '25-09-2009', '', True, False),
            (2, 'MgB2_E', 'Quantum Espresso', 'Total Energy', 'Electron simualtion', 'MgB2', '26-09-2009', '', True, False),
            (3, 'MgB2_MP', 'Quantum Espresso', 'Multi-Phonon', 'Multy-Phonon simualtion', 'MgB2', '27-09-2009', '', True, False),
            (4, 'Ni_SP', 'Quantum Espresso', 'Single-Phonon', 'Single-Phonon simualtion', 'Ni', '25-09-2009', '', False, True),
            (5, 'Al_E', 'Quantum Espresso', 'Total Energy', 'Electron simualtion', 'Al', '26-09-2009', '', False, True),
            (6, 'Si_MP', 'Quantum Espresso', 'Multi-Phonon', 'Multy-Phonon simualtion', 'Si', '27-09-2009', '', False, True)
            )

class Simulation(Table):

    name = "simulation"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    sname = pyre.db.varchar(name="name", length=128)
    sname.meta['tip'] = ""

    package = pyre.db.varchar(name="package", length=128)
    package.meta['tip'] = ""

    type = pyre.db.varchar(name="type", length=128)
    type.meta['tip'] = ""

    description = pyre.db.varchar(name="description", length=1024)
    description.meta['tip'] = ""

    formula = pyre.db.varchar(name="formula", length=32)
    formula.meta['tip'] = ""

    created = pyre.db.varchar(name="created", length=16)
    created.meta['tip'] = ""

    modified = pyre.db.varchar(name="modified", length=16)
    modified.meta['tip'] = ""

    isFavorite = pyre.db.boolean(name="isFavorite", default=False)
    isFavorite.meta['tip'] = ""

    isExample = pyre.db.boolean(name="isExample", default=True)
    isExample.meta['tip'] = ""


def inittable(db):
    def simulation(id, sname, package, type, description, formula, created, modified, isFavorite, isExample):
        r               = Simulation()
        r.id            = id
        r.sname         = sname
        r.package       = package
        r.type          = type
        r.description   = description
        r.formula       = formula
        r.created       = created
        r.modified      = modified
        r.isFavorite    = isFavorite
        r.isExample     = isExample
        return r


    records = []
    for e in examples:
        records.append(simulation(id        = e[0],
                                  sname     = e[1],
                                  package   = e[2],
                                  type      = e[3],
                                  description = e[4],
                                  formula   = e[5],
                                  created   = e[6],
                                  modified  = e[7],
                                  isFavorite = e[8],
                                  isExample = e[9]))

    for r in records: db.insertRow( r )
    return


def test():
    for e in examples:
        s = ""
        for v in e:
            s += "%s " % v
        print s

if __name__ == "__main__":
    test()



__date__ = "$Oct 5, 2009 8:53:44 AM$"


