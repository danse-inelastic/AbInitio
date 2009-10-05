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

class Simulation(Table):

    name = "simulation"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    name = pyre.db.varchar(name="name", length=128)
    name.meta['tip'] = ""

    kind = pyre.db.varchar(name="kind", length=128)
    kind.meta['tip'] = ""

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

    isFavorite = pyre.db.varchar(name="isFavorite", length=16)
    isFavorite.meta['tip'] = ""

    isExample = pyre.db.varchar(name="isExample", length=16)
    isExample.meta['tip'] = ""


def inittable(db):
    def simulation(id):
        r           = Simulation()
        r.id        = id
        return r

    records = [
        simulation( 1 )
        ]
    for r in records: db.insertRow( r )
    return


__date__ = "$Oct 5, 2009 8:53:44 AM$"


