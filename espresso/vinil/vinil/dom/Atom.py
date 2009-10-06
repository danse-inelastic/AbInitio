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

class Atom(Table):

    name = "atom"
    import pyre.db

    id          = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    matterId    = pyre.db.varchar(name="matterId", length=8)
    matterId.constraints = 'REFERENCES matter (id)'
    matterId.meta['tip'] = ""

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"

    name = pyre.db.varchar(name="name", length=256, default='')
    name.meta['tip'] = "name"

    mass = pyre.db.double(name="mass", dafault=0.0)
    mass.meta['tip'] = "mass"

    massUnit = pyre.db.varchar(name="massUnit", length=32, default='')
    massUnit.meta['tip'] = "massUnit"

    pseudoPotential = pyre.db.varchar(name="pseudoPotential", length=256, default='')
    pseudoPotential.meta['tip'] = "pseudoPotential"

    position = pyre.db.varchar(name="position", length=256, default='')
    position.meta['tip'] = "position"


# Do I need generate default jobs? For debugging, yes!
"""
def inittable(db):
    def atom(params):
        r           = Atom()
        r.id        = params['id']
        r.matterId  = params['matterId']
        return r

    records = [
        atom( {"id": 1, "matterId": 1})
        ]
    for r in records: db.insertRow( r )
    return
"""

__date__ = "$Oct 5, 2009 8:46:22 AM$"


