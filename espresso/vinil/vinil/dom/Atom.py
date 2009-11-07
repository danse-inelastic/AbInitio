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
Atom  - table that contains atom data (not sure if I need it)
"""

from vinil.components.DBTable import DBTable

class Atom(DBTable):

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

    aname = pyre.db.varchar(name="aname", length=256, default='')
    aname.meta['tip'] = "aname"

    mass = pyre.db.double(name="mass", default=0.0)
    mass.meta['tip'] = "mass"

    massUnit = pyre.db.varchar(name="massUnit", length=32, default='')
    massUnit.meta['tip'] = "massUnit"

    pseudoPotential = pyre.db.varchar(name="pseudoPotential", length=256, default='')
    pseudoPotential.meta['tip'] = "pseudoPotential"

    position = pyre.db.varchar(name="position", length=256, default='')
    position.meta['tip'] = "position"


# Default records
defaults    = ()

# Init tables
def inittable(clerk):
    for params in defaults:
        r   =Atom()
        r.setClerk(clerk)
        r.createRecord(params)

__date__ = "$Oct 5, 2009 8:46:22 AM$"


