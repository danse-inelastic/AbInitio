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
VTable  - not sure why I need this table
"""

from vinil.components.DBTable import DBTable

class VTable(DBTable):

    name = "vtable"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    tname = pyre.db.varchar(name="tname", length=256, default='')
    tname.meta['tip'] = "tname"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"


# Default records
defaults    = ()

# Init tables
def inittable(clerk):
    for params in defaults:
        r   = VTable()
        r.setClerk(clerk)
        r.createRecord(params)

__date__ = "$Oct 5, 2009 9:22:07 AM$"


