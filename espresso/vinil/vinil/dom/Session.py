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
Session  - table that contains session data
"""

from vinil.components.DBTable import DBTable

class Session(DBTable):

    name = "session"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    sname = pyre.db.varchar(name="sname", length=256, default='')
    sname.meta['tip'] = "sname"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"

# Default records
defaults    = ()

# Init tables
def inittable(clerk):
    for params in defaults:
        r   =Session()
        r.setClerk(clerk)
        r.createRecord(params)

__date__ = "$Oct 5, 2009 8:11:18 PM$"


