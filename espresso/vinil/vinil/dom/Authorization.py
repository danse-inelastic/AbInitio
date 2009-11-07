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
Authorization  - table that contains authorization data
"""

from vinil.components.DBTable import DBTable

class Authorization(DBTable):

    name = "authorization"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    serverId    = pyre.db.varchar(name="serverId", length=8, default='')
    #serverId.constraints = 'REFERENCES server (id)'
    serverId.meta['tip'] = "serverId"

    userId    = pyre.db.varchar(name="userId", length=8, default='')
    #userId.constraints = 'REFERENCES user (id)'
    userId.meta['tip'] = "userId"

    aname = pyre.db.varchar(name="aname", length=256, default='')
    aname.meta['tip'] = "aname"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"


# Default records
defaults    = ()

# Init tables
def inittable(clerk):
    for params in defaults:
        r   =Authorization()
        r.setClerk(clerk)
        r.createRecord(params)

__date__ = "$Oct 6, 2009 12:18:36 PM$"


