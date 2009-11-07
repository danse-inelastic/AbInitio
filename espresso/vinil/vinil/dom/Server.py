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
Server  - table that contains server data
"""

from vinil.components.DBTable import DBTable

class Server(DBTable):

    name = "server"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    sname = pyre.db.varchar(name="sname", length=256)
    sname.meta['tip'] = "sname"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"

    ipAddress = pyre.db.varchar(name="ipAddress", length=64, default='')
    ipAddress.meta['tip'] = "ipAddress"


# Default records
defaults    = ()

# Init tables
def inittable(clerk):
    for params in defaults:
        r   = Server()
        r.setClerk(clerk)
        r.createRecord(params)

__date__ = "$Oct 5, 2009 8:10:57 PM$"


