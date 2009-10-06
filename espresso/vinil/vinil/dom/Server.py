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

class Server(Table):

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


"""
# For debugging
def inittable(db):
    def server(params):
        r           = Server()
        r.id        = params['id']
        return r

    records = [
        server( {"id": 1})
        ]
    for r in records: db.insertRow( r )
    return
"""

__date__ = "$Oct 5, 2009 8:10:57 PM$"


