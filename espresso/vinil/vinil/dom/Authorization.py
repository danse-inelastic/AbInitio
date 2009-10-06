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

    serverId    = pyre.db.varchar(name="serverId", length=8, default='')
    serverId.constraints = 'REFERENCES server (id)'
    serverId.meta['tip'] = "serverId"

    userId    = pyre.db.varchar(name="userId", length=8, default='')
    userId.constraints = 'REFERENCES user (id)'
    userId.meta['tip'] = "userId"

    aname = pyre.db.varchar(name="aname", length=256, default='')
    aname.meta['tip'] = "aname"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"


"""
# For debugging
def inittable(db):
    def server(params):
        r           = Server()
        r.id        = params['id']
        r.serverId  = params['serverId']
        r.userId    = params['userId']
        return r

    records = [
        server( {"id": 1, "serverId": 1, "userId": 1})
        ]
    for r in records: db.insertRow( r )
    return
"""

__date__ = "$Oct 6, 2009 12:18:36 PM$"


