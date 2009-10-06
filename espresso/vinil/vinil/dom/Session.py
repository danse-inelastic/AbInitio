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

class Session(Table):

    name = "session"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    sname = pyre.db.varchar(name="sname", length=256, default='')
    sname.meta['tip'] = "sname"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"


"""
# For debugging
def inittable(db):
    def session(params):
        r           = Session()
        r.id        = params['id']
        return r

    records = [
        session( {"id": 1} )
        ]
    for r in records: db.insertRow( r )
    return"""
"""

__date__ = "$Oct 5, 2009 8:11:18 PM$"


