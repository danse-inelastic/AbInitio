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

from pyre.db.Table import Table as base

class VTable(base):

    name = "vtable"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    tname = pyre.db.varchar(name="tname", length=256, default='')
    tname.meta['tip'] = "tname"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"


def inittable(db):
    def vtable(params):
        r           = VTable()
        r.id        = params['id']
        return r

    records = [
        vtable( {"id": 1} )
        ]
    for r in records: db.insertRow( r )
    return


__date__ = "$Oct 5, 2009 9:22:07 AM$"


