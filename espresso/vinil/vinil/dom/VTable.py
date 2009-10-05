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

    id = pyre.db.varchar(name="id", length=10)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"


def inittable(db):
    def vtable(id):
        r           = VTable()
        r.id        = id
        return r

    records = [
        vtable( 1 )
        ]
    for r in records: db.insertRow( r )
    return


__date__ = "$Oct 5, 2009 9:22:07 AM$"


