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

class User(Table):

    name = "user"
    import pyre.db

    id = pyre.db.varchar(name="id", length=10)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"


def inittable(db):
    def user(id):
        r           = User()
        r.id        = id
        return r

    records = [
        user( 1 )
        ]
    for r in records: db.insertRow( r )
    return

__date__ = "$Oct 5, 2009 8:11:34 PM$"

