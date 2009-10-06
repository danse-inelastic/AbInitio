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

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    username = pyre.db.varchar(name="username", length=256, default='')
    username.meta['tip'] = "username"

    firstName = pyre.db.varchar(name="firstName", length=256, default='')
    firstName.meta['tip'] = "firstName"

    lastName = pyre.db.varchar(name="lastName", length=256, default='')
    lastName.meta['tip'] = "lastName"

    email  = pyre.db.varchar(name="email ", length=256, default='')
    email .meta['tip'] = "email "

    affiliation = pyre.db.varchar(name="affiliation", length=512, default='')
    affiliation.meta['tip'] = "affiliation"

    password = pyre.db.varchar(name="password", length=256, default='')
    password.meta['tip'] = "password"

"""
# For debugging
def inittable(db):
    def user(params):
        r           = User()
        r.id        = params['id']
        return r

    records = [
        user( {"id": 1} )
        ]
    for r in records: db.insertRow( r )
    return
"""

__date__ = "$Oct 5, 2009 8:11:34 PM$"


