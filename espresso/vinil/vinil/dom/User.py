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

defaults    = ({"id": 1, "username": "dexity", "firstName": "Alex",
               "lastName": "Dementsov", "email": "somemail@gmail.com",
               "affiliation": "CalTech", "password": "5f4dcc3b5aa765d61d8327deb882cf99"}, # 'passowrd' -> md5
              )

# For some reason when I try to use "user" name for database table, it returns error
# Probably conflict with internal code?
from pyre.db.Table import Table

class User(Table):

    name = "users"  # Table name
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

    timeCreated = pyre.db.varchar(name="timeCreated", length=16, default='')
    timeCreated.meta['tip'] = "timeCreated"

    timeModified = pyre.db.varchar(name="timeModified", length=16, default='')
    timeModified.meta['tip'] = "timeModified"

    password = pyre.db.varchar(name="password", length=256, default='')
    password.meta['tip'] = "password"


    def updateRecord(self, director, params):
        """
        Updates configuration row (even if key in params is not present).
        'id' ans 'timeCreated' are not updated!
        """
        self.username       = setname(params, self, 'username')
        self.firstName      = setname(params, self, 'firstName')
        self.lastName       = setname(params, self, 'lastName')
        self.email          = setname(params, self, 'email')
        self.affiliation    = setname(params, self, 'affiliation')
        self.timeModified   = timestamp()
        self.password       = setname(params, self, 'password')

        director.clerk.updateRecord(self)   # Update record


    def createRecord(self, director, params):
        """Inserts user row """

        self.id             = newid(director)
        self.username       = setname(params, self, 'username')
        self.firstName      = setname(params, self, 'firstName')
        self.lastName       = setname(params, self, 'lastName')
        self.email          = setname(params, self, 'email')
        self.affiliation    = setname(params, self, 'affiliation')
        self.timeCreated    = timestamp()
        self.timeModified   = timestamp()
        self.password       = setname(params, self, 'password')

        director.clerk.insertRecord(self)


    def deleteRecord(self, director):
        """Deletes record"""
        director.clerk.deleteRecord(self, id=self.id)


# For debugging
def inittable(db):

    def user(params):
        r                = User()
        r.id             = params['id']
        r.username       = params['username']
        r.firstName      = params['firstName']
        r.lastName       = params['lastName']
        r.email          = params['email']
        r.affiliation    = params['affiliation']
        r.password       = params['password']

        return r

    records = []
    for e in defaults:
        records.append(user(e))

    for r in records: db.insertRow( r )
    return

def test():
    for e in defaults:
        s = ""
        for v in e.keys():
            s += "%s: %s " % (v, e[v])
        print s

if __name__ == "__main__":
    test()


__date__ = "$Oct 5, 2009 8:11:34 PM$"


