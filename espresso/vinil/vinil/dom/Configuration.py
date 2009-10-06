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

class Configuration(Table):

    name = "configuration"
    import pyre.db

    id = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    simulationId    = pyre.db.varchar(name="simulationId", length=8)
    simulationId.constraints = 'REFERENCES simulation (id)'
    simulationId.meta['tip'] = "simulationId"

    name = pyre.db.varchar(name="name", length=1024, default='')
    name.meta['tip'] = "name"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"

    date = pyre.db.varchar(name="date", length=16, default='')
    date.meta['tip'] = "date"

    text = pyre.db.varchar(name="text", length=8192, default='')
    text.meta['tip'] = "text"


"""
# For debugging
def inittable(db):
    def configuration(id):
        r           = Configuration()
        r.id        = params['id']
        r.simulationId  = params['simulationId']
        return r

    records = [
        configuration( {"id": 1, "simulationId": 1} )
        ]
    for r in records: db.insertRow( r )
    return
"""

__date__ = "$Oct 5, 2009 8:58:32 AM$"


