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

class Matter(Table):

    name = "matter"
    import pyre.db

    id = pyre.db.varchar(name="id", length=10)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    simulationId    = pyre.db.varchar(name="simulationId", length=8)
    simulationId.constraints = 'REFERENCES simulation (id)'
    simulationId.meta['tip'] = "simulationId"

    name = pyre.db.varchar(name="name", length=1024, default='')
    name.meta['tip'] = "name"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"


"""
# For debugging
def inittable(db):
    def matter(params):
        r               = Matter()
        r.id            = params['id']
        r.simulationId  = params['simulationId']
        return r

    records = [
        matter( {"id": 1, "simulationId": 1} )
        ]
    for r in records: db.insertRow( r )
    return
"""

__date__ = "$Oct 5, 2009 8:50:39 AM$"


