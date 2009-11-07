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

"""
Matter  - table that contains data related to matter
"""

from vinil.components.DBTable import DBTable

class Matter(DBTable):

    name = "matter"
    import pyre.db

    id = pyre.db.varchar(name="id", length=10)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    simulationId    = pyre.db.varchar(name="simulationId", length=8)
    simulationId.constraints = 'REFERENCES simulation (id)'
    simulationId.meta['tip'] = "simulationId"

    mname = pyre.db.varchar(name="mname", length=1024, default='')
    mname.meta['tip'] = "mname"

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = "description"


# Default records
defaults    = ()

# Init tables
def inittable(clerk):
    for params in defaults:
        r   = Matter()
        r.setClerk(clerk)
        r.createRecord(params)
        
__date__ = "$Oct 5, 2009 8:50:39 AM$"


