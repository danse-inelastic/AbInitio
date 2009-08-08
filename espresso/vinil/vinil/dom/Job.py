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

class Job(Table):

    name = "job"
    import pyre.db
    
    id = pyre.db.varchar(name="id", length=10)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    type = pyre.db.varchar(name="type", length=100)
    type.meta['tip'] = "Type of the simulation (electron/phonon)"

    status = pyre.db.varchar(name="status", length=100)
    status.meta['tip'] = "Status of the job (not started/running/finished)"

    created = pyre.db.varchar(name="created", length=100)
    created.meta['tip'] = "Date and time when the simulation started"

    config = pyre.db.varchar(name="config", length=5000)
    config.meta['tip'] = "Configuration text"

    pass # end of DbObject

__date__ = "$Jul 29, 2009 8:31:54 PM$"


