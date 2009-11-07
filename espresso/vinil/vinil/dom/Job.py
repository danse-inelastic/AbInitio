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

# "userId": 1 -> dexity

from vinil.utils.const import STATES
from vinil.utils.utils import timestamp
from vinil.components.DBTable import DBTable

class Job(DBTable):
    # 'name' attribute should be present in every class table.
    name = "job"
    import pyre.db
    
    id          = pyre.db.varchar(name="id", length=8)
    id.constraints = 'PRIMARY KEY'
    id.meta['tip'] = "the unique id"

    userId      = pyre.db.varchar(name="userId", length=8)
    userId.constraints = 'REFERENCES users (id)'
    userId.meta['tip'] = ""

    simulationId = pyre.db.varchar(name="simulationId", length=8)
    simulationId.constraints = 'REFERENCES simulation (id)'
    simulationId.meta['tip'] = ""

    serverId    = pyre.db.varchar(name="serverId", length=8)
    serverId.constraints = 'REFERENCES server (id)'
    serverId.meta['tip'] = ""

    description = pyre.db.varchar(name="description", length=1024, default='')
    description.meta['tip'] = ""

    status = pyre.db.varchar(name="status", length=64, default='')
    status.meta['tip'] = ""

    timeSubmitted = pyre.db.varchar(name="timeSubmitted", length=16, default='')
    timeSubmitted.meta['tip'] = "timeSubmitted"

    timeStarted = pyre.db.varchar(name="timeStarted", length=16, default='')
    timeStarted.meta['tip'] = "timeStarted"

    timeRestarted = pyre.db.varchar(name="timeRestarted", length=16, default='')
    timeRestarted.meta['tip'] = "timeRestarted"

    timeCompleted = pyre.db.varchar(name="timeCompleted", length=16, default='')
    timeCompleted.meta['tip'] = "timeCompleted"

    exitCode = pyre.db.integer(name="exitCode", default=-1)
    exitCode.meta['tip'] = "exitCode"

    numberProcessors = pyre.db.integer(name="numberProcessors", default=0)
    numberProcessors.meta['tip'] = "numberProcessors"

    errorFilename = pyre.db.varchar(name="errorFilename", length=256, default='stderr.log')
    errorFilename.meta['tip'] = "errorFilename"

    outputFilename = pyre.db.varchar(name="outputFilename", length=256, default='stdout.log')
    outputFilename.meta['tip'] = "outputFilename"

    statusMessage = pyre.db.varchar(name="statusMessage", length=256, default='')
    statusMessage.meta['tip'] = "statusMessage"


# Default records
defaults    = ({"id": 1, "userId": 1, "simulationId": 4, "description": "",
               "status": STATES["C"], "timeCompleted": timestamp() + 60, "exitCode": 0,
                "numberProcessors": 8},
               {"id": 2, "userId": 1, "simulationId": 5, "description": "",
               "status": STATES["C"], "timeCompleted": timestamp() + 60, "exitCode": 0,
                "numberProcessors": 8},
               {"id": 3, "userId": 1, "simulationId": 6, "description": "",
               "status": STATES["R"], "timeCompleted": timestamp() + 60, "exitCode": 0,
                "numberProcessors": 8})

# Init tables
def inittable(clerk):
    for params in defaults:
        r   = Job()
        r.setClerk(clerk)
        r.createRecord(params)


def test():
    for e in defaults:
        s = ""
        for v in e.keys():
            s += "%s: %s " % (v, e[v])
        print s

if __name__ == "__main__":
    test()


# ******************* DEAD CODE ******************
#   {"id": 1, "userId": 1, "simulationId": 4, "description": "",
#   "status": STATES["C"], "timeSubmitted":, timeStarted: ,
#    timeRestarted: , timeCompleted: , exitCode: ,
#    numberProcessors: , errorFilename: , outputFilename: ,
#    statusMessage


__date__ = "$Jul 29, 2009 8:31:54 PM$"

