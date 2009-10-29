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
JobManager manages jobs on the remote cluster

Specification:
1. No database involved - just config files!

Improvements:
1. Output the status of the simulation (as in ITaskApp.py)
2. [?] Log the output
3. [?] Store files and directories manipulation in Distributed Data Storage

Questions:
What is the ComputationResultsRetriever.py for?

TODO:
[?] Set limitation for the simulation time (start, finish time)
"""

import os
from pyre.applications.Script import Script

class JobManager(Script):

    class Inventory(Script.Inventory):
        import pyre.inventory        
        import jobmanager.components

        csaccessor  = pyre.inventory.facility(name='csaccessor', factory = jobmanager.components.ssher)
        csaccessor.meta['tip'] = 'computing server accessor'

        settings     = pyre.inventory.str(name="settings", default=None)
        settings.meta['label'] = 'Settings for simulation parameters'

    # Main method!
    def main(self):
        from jobmanager.components.Worker import Worker
        
        d   = Worker(self)      # Need to pass parameters?
        d.run()
        return

    def __init__(self, name=None):
        super(JobManager, self).__init__(name=name)

    def _configure(self):
        super(JobManager, self)._configure()
        
        self.csaccessor = self.inventory.csaccessor
        self.settings   = self.inventory.settings
        if self.settings is None:
            print """
Settings configuration file should be provided!
Usage: jm.py --settings=<filename>
"""
            raise

    def _init(self):
        super(JobManager, self)._init()

    def _getPrivateDepositoryLocations(self):
        return ['../config']


if __name__ == "__main__":
    app = JobManager(name="main")
    app.run()


__date__ = "$Oct 21, 2009 7:11:53 AM$"






# ********** DEAD CODE ****************************************

#        Set dds and csaccessor parameters
#        dds         = pyre.inventory.facility(name="dds", factory=jobmanager.components.dds)
#        dds.meta['tip'] = "the component manages data files"

#        server     = pyre.inventory.str(name="server", default=thisfile+"/../../config/foxtrot.danse.us.conf")
#        server.meta['label'] = 'Settings related to server'
#
#        servername    = pyre.inventory.str(name="servername", default="localhost")
#        servername.meta['label'] = 'Computation server'
#        servername.meta['tip'] = ('Please choose the server on which the job will be run')

#        serverip      = pyre.inventory.str(name="serverip", default="127.0.0.1")
#
#        numprocessors = pyre.inventory.str( 'numprocessors', default = 1 )
#        numprocessors.meta['label'] = 'Number of processors'
#        numprocessors.meta['tip'] = ('Please input the number of processors')
#        numprocessors.meta['tiponerror'] = ('Please enter a positive integer')
#
#        #?
#        walltime = pyre.inventory.str( 'walltime', default = 10)
#        walltime.meta['label'] = 'Time limit (hours)'
#        walltime.meta['tip'] = ('Please input a limit on the time your job will run. (Unit: hours)')
#        walltime.meta['tiponerror'] = ('Please enter a positive integer')


