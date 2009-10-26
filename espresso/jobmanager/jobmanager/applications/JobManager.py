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

# No database involved - just config files!

# Improvements:
# 1. Output the status of the simulation (as in ITaskApp.py)
# 2. (?) Log the status

# TODO:
# Where to put ComputationResultsRetriever.py?

from pyre.applications.Script import Script

class JobManager(Script):

    class Inventory(Script.Inventory):
        import pyre.inventory        
        import jobmanager.components

        # Set dds and csaccessor parameters
        dds         = pyre.inventory.facility(name="dds", factory=jobmanager.components.dds)
        dds.meta['tip'] = "the component manages data files"

        csaccessor  = pyre.inventory.facility(name='csaccessor', factory = jobmanager.components.ssher)
        csaccessor.meta['tip'] = 'computing server accessor'

        serverip      = pyre.inventory.str(name="serverip", default="127.0.0.1")
        
        servername    = pyre.inventory.str(name="servername", default="localhost")
        servername.meta['label'] = 'Computation server'
        servername.meta['tip'] = ('Please choose the server on which the job will be run')

        numprocessors = pyre.inventory.str( 'numprocessors', default = 1 )
        numprocessors.meta['label'] = 'Number of processors'
        numprocessors.meta['tip'] = ('Please input the number of processors')
        numprocessors.meta['tiponerror'] = ('Please enter a positive integer')

        #?
        walltime = pyre.inventory.str( 'walltime', default = 10)
        walltime.meta['label'] = 'Time limit (hours)'
        walltime.meta['tip'] = ('Please input a limit on the time your job will run. (Unit: hours)')
        walltime.meta['tiponerror'] = ('Please enter a positive integer')


    # Here is where the essence of the script goes
    def main(self):
        # Start time, Finish time
        from jobmanager.components.Worker import Worker

        # Need to pass parameters of the
        d   = Worker(self)
        d.run()
        return

    def __init__(self, name=None):
        super(JobManager, self).__init__(name=name)

    def _configure(self):
        super(JobManager, self)._configure()
        
        self.dds        = self.inventory.dds
        self.dds.director = self
        self.csaccessor = self.inventory.csaccessor
        self.servername     = self.inventory.servername

    def _init(self):
        super(JobManager, self)._init()

    def _getPrivateDepositoryLocations(self):
        return ['../config']


if __name__ == "__main__":
    app = JobManager(name="main")
    app.run()


__date__ = "$Oct 21, 2009 7:11:53 AM$"


