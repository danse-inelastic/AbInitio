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

# Improvements:
# 1. Output the status of the simulation (as in ITaskApp.py)
# 2. (?) Log the status

from pyre.applications.Application import Application

class JobManager(Application):
    
    class Inventory(Application.Inventory):
        import pyre.inventory
        import jobmanager.components

        # Set dds and csaccessor parameters
        dds = pyre.inventory.facility(name="dds", factory=jobmanager.components.dds)
        dds.meta['tip'] = "the component manages data files"

        csaccessor = pyre.inventory.facility(name='csaccessor', factory = jobmanager.components.ssher)
        csaccessor.meta['tip'] = 'computing server accessor'

    def main(self):
        print "Hello world!"

    def __init__(self, name):
        super(JobManager, self).__init__(self, name)

    def _configure(self):
        super(JobManager, self)._configure()
        self.dds = self.inventory.dds
        self.dds.director = self
        self.csaccessor = self.inventory.csaccessor

    def _init(self):
        super(JobManager, self)._init()




__date__ = "$Oct 21, 2009 7:11:53 AM$"


