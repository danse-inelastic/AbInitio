#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                  Jiao Lin
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

#    def run(self):
#        # 1. Read configuration
#        # 2. Create thread (syncronous)
#        # 3. Execute command on remote cluster
#        # 4. Get status of the task
"""
History: Adapter from submitjob.odb
"""

import os
from jobmanager.utils.Server import Server
from pyre.components.Component import Component
import ConfigParser

class Worker(Component):

    class Inventory(Component.Inventory):
        pass

    def __init__(self, director):
        self._director  = director
        self._jobid     = director.jobid
        self._ssher     = director.csaccessor
        self._settings  = ConfigParser.ConfigParser()
        self._settings.read(director.settings)
        self._localpath     = self._settings.get("simulation", "local-path")
        self._remotepath    = self._settings.get("simulation", "remote-path")
        self._input         = self._settings.get("simulation", "input-file")
        self._output        = self._settings.get("simulation", "output-file")
        self._username      = self._settings.get("user", "username")
        self._jobname       = self._settings.get("simulation", "job-name")

    # Add additional parameter to specify the task
    def run(self):
        try:
            self.prepare(self._jobid)
            self.schedule(self._jobid)
            self.retrieveResults(self._jobid)   # specific for this example
            self.clean(self._jobid)
        except Exception, e:
            import traceback
            #self._debug.log('submission of Job failed. %s' % traceback.format_exc())
            raise

    def prepare(self, jobid):
        serverA = Server(None, None, self._username)
        serverB = Server(self._director.servername, None, self._username)   #
        self._ssher.mkdir(serverB, self._remotepath+"/temp" ) # Take output directory from config file
        self._ssher.copy(serverA, self._localpath+"/"+self._input, serverB, self._remotepath)
        self._ssher.copy(serverA, self._localpath+"/Ni.pbe-nd-rrkjus.UPF", serverB, self._remotepath)

    def schedule(self, jobid):
        from jobmanager.components.Scheduler import Scheduler
        s   = Scheduler(self._director)     # job,

        self._jobid = s.schedule()


    def retrieveResults(self, jobid):
        serverA = Server(None, None, self._username)
        serverB = Server(self._director.servername, None, self._username)   #

        #dir     = self._localpath+"/Ni.4969"
        dir     = self._localpath+"/%s.%s" % (self._jobname, self._jobid)
        os.mkdir(dir)
        import time
        time.sleep(1)
        self._ssher.copy(serverB, self._remotepath+"/"+self._input, serverA, dir)
        self._ssher.copy(serverB, self._remotepath+"/"+self._output, serverA, dir)


    def clean(self, jobid):
        """Cleans up after simulation (currently just removes the simulation directory) """
        serverB = Server(self._director.servername, None, self._username)
        self._ssher.rmdir(serverB, self._remotepath )


if __name__ == "__main__":
    """Cannot be called directly"""
    #w   = Worker(None)
    #w.run()



# **************** DEAD CODE *************************

#    def run(self, task):
#
#        director = self.director
#        clerk = director.clerk
#        db = clerk.db
#
#        director.declareProgress(0.1, 'Verifying job ...')
#
#        job = task.beneficiary.dereference(db)
#        id = job.id
#
#        state = job.state
#        if state not in ['created', 'submissionfailed']:
#            raise RuntimeError, "Job %s not suitable for submission: %s" % (id, state)
#
#        job.state = 'submitting'
#        clerk.updateRecordWithID(job)
#
#        try:
#            director.declareProgress(0.2, 'Verifying computation ...')
#            computation = job.computation
#            if not computation:
#                raise RuntimeError, 'computation is not specified for Job: %s' % (id,)
#
#            director.declareProgress(0.3, 'preparing the job ...')
#            self.prepare(job)
#
#            director.declareProgress(0.7, 'scheduling the job ...')
#            self.schedule(job)
#
#            director.declareProgress(1.0, 'done')
#
#        except Exception, e:
#            job.state = 'submissionfailed'
#            errmsg = '%s: %s' % (e.__class__.__name__, e)
#            job.error = errmsg
#            clerk.updateRecordWithID(job)
#
#            import traceback
#            self._debug.log('submission of Job %s failed. %s' % (
#                id, traceback.format_exc()) )
#
#            raise
#
#        return
#
#
#    def prepare(self, job):
#        director = self.director
#        dds = director.dds
#        clerk = director.clerk
#
#        jobpath = dds.abspath(job)
#        computation = clerk.dereference(job.computation)
#
#        director.declareProgress(0.4, 'Preparing: building job ...')
#        from vnfb.components.job_utils import buildjob
#        files, deps = buildjob(
#            computation, db=clerk.db, dds=dds, path=jobpath, director=director)
#        for f in files: dds.remember(job, f)
#
#        # make job related files available on the server
#        server = clerk.dereference(job.server)
#        director.declareProgress(0.5, 'Preparing: transfering job to the server %s ...' % server.short_description)
#        dds.make_available(job, server=server, files=files)
#
#        # make dependencies available on the server
#        director.declareProgress(0.6, 'Preparing: setting up dependencies ...')
#        for dep in deps: self.prepare_dependency(dep, job)
#        return
#
#
#    def prepare_dependency(self, dep, job):
#        director = self.director
#        dds = director.dds
#        clerk = director.clerk
#
#        if dds.is_available(dep): dds.remember(dep)
#
#        server = clerk.dereference(job.server)
#        director.declareProgress(0.65, 'Preparing: setting up dependency %s###%s ...' % (dep.__class__.__name__, dep.id) )
#        dds.make_available(dep, server=server)
#        return
#
#
#    def schedule(self, job):
#        from vnf.components.Scheduler import schedule
#        return schedule(job, self.director)
# Rename building job to running
# Call jobmanager.simulation.QE

#def buildjob() #computation, db=None, dds=None, path=None, director=None):
#    name = computation.__class__.__name__.lower()
#    builder = director.retrieveComponent(
#        name,
#        factory="job_builder", args=[name, path],
#        vault=['job_builders'])
#    files = builder.build(computation, db=db, dds=dds)
#    deps = builder.getDependencies()
#    return files, deps



__date__ = "$Oct 25, 2009 3:53:41 PM$"


