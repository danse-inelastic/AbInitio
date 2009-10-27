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


from pyre.components.Component import Component

class Worker(Component):

    class Inventory(Component.Inventory):
        pass

    def __init__(self, director):
        self._director  = director

    # Add additional parameter to specify the task
    def run(self):
        job = None  # Populate job from task

        try:
            self.prepare(job)
            #self.schedule(job)
        except Exception, e:
            import traceback
            #self._debug.log('submission of Job failed. %s' % traceback.format_exc())
            raise

    def prepare(self, job):
#        dds = self._director.dds
#
#        jobpath = dds.abspath(job)

#        #At this point you should know about the simulation you want to run
#        from jobmanager.simulations.QE import QE
#        files, deps = buildjob(computation, db=clerk.db, dds=dds, path=jobpath, director=director)
#        for f in files:
#            dds.remember(job, f)
#
#        # make job related files available on the server
#        dds.make_available(job, server=server, files=files)

        # make dependencies available on the server

        #for dep in deps:
        #    self.prepare_dependency(dep, job)
        return

    def makeAvailable(self):
        pass

    def prepare_dependency(self, dep, job):
#        dds = self._director.dds
#
#        if dds.is_available(dep):
#            dds.remember(dep)
#
#        server = None
#        dds.make_available(dep, server=server)
        return


    def schedule(self, job):
        from jobmanager.components.Scheduler import Scheduler
        s   = Scheduler(self._director)     # job,
        return s.schedule()


if __name__ == "__main__":
    w   = Worker(None)
    w.run()



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


