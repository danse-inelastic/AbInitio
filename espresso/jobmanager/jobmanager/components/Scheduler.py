# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                   Jiao Lin
#                      California Institute of Technology
#                      (C) 2007-2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

"""
History: Adopted from Scheduler.py
"""

import os
import ConfigParser
from jobmanager.utils.Server import Server
from jobmanager.components.Torque import Torque

class Scheduler:

    def __init__(self, director):
        self._director  = director
        self._state      = None

        self.settings    = ConfigParser.ConfigParser()
        self.settings.read(self._director.settings)

        self._initScheduler()

    #TODO: Consider case when parameters in cmd are empty!
    # Rename parameters
    def schedule(self):
        simtype     = self.settings.get("simulation", "job-type")
        npool       = int(self.settings.get("server", "npool"))
        input       = self.settings.get("simulation", "input-file")
        output      = self.settings.get("simulation", "output-file")


        """ E.g.: pw.x -npool 8 -inp  ni.scf.in > ni.scf.out"""
        cmd     = "%s -npool %d -inp  %s > %s" % (simtype, npool, input, output)

        jobid   = self._scheduler.submit(cmd)
#        status  = self.status(jobid)
#
#        print "Simulation started: %s" % status['time_start']
#
#        self._state = status['state']
#        import time
##        print "State: %s, Time started: %s" % (status['state'], status['time_start'])
##
##        status  = self.cancel(jobid)
##        print "State: %s, Time started: %s" % (status['state'], status['time_start'])
#
#        while (status['state'] != 'finished'):
#            print "State: %s, Time: %s" % (status['state'], time.ctime())
#            import time
#            time.sleep(3)
#            status  = self._scheduler.status(jobid)
#            self._state = status['state']
#
#        print "State: %s, Time: %s" % (status['state'], time.ctime())

        return jobid

    def _initScheduler(self):
        servername  = self.settings.get("server", "server-name")
        username    = self.settings.get("user", "username")
        simpath     = self.settings.get("simulation", "remote-path")
        server      = Server(servername, None, username)

        launch = lambda cmd: self._director.csaccessor.execute(cmd, server, simpath, suppressException=True)
        self._scheduler       = Torque(launch, self._director) #launcher) launcher    = self._director.csaccessor.execute


    def status(self, jobid):
        "check status of a job"
        # See for more ideas the old version

        if self._scheduler:
            return self._scheduler.status(jobid)

        return None
    
    def trace(self, jobid):
        pass

    def cancel(self, jobid ):
        "cancel a job"
        status  = self.status(jobid) # Bad!

        if status['state'] != 'running':
            return status

        self._scheduler.delete( jobid )
        
        return self.status(jobid)

if __name__ == "__main__":
    s   = Scheduler(None)
    s.schedule()



# ****** DEAD CODE *****************

#    def schedule(self, job, director ):
#        # copy local job directory to server
#        server          = director.server
#        server_jobpath  = director.dds.abspath(job, server=server)
#
#        # the server
#        server = job.server.dereference(director.clerk.db)
#
#        # the scheduler
#        scheduler = schedulerfactory( server )
#        launch = lambda cmd: director.csaccessor.execute(
#            cmd, server, server_jobpath, suppressException=True)
#        scheduler = scheduler(launch, prefix = 'source ~/.vnf' )
#
#        # submit job through scheduler
#        walltime = job.walltime
#        from pyre.units.time import hour
#        walltime = walltime*hour
#        id1 = scheduler.submit( 'cd %s && sh run.sh' % server_jobpath, walltime=walltime )
#
#        # write id to the remote directory
#        director.csaccessor.execute('echo "%s" > jobid' % id1, server, server_jobpath)
#
#        # update job db record
#        job.id_incomputingserver = id1
#        job.state = 'submitted'
#        import time
#        job.time_start = time.ctime()
#        director.clerk.updateRecordWithID(job)
#
#        return
#
#
#    def check(self, job, director ):
#        "check status of a job"
#
#        if job.state in ['finished', 'failed', 'terminated', 'cancelled']:
#            return job
#
#        oldstate = job.state
#
#        #scheduler
#        server = director.clerk.dereference(job.server)
#        scheduler = schedulerfactory( server )
#
#        #remote job path
#        server_jobpath = director.dds.abspath(job, server=server)
#
#        #
#        launch = lambda cmd: director.csaccessor.execute(
#            cmd, server, server_jobpath, suppressException=True)
#
#        scheduler = scheduler(launch, prefix = 'source ~/.vnf' )
#
#        jobstatus = scheduler.status( job.id_incomputingserver )
#
#        for k,v in jobstatus.iteritems():
#            setattr(job, k, v)
#            continue
#
#        director.clerk.updateRecordWithID( job )
#
#        newstate = job.state
#
#        if oldstate != newstate:
#            # alert user
#            user = director.clerk.getUser(job.creator)
#
#            from vnf.components.misc import announce
#            announce(director, 'job-state-changed', job, user)
#
#        return job
#
#
#    def cancel(self, job, director ):
#        "cancel a job"
#
#        if job.state not in ['running']:
#            return job
#
#        oldstate = job.state
#
#        #scheduler
#        server = director.clerk.dereference(job.server)
#        scheduler = schedulerfactory( server )
#
#        #remote job path
#        server_jobpath = director.dds.abspath(job, server=server)
#
#        #
#        launch = lambda cmd: director.csaccessor.execute(
#            cmd, server, server_jobpath, suppressException=True)
#
#        scheduler = scheduler(launch, prefix = 'source ~/.vnf' )
#
#        scheduler.delete( job.id_incomputingserver )
#
#        job.state = 'cancelled'
#        director.clerk.updateRecordWithID( job )
#
#        newstate = job.state
#
#        if oldstate != newstate:
#            # alert user
#            user = director.clerk.getUser(job.creator)
#
#            from vnf.components.misc import announce
#            announce(director, 'job-state-changed', job, user)
#
#        return job

#
#
## Creates scheduler (torque)
#def schedulerfactory( server ):
#    'obtain scheduler factory'
#    #right now, scheduler info is saved in db record of the server
#    scheduler = server.scheduler
#    if scheduler in [ None, '', 'None' ]:
#        raise RuntimeError, "scheduler not specified"
#
#    from vnf.clusterscheduler import scheduler as factory
#    try: scheduler = factory( scheduler )
#    except: raise NotImplementedError, 'scheduler %r' % scheduler
#    return scheduler

