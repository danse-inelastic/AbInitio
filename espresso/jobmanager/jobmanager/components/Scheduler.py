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

class Server:
    def __init__(self, address, port, username):
        self.address     = address
        self.port        = port
        self.username    = username

class Scheduler:

    def __init__(self, director):
        self._director   = director
        self.settings           = ConfigParser.ConfigParser()
        self.settings.read(self._director.settings)

    def schedule(self):
        from jobmanager.components.Torque import Torque
        servername  = self.settings.get("server", "serverName")
        username    = self.settings.get("user", "username")
        simpath     = self.settings.get("simulation", "simPath")
        simtype     = self.settings.get("simulation", "simType")
        input       = self.settings.get("simulation", "inputFile")
        output      = self.settings.get("simulation", "outputFile")
        npool       = int(self.settings.get("server", "npool"))

        server  = Server(servername, None, username)
        launch = lambda cmd: self._director.csaccessor.execute(
            cmd, server, simpath, suppressException=True)
        s       = Torque(launch, self._director) #launcher) launcher    = self._director.csaccessor.execute

        """ E.g.: pw.x -npool 8 -inp  ni.scf.in > ni.scf.out"""
        cmd     = "%s -npool %d -inp %s > %s" % (simtype, npool, input, output) # TODO: Consider case when parameters are empty!
        print cmd
        cmd     = "cd %s && sh run.sh" % simpath
        print s.submit(cmd)


        #self._director.csaccessor.execute("bash run.sh", server, remotepath = "/home/dexity/espresso/Ni")
        #self._director.csaccessor.execute(cmd, server, remotepath = "/home/dexity/espresso/Ni")


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

if __name__ == "__main__":
    s   = Scheduler(None)
    s.schedule()

# version
__id__ = "$Id$"

# End of file 
