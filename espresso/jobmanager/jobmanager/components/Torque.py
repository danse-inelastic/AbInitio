# -*- Python -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                                   Jiao Lin
#                      California Institute of Technology
#                        (C) 2007  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


"""
This is a wrapper of torque commands.

It is done by firing the commands and then parse the outputs.
It is not a good implementation strategy because output formats
could be changed constantly.
Better way to do this is to use python bindings of torque.
See several candidates at this moment, all without enoughh documentation:

http://www-unix.mcs.anl.gov/openpbs/patches/pbs_python/README.txt
https://subtrac.sara.nl/oss/pbs_python
http://code.google.com/p/py-pbs/

History: Adopted from 'clusterscheduler/torque.py'
"""


STATES = {
    'C': 'finished',
    'R': 'running',
    'Q': 'queued',
    'E': 'exiting', #after having run
    'H': 'onhold',
    'W': 'waiting',
    'S': 'suspend',
    }

def _state( state ):
    r = STATES[state]
    if r:
        return r
    return 'unknown state: %s' % state

# Useful if you want to limit time
def _walltime_str(time):
    seconds = int(time/second)%60
    seconds = '%0*d' % (2, seconds)

    mins = int(time/minute)%60
    mins = '%0*d' % (2, mins)

    hours = int(time/hour)
    hours = str(hours)
    return ':'.join([hours, mins, seconds])



import ConfigParser
from pyre.units.time import hour, minute, second    # Useful if you want to limit your simulation!

class SchedulerDaemonNotStarted(Exception):
    pass


class Torque:
    def __init__(self, launcher, director, prefix = None, outputstr_maxlen = 2048):
        self._director          = director
        self._envs              = []
        self._jobid             = 0
        self.prefix             = prefix
        self.launcher           = launcher          # Usually it is just os.system
        self.outputstr_maxlen   = outputstr_maxlen
        self.outfilename        = 'STDOUT.log'
        self.errfilename        = 'STDERR.log'


        self.server             = ConfigParser.ConfigParser()
        self.server.read(self._director.server)
        self.settings           = ConfigParser.ConfigParser()
        self.settings.read(self._director.settings)
        return

    def submit(self, cmd):  # Unlimited time (walltime=1*hour)
        executable  = self.settings.get("server", "executable")
        params      = self.settings.get("server", "params")
        simpath     = self.settings.get("simulation", "simPath")
        simname     = self.settings.get("simulation", "simName")
        nodes       = int(self.settings.get("server", "numNodes"))
        ppn         = int(self.settings.get("server", "procPerNode"))

        """E.g.: echo 'mpirun  --mca btl openib,sm,self pw.x -npool 8 -inp  ni.scf.in > ni.scf.out' | qsub -d /home/dexity/espresso/Ni  -V -N Ni -l nodes=8:ppn=12 -"""

        str     = self.setEnv()
        str     += self.addModules()
        str     += "echo '%s %s %s' | qsub -d %s -V -N %s -l nodes=%d:ppn=%d" % (executable, params, cmd, simpath, simname, nodes, ppn)

        failed, output, error = self._launch( str )

        if failed:
            if error.find( 'check pbs_server daemon' ) != -1:
                raise SchedulerDaemonNotStarted, "pbs_server"
            msg = "error in executing cmds %s. output: %s, error: %s" % ( cmd, output, error )
            raise RuntimeError, msg

        output      = output.strip()
        self._jobid = output.split(".")[0]

        return self._jobid

    # foxtrot specific
    def setEnv(self):
        return """
for i in /etc/profile.d/*.sh; do
       	if [ -r "$i" ]; then
               	. $i
       	fi
done
    """
    
    # foxtrot specific
    def addModules(self):
        modules    = self.server.get("modules", "modules-espresso")
        s   = 'module add %s\n' % modules
        return s
        
    def delete(self, jobid):
        """Deletes job in torque specified by jobid"""
        cmd     = 'qdel %s' % jobid
        failed, output, error  = self._launch( cmd )
        if failed:
            msg = "error in executing cmds %s. output: %s, error: %s" % (
                cmd, output, error )
            raise RuntimeError, msg
        return


    # REDO!!!   - Looks UGLY
    def status( self, jobid ):
        cmd     = 'qstat -f %s' % jobid
        failed, output, error  = self._launch( cmd )
        if failed:
            if error.find( 'Unknown Job Id' ) != -1:
                return self.statusByTracejob( jobid )
            
            msg = "error in executing cmds %s. output: %s, error: %s" % (
                cmds, output, error )
            raise RuntimeError, msg

        lines = output.split( '\n' )
        lines = lines[1:] # first line removed
        if len(lines) == 0: return self.statusByTracejob( jobid )
        d = {}
        for line in lines:
            try:
                k,v = line.split( '=' )
            except:
                continue
            d[ k.strip() ] = v.strip()
            continue

        state = d['job_state']
        import time
        start_time = d.get('start_time') or time.ctime()

        ret = {
            'remote_outputfilename' : self.outfilename,
            'remote_errorfilename'  : self.errfilename,
            'state'                 : _state( state ),
            'time_start'            : start_time
            }

        if ret['state'] == 'finished':
            output, error = self._readoutputerror(self.outfilename, self.errfilename )

            ret.update(
                { 'exit_code'       : d['exit_status'],
                  'time_completion' : d['mtime'],
                  'output'          : output,
                  'error'           : error
                  } )
            pass

        return ret

    # REDO!!!   - Looks UGLY
    def statusByTracejob( self, jobid ):

        d = {}

        tag = 'Exit_status'
        try:
            words = self._tracejob_search( jobid, tag )
        except self.TracejobFailed:
            # this job must have been terminated for a long time
            return self.unknownTerminatedStatus(jobid)

        status = words[3]
        key, value = status.split( '=' )
        assert key.lower() == 'exit_status'
        d [ 'exit_code' ] =  value

        tag = 'job was terminated'
        words = self._tracejob_search( jobid, tag )
        d[ 'time_completion' ] = ' '.join( words[0:2] )

        output, error = self._readoutputerror(
            self.outfilename, self.errfilename )

        d.update( {
            'output': output,
            'error': error,
            'state': 'terminated',
            } )

        return d

    # REDO!!!   - Looks UGLY
    def unknownTerminatedStatus(self, jobid):
        d = {}
        d['exit_code'] = '999999'
        d['state'] = 'terminated'
        output, error = self._readoutputerror(self.outfilename, self.errfilename )

        d.update( {
            'output': output,
            'error': error,
            } )
        return d

    # REDO!!!   - Looks UGLY
    def _readoutputerror(self, outputfilename, errorfilename ):
        return self._read( outputfilename ), self._read( errorfilename )

    # REDO!!!   - Looks UGLY
    def _read(self, filename):
        'read file in the remote job directory'
        cmd     = 'tail %r' % filename
        failed, output, error = self._launch( cmd )
        if failed:
            msg = "error in executing cmds %s. output: %s, error: %s" % (
                cmds, output, error )
            raise RuntimeError, msg
        
        maxlen = self.outputstr_maxlen
        return output[-maxlen+1:]

    # REDO!!!   - Looks UGLY
    def _tracejob_search(self, jobid, tag):
        jobid   = jobid.strip()
        cmd     =  'tracejob %s | grep %r' % (jobid, tag) 

        failed, output, error = self._launch( cmd )

        if failed:
            msg = "error in executing cmds %s. output: %s, error: %s" % (
                cmds, output, error )
            raise self.TracejobFailed, msg

        # remove trailing \n to make parsing easier
        if output.endswith( '\n' ): output = output[:-1]
        lines = output.split( '\n' )
        words = lines[-1].split( )
        #debug.log( 'words: %s' % words )
        return words

    class TracejobFailed(Exception): pass

    def _launch(self, script):
        return self.launcher( script )

def test():
    import os
    s = Torque( os.system, None )
    print s.submit( 'ls' )
    return


def main():
    test()
    return


if __name__ == '__main__':
    main()





# ****** DEAD CODE *****************


#    def submit( self, cmd, walltime=1*hour ):
#        walltime = _walltime_str(walltime)
#
#        cmds = [ r'echo \"%s\" | qsub -l walltime=%s -o %s -e %s' % (
#            cmd, walltime, self.outfilename, self.errfilename) ]
#        failed, output, error = self._launch( cmds )
#        if failed:
#            if error.find( 'check pbs_server daemon' ) != -1:
#                raise SchedulerDaemonNotStarted, "pbs_server"
#            msg = "error in executing cmds %s. output: %s, error: %s" % (
#                cmds, output, error )
#            raise RuntimeError, msg
#        return output.strip()

#    def status( self, jobid ):
#        cmds = [ 'qstat -f %s' % (jobid,) ]
#        failed, output, error  = self._launch( cmds )
#        if failed:
#            if error.find( 'Unknown Job Id' ) != -1:
#                return self.statusByTracejob( jobid )
#            msg = "error in executing cmds %s. output: %s, error: %s" % (
#                cmds, output, error )
#            raise RuntimeError, msg
#
#        lines = output.split( '\n' )
#        lines = lines[1:] # first line removed
#        if len(lines) == 0: return self.statusByTracejob( jobid )
#        d = {}
#        for line in lines:
#            try:
#                k,v = line.split( '=' )
#            except:
#                continue
#            d[ k.strip() ] = v.strip()
#            continue
#
#        #errorpath = d['Error_Path']
#        #dummy, errorfilename = os.path.split(errorpath)
#        #assert errorfilename == self.errfilename, '%r != %r' % (errorfilename, self.errfilename)
#        errorfilename = self.errfilename
#
#        #outputpath = d['Output_Path']
#        #dummy, outputfilename = os.path.split(outputpath)
#        #assert outputfilename == self.outfilename, '%r != %r' % (outputfilename, self.outfilename)
#        outputfilename = self.outfilename
#
#        state = d['job_state']
#        import time
#        start_time = d.get('start_time') or time.ctime()
#
#        ret = {
#            'remote_outputfilename': outputfilename,
#            'remote_errorfilename': errorfilename,
#            'state': _state( state ),
#            'time_start': start_time,
#            }
#
#        if ret['state'] == 'finished':
#            output, error = self._readoutputerror(
#                outputfilename, errorfilename )
#            ret.update(
#                { 'exit_code': d['exit_status'],
#                  'time_completion': d['mtime'],
#                  'output': output,
#                  'error': error,
#                  } )
#            pass
#
#        return ret
#
#
#    def statusByTracejob( self, jobid ):
#
#        d = {}
#
#        tag = 'Exit_status'
#        try:
#            words = self._tracejob_search( jobid, tag )
#        except self.TracejobFailed:
#            # this job must have been terminated for a long time
#            return self.unknownTerminatedStatus(jobid)
#
#        status = words[3]
#        key, value = status.split( '=' )
#        assert key.lower() == 'exit_status'
#        d [ 'exit_code' ] =  value
#
#        tag = 'job was terminated'
#        words = self._tracejob_search( jobid, tag )
#        d[ 'time_completion' ] = ' '.join( words[0:2] )
#
#        output, error = self._readoutputerror(
#            self.outfilename, self.errfilename )
#
#        d.update( {
#            'output': output,
#            'error': error,
#            'state': 'terminated',
#            } )
#
#        return d
#
#
#    def unknownTerminatedStatus(self, jobid):
#        d = {}
#        d['exit_code'] = '999999'
#        d['state'] = 'terminated'
#        output, error = self._readoutputerror(
#            self.outfilename, self.errfilename )
#
#        d.update( {
#            'output': output,
#            'error': error,
#            } )
#        return d
#
#
#    def _readoutputerror(self, outputfilename, errorfilename ):
#        return self._read( outputfilename ), self._read( errorfilename )
#
#
#    def _read(self, filename):
#        'read file in the remote job directory'
#        cmds = [ 'tail %r' % (filename,) ]
#        failed, output, error = self._launch( cmds )
#        if failed:
#            msg = "error in executing cmds %s. output: %s, error: %s" % (
#                cmds, output, error )
#            raise RuntimeError, msg
#        maxlen = self.outputstr_maxlen
#        return output[-maxlen+1:]
#
#
#    def _tracejob_search(self, jobid, tag):
#        jobid = jobid.strip()
#        cmds = [ 'tracejob %s | grep %r' % (jobid, tag) ]
#
#        failed, output, error = self._launch( cmds )
#
#        if failed:
#            msg = "error in executing cmds %s. output: %s, error: %s" % (
#                cmds, output, error )
#            raise self.TracejobFailed, msg
#
#        # remove trailing \n to make parsing easier
#        if output.endswith( '\n' ): output = output[:-1]
#        lines = output.split( '\n' )
#        words = lines[-1].split( )
#        debug.log( 'words: %s' % words )
#        return words
#

# version
__id__ = "$Id$"

# End of file 
