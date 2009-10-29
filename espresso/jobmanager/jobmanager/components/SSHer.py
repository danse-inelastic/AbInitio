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

import os
from jobmanager.utils.spawn import spawn
from jobmanager.utils import LOCALHOST
from pyre.components.Component import Component

"""
Performs ssh connection, executes commands remotely and creates ssh tunnel to remote server
Parameters can be configured from ssher.pml
"""

class SSHer(Component):

    class Inventory(Component.Inventory):
        import pyre.inventory
        known_hosts = pyre.inventory.str( 'known_hosts' )
        private_key = pyre.inventory.str( 'private_key' )
        user        = pyre.inventory.str('user', default = '')


    def __init__(self, *args, **kwds):
        super(SSHer, self).__init__(*args, **kwds)
        return


#    def copy(self, serverA, pathA, serverB, pathB):
#        """copy recursively from serverA to serverB"""
#
#        # Copy from remote server to remote server
#        if not self._isLocalhost(serverA) and not self._isLocalhost(serverB):
#            return self._copyR2R(serverA, pathA, serverB, pathB)
#
#        # Copy from one directory to other directory on local server
#        if self._isLocalhost(serverA) and self._isLocalhost(serverB):
#            import shutil
#            shutil.copy(pathA, pathB)
#
#        # Copy from local server to remote server
#        if self._isLocalhost(serverA):
#            self._copyL2R(pathA, serverB, pathB)
#
#        # Copy from remote server to local server
#        if self._isLocalhost(serverB):
#            if os.path.isdir(pathB):
#                dir = pathB
#                newfilename = None
#            else:
#                dir = os.path.dirname(pathB)
#                newfilename = os.path.basename(pathB)
#
#            self.getFile(serverA, pathA, localdir=dir, newfilename=newfilename)
#
#        return
#
#
#    def pushdir( self, path, server, remotepath ):
#        """
#        Pushes a local directory to remote server:
#            1. Archives the "path" directory to tar file
#            2. Copy the tar file to remote server
#            3. Untars the atchived file on remote server
#
#        E.g.: pushdir('/a/b/c', server1, '/a1/b1')
#            Directory /a/b/c will be copied to server1 and become /a1/b1/c
#        """
#
#        address     = server.address
#        port        = server.port
#        username    = server.username
#        known_hosts = self.inventory.known_hosts
#        private_key = self.inventory.private_key
#
#        """tar -czf - <sourcepath> | ssh <remotehost> "(cd <remotepath>; tar -xzf -)" """
#        pieces = [
#            'tar',
#            '-czf',
#            '-',
#            path,
#            '|',
#            'ssh',
#            "-o 'StrictHostKeyChecking=no'",
#            ]
#
#        if port:
#            pieces.append( '-p %s' % port )
#
#        if known_hosts:
#            pieces.append( "-o 'UserKnownHostsFile=%s'" % known_hosts )
#
#        if private_key:
#            pieces.append( "-i %s" % private_key )
#
#        pieces += [
#            '%s@%s' % (username, address),
#            '"(cd %s; tar -xzf -)"' % remotepath,
#            ]
#
#        cmd = ' '.join(pieces)
#
#        failed, output, error = spawn( cmd )
#        if failed:
#            msg = '%r failed: %s' % ( cmd, error )
#            raise RemoteAccessError, msg
#        return
#
#
#    def getFile( self, server, remotepath, localdir, newfilename=None):
#        """retrieve file from remote server to local path"""
#
#        address     = server.address
#        port        = server.port
#        username    = server.username
#        known_hosts = self.inventory.known_hosts
#        private_key = self.inventory.private_key
#
#        pieces = [
#            'scp',
#            "-o 'StrictHostKeyChecking=no'",
#            ]
#
#        if port:
#            pieces.append( '-P %s' % port )
#
#        if known_hosts:
#            pieces.append( "-o 'UserKnownHostsFile=%s'" % known_hosts )
#
#        if private_key:
#            pieces.append( "-i %s" % private_key )
#
#        localpath = localdir
#        if newfilename:
#            localpath = os.path.join(localdir, newfilename)
#        pieces += [
#            '%s@%s:%s' % (username, address, remotepath),
#            '%s' % localpath,
#            ]
#
#        cmd = ' '.join(pieces)
#
#        #print "execute: %s" % cmd
#
#        failed, output, error = spawn( cmd )
#        if failed:
#            msg = '%r failed: %s' % (
#                cmd, error )
#            raise RemoteAccessError, msg
#
#        remotedir, filename = os.path.split( remotepath )
#
#        return os.path.join( localdir, filename )
#
#    # REFACTOR!!! Same as getFile() except passing '-r' (recursive) key
#    def getDir( self, server, remotepath, localdir ):
#        """retrieve a directory from remote server to local path"""
#
#        address     = server.address
#        port        = server.port
#        username    = server.username
#        known_hosts = self.inventory.known_hosts
#        private_key = self.inventory.private_key
#
#        pieces = [
#            'scp',
#            "-o 'StrictHostKeyChecking=no'",
#            ]
#
#        if port:
#            pieces.append( '-P %s' % port )
#
#        if known_hosts:
#            pieces.append( "-o 'UserKnownHostsFile=%s'" % known_hosts )
#
#        if private_key:
#            pieces.append( "-i %s" % private_key )
#
#        pieces += [
#            '-r',
#            '%s@%s:%s' % (username, address, remotepath),
#            '%s' % localdir,
#            ]
#
#        cmd = ' '.join(pieces)
#        self._info.log( 'execute: %s' % cmd )
#
#        failed, output, error = spawn( cmd )
#        if failed:
#            msg = '%r failed: %s' % (
#                cmd, error )
#            raise RemoteAccessError, msg
#
#        remotedir, filename = os.path.split( remotepath )
#
#        return os.path.join( localdir, filename )


    def execute( self, cmd, server, remotepath, suppressException = False ):
        """Execute command in the given directory of the given server"""

        address     = server.address
        port        = server.port
        username    = server.username
        known_hosts = self.inventory.known_hosts
        private_key = self.inventory.private_key

        rmtcmd = 'cd %s && %s' % (remotepath, cmd)

        pieces = [
            'ssh',
            "-o 'StrictHostKeyChecking=no'",
            ]

        if port:
            pieces.append( '-p %s' % port )

        if known_hosts:
            pieces.append( "-o 'UserKnownHostsFile=%s'" % known_hosts )

        if private_key:
            pieces.append( "-i %s" % private_key )

        pieces += [
            '%s@%s' % (username, address),
            '"%s"' % rmtcmd,
            ]

        cmd = ' '.join(pieces)

        #print "execute: %s" % cmd
        failed, output, error = spawn( cmd )

        if failed and not suppressException:
            msg = '%r failed: %s' % (
                cmd, error )
            raise RemoteAccessError, msg

        return failed, output, error

#    # Copy from remote server to remote server
#    def _copyR2R(self, serverA, pathA, serverB, pathB):
#        'push a remote file to another remote server'
#
#        if serverA == serverB:
#            pieces = [
#                'cp -r',
#                pathA,
#                pathB,
#                ]
#        elif _TunneledRemoteServer(serverA) or _TunneledRemoteServer(serverB):
#            raise NotImplementedError, 'serverA: %s, serverBd: %s' % (serverA, serverB)
#        else:
#            address2 = serverB.address
#            port2 = serverB.port
#            username2 = serverB.username
#
#            pieces = [
#                'scp',
#                ]
#            if port2:
#                pieces.append('-P %s' % port2)
#            pieces.append('-r %s' % pathA)
#            pieces.append('%s@%s:%s' % (username2, address2, pathB))
#
#        cmd = ' '.join(pieces)
#
#        self.execute(cmd, serverA, '')
#
#        return
#
#    # Copy from local server to remote server
#    def _copyL2R( self, path, server, remotepath ):
#        'push a local file to remote server'
#
#        address     = server.address
#        port        = server.port
#        username    = server.username
#        known_hosts = self.inventory.known_hosts
#        private_key = self.inventory.private_key
#
#        pieces = [
#            'scp',
#            "-o 'StrictHostKeyChecking=no'",
#            ]
#
#        if port:
#            pieces.append( '-P %s' % port )
#
#        if known_hosts:
#            pieces.append( "-o 'UserKnownHostsFile=%s'" % known_hosts )
#
#        if private_key:
#            pieces.append( "-i %s" % private_key )
#
#        pieces += [
#            path,
#            '%s@%s:%s' % (username, address, remotepath),
#            ]
#
#        cmd = ' '.join(pieces)
#
#        #self._info.log( 'execute: %s' % cmd )
#
#        failed, output, error = spawn( cmd )
#        if failed:
#            msg = '%r failed: %s' % (
#                cmd, error )
#            raise RemoteAccessError, msg
#
#        return
#
#    # Checks if server is a localhost
#    def _isLocalhost(self, server):
#        address     = server.address
#        port        = server.port
#        return (port in LOCALHOST["port"]) and (address in LOCALHOST["address"])
#
#    # Don't know why we need this
#    def _isTunneledRemoteServer(self, server):
#        address     = server.address
#        port        = server.port
#        if address not in LOCALHOST["address"]:
#            return False
#
#        if port not in LOCALHOST["port"]:
#            return False
#
#        return True
#
#    pass # end of SSHer


# version
__id__ = "$Id$"

# End of file 
