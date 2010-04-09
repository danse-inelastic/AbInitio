
import os
from twisted.python import log
from twisted.application.service import MultiService
from twisted.application.internet import TCPClient, TCPServer
#from twisted.application import strports

JOB_MESSAGE = " -> JobTracker"

from twisted.protocols.basic import LineReceiver
from twisted.internet.protocol import ClientFactory, ServerFactory
from twisted.internet import reactor

from cassandra.applications.Daemon import Daemon

class JobClient(LineReceiver):

    def connectionMade(self):
        log.msg("JobClient: got new client!")
        self.factory.client   = self    # Save the client


    def lineReceived(self, line):
        log.msg("JobClient: received", repr(line))
        proxy   = self.factory.proxy
        proxy.sendLine(line + JOB_MESSAGE)


class JobTracker(LineReceiver):

    def connectionMade(self):
        print "JobTracker: got new client!"
        self.clientFactory          = self.master.clientFactory()
        self.clientFactory.proxy    = self

#        self.clientFactory     = JobClientFactory()
#        self.clientFactory.client   = None
#        self.clientFactory.proxy    = self      # Important line!
#        reactor.connectTCP('localhost', TT_PORT, self.clientFactory)


    def connectionLost(self, reason):
        print "JobTracker: lost a client!"
        self.factory.client = None


    def lineReceived(self, line):
        print "JobTracker: received", repr(line)
        self.clientFactory.sendLine(line + JOB_MESSAGE)

#        message     = line + JOB_MESSAGE
#        client      = self.clientFactory.client
#        client.sendLine(message)



class TaskMaster(Daemon):   #MultiService):
#    debug = 0
#
    def __init__(self, basedir, configFileName = "taskmaster.cfg"):
        Daemon.__init__(self, basedir, configFileName)
#        MultiService.__init__(self)
#        self._basedir                   = basedir
#        self._configFileName            = configFileName
        self._clientFactory             = ClientFactory()  # Rename? # Later use: pb.PBServerFactory(p)
        self._clientFactory.protocol    = JobClient
        self._serverFactory             = ServerFactory()
        self._serverFactory.protocol    = JobTracker
        self._serverFactory.master      = self


    def startService(self):
        "Starts service"
        Daemon.startService(self)

#        self._loadConfig()
        self._setParams()
        self._listen()
#        self._connect()


#    def stopService(self):
#        "Stops service"
#        MultiService.stopService(self)


    def clientFactory(self):
        return self._clientFactory


    def serverFactory(self):
        return self._serverFactory


#    def _loadConfig(self):
#        configFile = os.path.join(self._basedir, self._configFileName)
#
#        log.msg("loading configuration from %s" % configFile)
#
#        try:
#            f = open(configFile, "r")
#        except IOError, e:
#            log.msg("unable to open config file '%s'" % configFile)
#            log.err(e)
#            return
#
#        try:
#            self._parseConfig(f)
#        except:
#            log.msg("error during parseConfig")
#            log.err()
#        f.close()
#
#
#    def _parseConfig(self, file):
#        "Parsing configuration file"
#        localDict = {'basedir': os.path.expanduser(self._basedir)}
#        try:
#            exec file in localDict
#        except:
#            log.msg("error while parsing config file")
#            raise
#
#        try:
#            config = localDict['CassandraConfig']
#        except KeyError:
#            log.err("missing config dictionary")

    def _setParams(self):
        self._setConfig()
        if not self._config:
            return

        self._masterPort    = self._config["masterPort"]
        self._workerPort    = self._config["workerPort"]
        self._workerHost    = self._config["workerHost"]


    def _connect(self):
        "Connect to the worker"
        masterClient = TCPClient(self._workerHost, self._workerPort, self._clientFactory)
        masterClient.setServiceParent(self)


    def _listen(self):
        "Set service to listen on tcp port"
        masterService = TCPServer(self._masterPort, self._serverFactory)
        masterService.setServiceParent(self)
        
#        # _masterPort supposed to be a strports specification
#        if type(self._masterPort) is int:
#            self._masterPort = "tcp:%d" % self._masterPort
#
#        if self._masterPort:
#            self.clientPort = strports.service(self._masterPort, self._clientFactory)
#            self.clientPort.setServiceParent(self)
#            log.msg("TaskMaster listening on port %s" % self._masterPort)


# DEAD CODE


#class JobClientFactory(ClientFactory):
#    protocol       = JobClient
#
#    def clientConnectionFailed(self, connector, reason):
#        print 'connection failed:', reason.getErrorMessage()
#        reactor.stop()
#
#    def clientConnectionLost(self, connector, reason):
#        print 'connection lost:', reason.getErrorMessage()
#        reactor.stop()
