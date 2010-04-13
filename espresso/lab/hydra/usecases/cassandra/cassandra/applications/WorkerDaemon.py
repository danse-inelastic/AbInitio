

from twisted.python import log
from twisted.application.internet import TCPServer
from twisted.protocols.basic import LineReceiver
from twisted.internet.protocol import ServerFactory

from cassandra.applications.Daemon import Daemon





class WorkerDaemon(Daemon):

    def __init__(self, basedir, configFileName = "taskmaster.cfg"):
        Daemon.__init__(self, basedir, configFileName)
        self._worker             = ServerFactory()
#        self._serverFactory.protocol    = TaskTracker

        log.msg("creating Worker")


    def startService(self):
        "Starts service"
        Daemon.startService(self)
        log.msg("starting Worker service")

        

        self._setParams()
        self._connect()


    def _setParams(self):
        self._setConfig()
        if not self._config:
            return

        self._masterPort    = self._config["masterPort"]
        self._masterHost    = self._config["workerHost"]

        log.msg("settings parameters from %s config" % self._configFileName)


    def _connect(self):
        "Connect to the taskmaster"
        workerClient = TCPClient(self._masterHost, self._masterPort, self._worker)
        workerClient.setServiceParent(self)
        log.msg("Worker connects to TaskMaster")


#    def _listen(self):
#        "Set service to listen on tcp port"
#        workerService = TCPServer(self._workerPort, self._serverFactory)
#        workerService.setServiceParent(self)
#        log.msg("TaskMaster is listening on port: %d" % self._workerPort)


#"""
#Starts as a daemon. Receives a line and sends back the attached task message
#"""
#
#class TaskTracker(LineReceiver):
#    def connectionMade(self):
#        log.msg("TaskTracker: got new client!")
#        self.factory.server = self  # Save the server
#
#
#    def connectionLost(self, reason):
#        log.msg("TaskTracker: lost a client!")
#        self.factory.server = None
#
#
#    def lineReceived(self, line):
#        log.msg("TaskTracker: received", repr(line))
#        self.sendLine(line + TASK_MESSAGE)
