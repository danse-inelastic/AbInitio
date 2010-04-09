
#from twisted.python import log
#from twisted.application.service import MultiService
from twisted.application.internet import TCPServer
from twisted.protocols.basic import LineReceiver
from twisted.internet.protocol import ServerFactory



from cassandra.applications.Daemon import Daemon

TASK_MESSAGE = " -> TaskTracker"


"""
Starts as a daemon. Receives a line and sends back the attached task message
"""

class TaskTracker(LineReceiver):
    def connectionMade(self):
        print "TaskTracker: got new client!"
        self.factory.server = self  # Save the server


    def connectionLost(self, reason):
        print "TaskTracker: lost a client!"
        self.factory.server = None


    def lineReceived(self, line):
        print "TaskTracker: received", repr(line)
        self.sendLine(line + TASK_MESSAGE)



class Worker(Daemon):

    def __init__(self, basedir, configFileName = "taskmaster.cfg"):
        Daemon.__init__(self, basedir, configFileName)
        self._serverFactory             = ServerFactory()
        self._serverFactory.protocol    = TaskTracker


    def startService(self):
        "Starts service"
        Daemon.startService(self)

        self._setParams()
        self._listen()


    def _setParams(self):
        self._setConfig()
        if not self._config:
            return

        self._workerPort    = self._config["workerPort"]


    def _listen(self):
        "Set service to listen on tcp port"
        workerService = TCPServer(self._workerPort, self._serverFactory)
        workerService.setServiceParent(self)
