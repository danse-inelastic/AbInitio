
from twisted.python import log
from twisted.application.service import MultiService
from twisted.application.internet import TCPServer
#from twisted.application import strports

JOB_MESSAGE = " -> JobTracker"

from twisted.protocols.basic import LineReceiver
from twisted.internet.protocol import ServerFactory

class Worker(MultiService):
    debug = 0

    def __init__(self, basedir, configFileName = "worker.cfg"):
        MultiService.__init__(self)
        self._basedir                   = basedir
        self._configFileName            = configFileName

    