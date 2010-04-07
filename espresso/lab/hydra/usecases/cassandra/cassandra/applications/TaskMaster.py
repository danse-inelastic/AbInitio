
import os
from twisted.python import log
from twisted.application.service import MultiService
from twisted.internet.protocol import ServerFactory
from twisted.application import strports

class TaskMaster(MultiService):
    debug = 0

    def __init__(self, basedir, configFileName = "taskmaster.cfg"):
        MultiService.__init__(self)
        self._basedir           = basedir
        self._configFileName    = configFileName
        self._slavePort         = None  # Should rename?
        self._slaveFactory      = ServerFactory()  # Start with small # Later use: pb.PBServerFactory(p)
        

    def startService(self):
        "Starts service"
        MultiService.startService(self)

        self._loadConfig()
        self._listen()


    def _loadConfig(self):
        configFile = os.path.join(self._basedir, self._configFileName)

        log.msg("creating TaskMaster")
        log.msg("loading configuration from %s" % configFile)
        
        try:
            f = open(configFile, "r")
        except IOError, e:
            log.msg("unable to open config file '%s'" % configFile)
            log.err(e)
            return

        try:
            self._parseConfig(f)
        except:
            log.msg("error during parseConfig")
            log.err()
        f.close()


    def _parseConfig(self, file):
        "Parsing configuration file"
        localDict = {'basedir': os.path.expanduser(self._basedir)}
        try:
            exec file in localDict
        except:
            log.msg("error while parsing config file")
            raise

        try:
            config = localDict['CassandraConfig']
        except KeyError:
            log.err("missing config dictionary")

        self._portNum   = config["taskmasterPort"]
        log.msg("setting port number %s " % self._portNum)  # For fun


    def _listen(self):
        "Set service to listen on tcp port"
        # _portNum supposed to be a strports specification
        if type(self._portNum) is int:
            self._portNum = "tcp:%d" % self._portNum

        if self._portNum:
            self.slavePort = strports.service(self._portNum, self._slaveFactory)
            self.slavePort.setServiceParent(self)
            log.msg("TaskMaster listening on port %s" % self._portNum)