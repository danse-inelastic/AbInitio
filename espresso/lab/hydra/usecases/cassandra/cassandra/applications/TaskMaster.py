
import os
from twisted.python import log
from twisted.application.service import MultiService

class TaskMaster(MultiService):
    debug = 0

    def __init__(self, basedir, configFileName = "taskmaster.cfg"):
        MultiService.__init__(self)
        self._basedir           = basedir
        self._configFileName    = configFileName
        

    def startService(self):
        "Starts service"
        MultiService.startService(self)
        self.loadConfig()


    def loadConfig(self):
        configFile = os.path.join(self._basedir, self._configFileName)

        log.msg("Creating TaskMaster")
        log.msg("loading configuration from %s" % configFile)
        
        try:
            f = open(configFile, "r")
        except IOError, e:
            log.msg("unable to open config file '%s'" % configFile)
            log.msg("leaving old configuration in place")
            log.err(e)
            return

        try:
            self._parseConfig(f)
        except:
            log.msg("error during parseConfig")
            log.err()
            log.msg("The new config file is unusable, so I'll ignore it.")
            log.msg("I will keep using the previous config file instead.")
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

        open("/tmp/taskmaster-parse", "w").write("parseConfig %s " % config["taskmasterPort"])
        